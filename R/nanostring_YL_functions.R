

#' @title
#' pagebreak
#'
#' @description
#' A simple function to create a pagebreak in a rmarkdown report.
#' 
#' @return
#' The correct command for a pagebreak (latex/PDF or html)
#'
#' @author Yohan Lefol
pagebreak <- function() {
  if(knitr::is_latex_output())
    return("\\newpage")
  else
    return('<div style="page-break-before: always;" />')
}


#' @title
#' QC_plus_RUV_processing
#'
#' @description
#' This function reads in RCC (nanostring) files, performs a set of quality controls,
#' and finally processed the data using RUV (Remove Unwanted Variance).
#' Depending on the parameter settings the function will also remove/filter out
#' housekeeping genes which are found to correlate with the condition.
#' 
#' @param guide_vect A vector with the files to read in as values (including the path)
#' and the associated group of each files as the name in the vector.
#' @param ref_name The group name to be used as the reference/control
#' @param remove_failed_hk Boolean indicating if housekeeping genes correlating 
#' with the condition should be removed from the downstream analysis. By default
#' it is set to FALSE.
#' 
#' 
#' RUV_total(raw=raw,pData=QCdata,fData=fData,k=1,exclude = bad_hk)
#' @return
#' A summaryexperiment object with the raw and processed data, it also contains
#' the quality control results.
#'
#' @author Yohan Lefol
QC_plus_RUV_processing<-function(guide_vect,ref_name,remove_failed_HK=FALSE){
  files.RCC<-unname(guide_vect)
  
  #Create raw expression matrix
  raw_expression = as.data.frame(matrix(nrow = nrow(readRcc(files.RCC[1])$Code_Summary),
                                        ncol = length(files.RCC)+2))
  colnames(raw_expression)[1:2] = c('Gene','Class')
  raw_expression[,1:2] = readRcc(files.RCC[1])$Code_Summary[,c(2,1)]
  
  #Initialize the QC data
  QCdata = as.data.frame(matrix(nrow = length(files.RCC),ncol = 11))
  colnames(QCdata) = c('BCAC_ID','SampleID','Owner','Comments','Date','GeneRLF','SystemAPF','imagingQC',
                       'bindingDensityQC','limitOfDetectionQC','positiveLinearityQC')
  
  
  #General QC
  for (i in 1:length(files.RCC)){
    rcc = readRcc(files.RCC[i])
    raw = rcc$Code_Summary
    #Update raw_expression
    raw_expression[,i+2] = as.numeric(raw$Count)
    QCdata[i,2:7] = as.vector(rcc$Sample_Attributes)[1:6]
    QCdata$imagingQC[i] = imagingQC(rcc)
    QCdata$bindingDensityQC[i] = bindingDensityQC(rcc,.05,2.25)
    QCdata$limitOfDetectionQC[i] = limitOfDetectionQC(rcc)
    QCdata$positiveLinearityQC[i] = positiveLinQC(rcc)
  }
  
  
  #Prep housekeeping QC, requires matrix and classData (fData)
  raw = raw_expression[,-c(1:2)]
  fData = raw_expression[,c(1:2)]
  rownames(raw) = fData$Gene
  
  #identify housekeeping (hk) genes
  hkIdx <- fData$Gene[fData$Class == "Housekeeping"]
  
  #Formatting
  rownames(raw) = fData$Gene
  row.names(QCdata)=QCdata$SampleID
  colnames(raw)=QCdata$SampleID
  # rownames(QCdata) = colnames(raw)
  QCdata$Group<-names(guide_vect)
  
  #Extract values for housekeeping genes
  hk_raw = raw[hkIdx,]
  #Empty vector for pvalues
  pval = vector(length = nrow(hk_raw))
  
  #Give a pvalue to see if the housekeeping genes are associated to Groups or not
  #Here we would expect them to be not significant
  for (i in 1:nrow(hk_raw)){
    reg = glm.nb(as.numeric(hk_raw[i,]) ~ as.factor(QCdata$Group))
    pval[i] = coef(summary(reg))[2,4]
  }
  #Give a readable output for a HKs association to groups
  names(pval)<-row.names(hk_raw)
  bad_hk<-names(pval[pval<0.05])
  if(length(bad_hk)!=0){
    if(remove_failed_HK==TRUE){
      message('The following housekeeping genes were found to be correlated to the conditions and have\nbeen removed')
      message(paste(bad_hk,collapse = ', '))
    }else{
      message('The following housekeeping genes were found to be correlated to the conditions')
      message(paste(bad_hk,collapse = ', '))
      message("To remove them, set the 'remove_failed_HK' parameter to TRUE.")
    }

  }
  
  #Fetch the negative controls and extract values
  negative_controls<-row.names(raw)[startsWith(row.names(raw),'NEG')]
  neg_raw = raw[negative_controls,]
  #Calculate the limit of detection value
  lod = colMeans(neg_raw) - apply(neg_raw,2,sd)
  QCdata$LODvalues<-unname(lod)
  #Create vectors establishing how many genes of each samples fall below the lod
  raw$Code.Class<-fData$Class
  num_endogenous_blod = colSums(raw[raw$Code.Class == 'Endogenous',] < lod)[-ncol(raw)]
  percent_endogenous_blod=num_endogenous_blod/nrow(raw[raw$CodeClass=='Endogenous'])[-ncol(raw)]
  
  num_hk_blod = colSums(raw[raw$Code.Class == 'Housekeeping',] < lod)[-ncol(raw)]
  # hk_failling_lod<-row.names(raw[raw$Code.Class == 'Housekeeping',] < lod)
  # if(length(hk_failling_lod)>0){
  #   message('The following genes have LOD failed in some or several samples: ')
  #   message(message(paste(hk_failling_lod,collapse = ', ')))
  # }
  #Establish samples in reference group
  ref_samples<-row.names(QCdata)[QCdata$Group==ref_name]
  
  #Create endogenous over 75% quantile check - rule below
  # A sample should be flagged/marked if the percentage of endogenous genes below 
  # the limit of detection (LOD) is greater than the top quantile of this percentage (75%)
  # in a reference group.
  threshold<-unname(quantile(percent_endogenous_blod[ref_samples])['75%'])
  flag_endogenous_lod<-unname(percent_endogenous_blod>threshold)
  
  QCdata$EndogenousLod<-flag_endogenous_lod
  QCdata$EndogenousLod[row.names(QCdata)%in%ref_samples]<- 'No flag' #No flag if a reference
  QCdata$EndogenousLod[QCdata$EndogenousLod == FALSE]<- 'No flag'
  QCdata$EndogenousLod[QCdata$EndogenousLod == TRUE]<- 'Flag'
  
  
  #Flag if any HK genes fall below the detectable limit
  QCdata$HKLod_values<-num_hk_blod
  QCdata$HKLod[QCdata$HKLod_values == 0]<- 'No flag'
  QCdata$HKLod[QCdata$HKLod_values > 0]<- 'Flag'
  
  
  
  #Some final processing
  raw<-raw[,-ncol(raw)]
  row.names(fData)=fData$Gene
  
  #Run RUV on the data
  if(remove_failed_HK==TRUE){
    RUV_processed<-RUV_total(raw=raw,pData=QCdata,fData=fData,k=1,exclude = bad_hk)
  }else{
    RUV_processed<-RUV_total(raw=raw,pData=QCdata,fData=fData,k=1)
  }
  fData(RUV_processed$set)<-fData
  return(RUV_processed)
  
}

#' @title
#' create summary table
#'
#' @description
#' Function which creates a table which summarises the results for a specific patient.
#' The table is intended to be plotted alongside a plot for that same patient.
#' 
#' @param subset_df A dataframe containing the geoscore of the targeted patient
#' @param normalized_counts A matrix of normalized counts. Gene names and sample names must be included.
#' @param samp_sheet A sample sheet which associates samples names to groups
#' @param target_patient Character indicating which patient the table should be made for
#' @param target_genes The genes which will be shown in the table along with the geoscore
#' 
#' @return
#' A dataframe with the geoscore and normalized gene counts for the defined patient and genes.
#'
#' @author Yohan Lefol
create_summary_table<-function(subset_df,normalized_counts,samp_sheet,target_patient,target_genes){
  normalized_counts<-as.data.frame(normalized_counts)
  
  patient_counts<-normalized_counts[samp_sheet$Lane[startsWith(samp_sheet$Sample,target_patient)]]
  target_genes<-target_genes[target_genes %in% row.names(patient_counts)]
  for_table<-rbind(as.character(round(subset_df$geoscore_ISG[startsWith(subset_df$Sample,target_patient)],digits=2)),
                   patient_counts[target_genes,,drop=FALSE])
  row.names(for_table)=c('geoscore_ISG',target_genes)
  colnames(for_table)=subset_df$Sample[startsWith(subset_df$Sample,target_patient)]
  
  return(for_table)
}

#' @title
#' analysis_with_geneset
#'
#' @description
#' 
#' @param 
#' 
#' @return
#' 
#'
#' @author Yohan Lefol
analysis_with_geneset<-function(main_data,target_genes,sample_data, group_name='',save_loc=''){
  
  #Start by correcting column names in main_data - currently it represents
  #lanes, we want it to be sample IDs
  main_data<-main_data[,sample_data$Lane]#Sort as in the sample sheet
  colnames(main_data)=sample_data$Sample
  
  #Subset genes for the ones present in the data
  genes_not_in_data<-target_genes[! target_genes %in% row.names(main_data)]
  target_genes<-target_genes[ target_genes %in% row.names(main_data)]
  
  sub_dta<-main_data[target_genes,]
  
  #Formatting for boxplot
  df_long <- melt(sub_dta, variable.name = "Sample", value.name = "Expression")
  colnames(df_long)=c('Genes','Sample','Expression')
  patients<-sample_data$Sample[sample_data$Group=='Pasient']
  group_info <- data.frame(Sample = colnames(sub_dta), 
                           Group = ifelse(colnames(sub_dta) %in% patients, "Pasient", "Kontroll"))
  df_long <- merge(df_long, group_info, by = "Sample")
  
  #Add translated column name instead of renaming cause lazy
  df_long$Gruppe<-df_long$Group
  
  boxplot<-ggplot(df_long, aes(x = Sample, y = Expression, fill = Gruppe)) +
    geom_boxplot() +
    scale_fill_manual(values = c("Kontroll" = "#00bfc4", "Pasient" = "#f8766d")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
    labs(title = paste0("Gene Expression per Sample for group: ",group_name), x = "", y = "genuttrykk")
  
  #save the plot
  ggsave(filename = paste0(save_loc,'summary_boxplot.png'),plot = boxplot, dpi=600)
  
  #Calculate geoscore
  nth_root <- function(x, n) {
    x^(1/n)
  }
  geo_score<-nth_root(apply(sub_dta, 2, prod),length(target_genes))
  
  #Create group vector
  groups<-sample_data$Group
  names(groups)<-sample_data$Sample

  sub_dta<-rbind(unname(geo_score),sub_dta)
  row.names(sub_dta)[1]<-'geoscore'
  
  #Add NAs for genes not in the data  
  NA_df<-data.frame(matrix(NA, nrow = length(genes_not_in_data), ncol = ncol(sub_dta)))
  row.names(NA_df)=genes_not_in_data
  colnames(NA_df)=colnames(sub_dta)
  
  #Put everything together
  results_df<-rbind(sub_dta,NA_df)
  
  missing_genes<-row.names(results_df)[is.na(results_df[1])]
  message('The following genes were not found in the nanostring data:')
  message(paste(missing_genes,collapse = ', '))
  
  return(list('boxplot'=boxplot,'results_df'=results_df))
  
}

#' @title
#' 
#'
#' @description
#' 
#' @param 
#' 
#' @return
#' 
#'
#' @author Yohan Lefol
plots_and_table_wrapper<-function(res_list,samp_sheet,norm_counts,target_genes,save_loc){
  
  #Have one entry per sample, including multi-sampled patients
  all_patients<-samp_sheet$Sample[samp_sheet$Group=='Pasient']
  parts <- strsplit(all_patients, "-")
  parts_trimmed <- lapply(parts, function(x) x[1:min(2, length(x))])
  all_patients <- unique(sapply(parts_trimmed, paste, collapse = "-"))
  
  #Get and format the data
  geoscore_ISG<-t(res_list$results_df['geoscore',])
  geoscore_ISG<-as.data.frame(cbind(geoscore_ISG,samp_sheet[match(row.names(geoscore_ISG), samp_sheet$Sample),'Group']))
  colnames(geoscore_ISG)=c('geoscore_ISG','Group')
  geoscore_ISG$Sample<-row.names(geoscore_ISG)
  geoscore_ISG$geoscore_ISG<-as.numeric(geoscore_ISG$geoscore_ISG)
  
  #Add timepoint
  tp_vect<-c()
  for(i in geoscore_ISG$Sample){
    timepoint<-unlist(strsplit(i,'-'))[c(F,F,T)]
    if(is.na(timepoint)==T){
      timepoint=1
    }
    tp_vect<-c(tp_vect,timepoint)
  }
  geoscore_ISG$timepoint<-tp_vect
  
  #Create quantiles
  control_quantiles<-quantile(geoscore_ISG$geoscore_ISG[geoscore_ISG$Group=='Kontroll'])
  message('Control quantiles: ')
  for(quant in names(control_quantiles)){
    message(paste0(quant,': ',round(control_quantiles[quant],digits=2)))
  }
  #Create plots for all patients
  for(target_patient in all_patients){
    sub_df<-rbind(geoscore_ISG[geoscore_ISG$Group=='Kontroll',], geoscore_ISG[startsWith(geoscore_ISG$Sample,target_patient),])
  
    pat_plot<-create_quantile_plot(plot_dta=sub_df,quantiles=control_quantiles)

    found_table<-create_summary_table(sub_df,norm_counts,samp_sheet,target_patient,target_genes)
    
    # Table theme: smaller font
    theme <- ttheme_default(core = list(fg_params=list(cex=0.8)), 
                            colhead = list(fg_params=list(cex=0.7)),
                            rowhead = list(fg_params=list(cex=0.8)))
    
    # Bind the table and plot together, width alternates based on the number of samples
    tg <- tableGrob(found_table, theme = theme)
    if(ncol(found_table)>2){
      width_vect<-c(c(2,2))
    }else{
      width_vect<-c(c(2,1))
    }
    #Store the plot for saving
    combined_plot<-arrangeGrob(pat_plot, tg, ncol = 2, widths = width_vect)
    #Print the plot for the report
    grid.arrange(pat_plot, tg, ncol = 2, widths = width_vect)
    #save plot
    ggsave(filename = paste0(save_loc,target_patient,'.png'),plot = combined_plot,dpi=600)
    
    #Page break for markdown
    if(target_patient != all_patients[length(all_patients)]){#If it isn't the last patient
      cat('\\pagebreak')
    }
  }
}


#' @title
#' 
#'
#' @description
#' 
#' @param 
#' 
#' @return
#' 
#'
#' @author Yohan Lefol
create_quantile_plot<-function(plot_dta,quantiles){
  #Find yaxis upper limit
  max_score=round_any(max(plot_dta[,'geoscore_ISG']), 100, f=ceiling)
  if(max_score-max(plot_dta[,'geoscore_ISG'])<25){
    max_score=max_score+100
  }
  
  patient_plot<-ggplot(plot_dta, aes(x = timepoint, y = geoscore_ISG, color=Group, fill = Group)) +
    # Set the quantile area for the controls
    geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=quantiles[1],ymax=quantiles[2]),
              fill="#9fe3e5",color=NA)+
    geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=quantiles[2],ymax=quantiles[4]),
              fill="#00bfc4",color=NA)+
    geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=quantiles[4],ymax=quantiles[5]),
              fill="#9fe3e5",color=NA)+

    
    #Add the quantile delimination lines
    geom_hline(yintercept=quantiles[1], linetype="dashed", color = "black",linewidth = 0.1)+
    geom_hline(yintercept=quantiles[2], linetype="dashed", color = "black",linewidth = 0.1)+
    geom_hline(yintercept=quantiles[3], linetype="solid", color = "black",linewidth = 0.25)+
    geom_hline(yintercept=quantiles[4], linetype="dashed", color = "black",linewidth = 0.1)+
    geom_hline(yintercept=quantiles[5], linetype="dashed", color = "black",linewidth = 0.1)+
    
    annotate(geom = 'text', label = '0-25', x = Inf, y = quantiles[1]+(quantiles[2]-quantiles[1])/2, size=2, hjust=1, vjust=0.5)+
    annotate(geom = 'text', label = '25-50', x = Inf, y = quantiles[2]+(quantiles[3]-quantiles[2])/2, size=2, hjust=1, vjust=0.5)+
    annotate(geom = 'text', label = '50-75', x = Inf, y = quantiles[3]+(quantiles[4]-quantiles[3])/2, size=2, hjust=1, vjust=0.5)+
    annotate(geom = 'text', label = '75-100', x = Inf, y = quantiles[4]+(quantiles[5]-quantiles[4])/2, size=2, hjust=1, vjust=0.5)

  
  #Draw the line for the patient samples (if there are several sample for the patient)
  if(length(plot_dta$Sample[plot_dta$Group=='Pasient'])>1){
    patient_plot<-patient_plot+geom_line(data = subset(plot_dta, Group == "Pasient"),aes(group = 1, color = Group))
  }
  
  
  patient_plot <- patient_plot+      
    #Add the text that will indicate which dots belong to which patients
    geom_text_repel(data = subset(plot_dta, Group == "Pasient"), aes(label = Sample), 
                    size = 5, 
                    nudge_y = 5, 
                    segment.color = "gray50",
                    color='black') +
    #Add the data points for the patients
    geom_point(data = subset(plot_dta, Group == "Pasient"),shape=21,colour = "black",size = 3)+
    
    #Manage the colors, necessary to do both in order to have a black border around the dots.
    scale_color_manual(values = c("Kontroll" = "#00bfc4", "Pasient" = "#f8766d"),guide="none") +
    scale_fill_manual(values = c("Kontroll" = "#00bfc4", "Pasient" = "#f8766d"),guide="none")+
    
    #Some quality of life/aesthetic improvements
    scale_y_continuous(limits = c(0, max_score),expand = c(0,0))+
    
    #Labels
    labs(title = "",x = "", y = "geometrisk gjennomsnitt") +
    theme_minimal()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank())
  
  
  return(patient_plot)
}




# The functions below were extracted from the NanoStringQCPro package. These functions were
# required for our analysis. Since the NanoStringQCPro was no longer maintained we extracted
# the functions and stored them here. The credit for the two functions below are for the creators
# of the NanoStringQCPro. The functions and documentation below are an exact copy, no alterations
# have been made.
################################# NanoStringQCPro

##' @title          Read an .RCC file
##' @description    Parse an .RCC file into a list with each part of the file (Header,
##'                 Sample_Attributes, Lane_Attributes, Code_Summary, etc) stored as a vector
##'                 or data frame.
##'
##' @param  rcc                 Path to the .RCC file.
##' @param  removeSpikeInLabels Logical. If TRUE (the default), RNA \dQuote{spike-in} input labels (if any)
##'                             in the GeneName for positive and negative control probes will be removed.
##'
##' @return
##' A list where each element holds the contents of one part of the .RCC file
##' (Header, Sample_Attributes, Lane_Attributes, Code_Summary, etc) as a
##' vector or data frame.
##'
##' @export
##'
##' @examples
##' rcc <- system.file("extdata", "RCC", "20140604_C1-unstim_C1-unstim_01.RCC", package="NanoStringQCPro")
##' rcc.ls <- readRcc(rcc)
##'
##' @author Robert Ziman
##'
readRcc <- function(rcc, removeSpikeInLabels=TRUE)
{
  suppressWarnings(
    lines <- readLines(rcc)
  )
  lines.df <-
    data.frame(stringsAsFactors=FALSE,
               line   = lines
               ,tag    = rep("", length(lines))
    )
  
  tags <- c("Header", "Sample_Attributes", "Lane_Attributes", "Code_Summary", "Messages")
  for(i in 1:length(tags))
  {
    tag_start <- which(lines.df$line == paste0("<", tags[i], ">"))
    tag_end <- which(lines.df$line == paste0("</", tags[i], ">"))
    lines.df$tag[ tag_start:tag_end ] <- tags[i]
  }
  
  unrecognized_lines <- lines.df$line[ lines.df$tag == "" ]
  unrecognized_tags <- grep(value=TRUE, "^<", unrecognized_lines)
  unrecognized_tags <- grep(invert=TRUE, value=TRUE, "^</", unrecognized_tags)
  if (length(unrecognized_tags) > 0) {
    warning(paste0("Unrecognized tags: ", paste(collapse=" ", unrecognized_tags)))
  }
  
  lines.df <- lines.df[ (lines.df$tag != ""), ]
  lines.df <- lines.df[ !grepl("<", lines.df$line), ]
  lines.df <- lines.df[ !grepl("</", lines.df$line), ]
  
  ID_lnum <- grep("^ID,", lines.df$line)
  lines.df$line[ID_lnum] <- paste0( sub("_Attributes", "", lines.df$tag[ID_lnum]), lines.df$line[ID_lnum] )
  
  rcc.ls <- list()
  
  rcc.ls$File_path <- rcc
  
  for (attr_tag in c("Header", "Sample_Attributes", "Lane_Attributes")) {
    attr_lnum <- which(lines.df$tag == attr_tag)
    attr.strsplit <- strsplit(lines.df$line[attr_lnum], ",")
    attr.v <- vapply(attr.strsplit, function(x) {if (!is.na(x[2])) {x[2]} else {""}}, character(1))
    names(attr.v) <- vapply(attr.strsplit, function(x) {x[1]}, character(1))
    rcc.ls[[attr_tag]] <- attr.v
  }
  
  Code_Summary_lnum.all <- which(lines.df$tag == "Code_Summary")
  Code_Summary_lnum.header <- Code_Summary_lnum.all[1]
  Code_Summary_lnum.body <- Code_Summary_lnum.all[2:length(Code_Summary_lnum.all)]
  
  if (lines.df$line[Code_Summary_lnum.header] != "CodeClass,Name,Accession,Count") {
    stop(
      paste0("Unrecognized header line for Code_Summary section of \"", rcc.ls, "\"",
             " (expected \"CodeClass,Name,Accession,Count\" but got \"", lines.df$line[Code_Summary_lnum.header], "\")")
    )
  }
  
  Code_Summary_fields <- strsplit(lines.df$line[Code_Summary_lnum.body], ",")
  rcc.ls$Code_Summary <- data.frame(stringsAsFactors=FALSE
                                    ,CodeClass   = vapply(Code_Summary_fields, function(x) {x[1]}, character(1))
                                    ,Name        = vapply(Code_Summary_fields, function(x) {x[2]}, character(1))
                                    ,Accession   = vapply(Code_Summary_fields, function(x) {x[3]}, character(1))
                                    ,Count       = vapply(Code_Summary_fields, function(x) {x[4]}, character(1))
  )
  
  if(removeSpikeInLabels)
  {
    spikein <- getSpikeInInput(CodeClass=rcc.ls$Code_Summary$CodeClass, GeneName=rcc.ls$Code_Summary$Name)
    rcc.ls$Code_Summary$Name <- spikein$GeneName
  }
  
  return(rcc.ls)
}




##' @title
##' getSpikeInInput
##'
##' @description
##' Gets the RNA \dQuote{spike-in} input levels for positive and negative control probes from the label in their GeneName.
##' Note that this is a helper function for readRlf() and elsewhere and is not intended for external use.
##'
##' @param CodeClass        Character vector with code classes for each probe.
##' @param GeneName         Character vector with gene names/symbols for each probe.
##' @param nonCtrlProbeVal  Value to assign as the spike-in input for the non-control probes.
##'
##' @return
##' A data frame with the input CodeClass and GeneName but where the latter has been split into
##' two columns: one showing the GeneName for each probe with spike-in input labels removed
##' -- and another with the spike-in input levels.
##'
##' @author Robert Ziman
##
getSpikeInInput <- function(CodeClass, GeneName, nonCtrlProbeVal=NA)
{
  stopifnot(!is.null(CodeClass))
  stopifnot(!is.null(GeneName))
  stopifnot(length(CodeClass) == length(GeneName))
  
  GeneName.split <- strsplit(GeneName, '\\(')
  GeneName.without_parentheses <- vapply(GeneName.split, function(x) {x[1]}, character(1))
  SpikeInInput.char <- sub('\\)$', '', vapply(GeneName.split, function(x) {x[2]}, character(1)))
  
  ctrls <- (CodeClass %in% c("Positive", "Negative"))
  
  GeneName.adjusted <- GeneName
  GeneName.adjusted[ ctrls ] <- GeneName.without_parentheses[ ctrls ]     # Only strip the labels from the control probes (in some cases there will be GeneNames with parentheses!)
  
  SpikeInInput <- rep(nonCtrlProbeVal, length(GeneName))
  SpikeInInput[ ctrls ] <- as.numeric(SpikeInInput.char[ ctrls ])
  
  return(
    data.frame(stringsAsFactors=FALSE,
               CodeClass   = CodeClass,
               GeneName    = GeneName.adjusted,
               SpikeInInput = SpikeInInput
    )
  )
}





