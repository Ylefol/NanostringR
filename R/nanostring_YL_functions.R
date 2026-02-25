

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
#' @export
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
#' @param guide_vect A named character vector where values are RCC file paths
#' (including the path) and names are the group label for each file.
#' @param ref_name The group name to be used as the reference/control.
#' @param remove_failed_HK Boolean indicating if housekeeping genes correlating
#' with the condition should be removed from the downstream analysis. By default
#' it is set to FALSE.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{set}: A SeqExpressionSet with raw counts and RUV-estimated
#'     unwanted variance factors in the phenotype data.
#'   \item \code{vsd}: A DESeqTransform object with variance-stabilized,
#'     batch-corrected expression values. Use \code{assay(vsd)} to extract
#'     the normalized count matrix.
#' }
#'
#' @author Yohan Lefol
#' @export
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
#' QC_plus_nSolver_processing
#'
#' @description
#' This function reads in RCC (nanostring) files, performs a set of quality controls,
#' and normalizes the data using nSolver-equivalent normalization via NanoStringNorm.
#' Normalization applies positive control scaling (geometric mean) followed by
#' housekeeping gene scaling (geometric mean), matching the default nSolver workflow.
#' Depending on the parameter settings the function will also remove/filter out
#' housekeeping genes which are found to correlate with the condition.
#'
#' @param guide_vect A vector with the files to read in as values (including the path)
#' and the associated group of each file as the name in the vector.
#' @param ref_name The group name to be used as the reference/control.
#' @param remove_failed_HK Boolean indicating if housekeeping genes correlating
#' with the condition should be removed from the downstream analysis. By default
#' it is set to FALSE.
#' @param background_correction Boolean indicating whether to apply background
#' correction using mean + 2SD of negative controls, as per nSolver defaults.
#' Default is TRUE.
#' @param controls_only_norm Boolean. When TRUE, normalization factors (positive
#' control and housekeeping scaling) are computed using only the reference/control
#' samples and then applied to all samples. This ensures patient samples do not
#' influence each other's normalization and that each patient is always scaled
#' against the same fixed control reference, regardless of how many patient
#' samples are processed together. Default is FALSE (all samples contribute to
#' the global reference, matching standard NanoStringNorm behaviour).
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{normalized_counts}: A gene x sample matrix of nSolver-normalized counts
#'   \item \code{QCdata}: A data frame containing QC metrics for each sample
#' }
#'
#' @author Yohan Lefol
#' @export
QC_plus_nSolver_processing <- function(guide_vect, ref_name, remove_failed_HK = FALSE,
                                        background_correction = TRUE,
                                        controls_only_norm = FALSE) {
  files.RCC <- unname(guide_vect)

  # Create raw expression matrix
  first_code_summary <- readRcc(files.RCC[1])$Code_Summary
  raw_expression = as.data.frame(matrix(nrow = nrow(first_code_summary),
                                        ncol = length(files.RCC) + 3))
  colnames(raw_expression)[1:3] = c('Gene', 'Class', 'Accession')
  raw_expression[, 1] = first_code_summary$Name
  raw_expression[, 2] = first_code_summary$CodeClass
  raw_expression[, 3] = first_code_summary$Accession

  # Initialize the QC data
  QCdata = as.data.frame(matrix(nrow = length(files.RCC), ncol = 11))
  colnames(QCdata) = c('BCAC_ID', 'SampleID', 'Owner', 'Comments', 'Date', 'GeneRLF', 'SystemAPF',
                       'imagingQC', 'bindingDensityQC', 'limitOfDetectionQC', 'positiveLinearityQC')

  # General QC
  for (i in 1:length(files.RCC)) {
    rcc = readRcc(files.RCC[i])
    raw = rcc$Code_Summary
    raw_expression[, i + 3] = as.numeric(raw$Count)
    QCdata[i, 2:7] = as.vector(rcc$Sample_Attributes)[1:6]
    QCdata$imagingQC[i] = imagingQC(rcc)
    QCdata$bindingDensityQC[i] = bindingDensityQC(rcc, .05, 2.25)
    QCdata$limitOfDetectionQC[i] = limitOfDetectionQC(rcc)
    QCdata$positiveLinearityQC[i] = positiveLinQC(rcc)
  }

  # Prep housekeeping QC
  raw = raw_expression[, -c(1:3)]
  fData = raw_expression[, c(1:3)]
  rownames(raw) = fData$Gene

  # Identify housekeeping genes
  hkIdx <- fData$Gene[fData$Class == "Housekeeping"]

  # Formatting
  row.names(QCdata) = QCdata$SampleID
  colnames(raw) = QCdata$SampleID
  QCdata$Group <- names(guide_vect)

  # Test if housekeeping genes are associated with Groups
  hk_raw = raw[hkIdx, ]
  pval = vector(length = nrow(hk_raw))
  for (i in 1:nrow(hk_raw)) {
    reg = glm.nb(as.numeric(hk_raw[i, ]) ~ as.factor(QCdata$Group))
    pval[i] = coef(summary(reg))[2, 4]
  }
  names(pval) <- row.names(hk_raw)
  bad_hk <- names(pval[pval < 0.05])
  if (length(bad_hk) != 0) {
    if (remove_failed_HK == TRUE) {
      message('The following housekeeping genes were found to be correlated to the conditions and have\nbeen removed')
      message(paste(bad_hk, collapse = ', '))
    } else {
      message('The following housekeeping genes were found to be correlated to the conditions')
      message(paste(bad_hk, collapse = ', '))
      message("To remove them, set the 'remove_failed_HK' parameter to TRUE.")
    }
  }

  # Fetch negative controls and calculate LOD
  negative_controls <- row.names(raw)[startsWith(row.names(raw), 'NEG')]
  neg_raw = raw[negative_controls, ]
  lod = colMeans(neg_raw) - apply(neg_raw, 2, sd)
  QCdata$LODvalues <- unname(lod)

  # Count genes below LOD per sample
  raw$Code.Class <- fData$Class
  num_endogenous_blod = colSums(raw[raw$Code.Class == 'Endogenous', ] < lod)[-ncol(raw)]
  n_endogenous = nrow(raw[raw$Code.Class == 'Endogenous', ])
  percent_endogenous_blod = num_endogenous_blod / n_endogenous

  num_hk_blod = colSums(raw[raw$Code.Class == 'Housekeeping', ] < lod)[-ncol(raw)]

  # Flag samples where endogenous LOD percentage exceeds 75th quantile of reference group
  ref_samples <- row.names(QCdata)[QCdata$Group == ref_name]
  threshold <- unname(quantile(percent_endogenous_blod[ref_samples])['75%'])
  flag_endogenous_lod <- unname(percent_endogenous_blod > threshold)

  QCdata$EndogenousLod <- flag_endogenous_lod
  QCdata$EndogenousLod[row.names(QCdata) %in% ref_samples] <- 'No flag'
  QCdata$EndogenousLod[QCdata$EndogenousLod == FALSE] <- 'No flag'
  QCdata$EndogenousLod[QCdata$EndogenousLod == TRUE] <- 'Flag'

  QCdata$HKLod_values <- num_hk_blod
  QCdata$HKLod[QCdata$HKLod_values == 0] <- 'No flag'
  QCdata$HKLod[QCdata$HKLod_values > 0] <- 'Flag'

  # Remove Code.Class column added for LOD
  raw <- raw[, -ncol(raw)]

  if (controls_only_norm) {
    # Manual nSolver normalization anchored to control samples only.
    # Normalization references (positive control and housekeeping geometric means)
    # are computed exclusively from ref_samples, then applied to all samples.
    # This ensures patient samples never influence each other's normalization.
    counts_mat <- as.matrix(raw)

    # Step 1: Background correction — threshold at mean + 2SD of negatives per sample
    if (background_correction) {
      neg_rows <- rownames(counts_mat)[startsWith(rownames(counts_mat), 'NEG')]
      for (s in colnames(counts_mat)) {
        bg <- mean(counts_mat[neg_rows, s]) + 2 * sd(counts_mat[neg_rows, s])
        counts_mat[counts_mat[, s] < bg, s] <- bg
      }
    }

    # Step 2: Positive control normalization — reference from controls only
    pos_rows <- fData$Gene[fData$Class == 'Positive']
    pos_geo_means <- apply(counts_mat[pos_rows, , drop = FALSE], 2,
                           function(col) exp(mean(log(pmax(col, 0.5)))))
    ref_pos <- mean(pos_geo_means[ref_samples])
    counts_mat <- sweep(counts_mat, 2, ref_pos / pos_geo_means, FUN = '*')

    # Step 3: Housekeeping normalization — reference from controls only
    hk_rows <- fData$Gene[fData$Class == 'Housekeeping']
    if (remove_failed_HK && length(bad_hk) > 0) {
      hk_rows <- hk_rows[!hk_rows %in% bad_hk]
    }
    hk_geo_means <- apply(counts_mat[hk_rows, , drop = FALSE], 2,
                          function(col) exp(mean(log(pmax(col, 0.5)))))
    ref_hk <- mean(hk_geo_means[ref_samples])
    counts_mat <- sweep(counts_mat, 2, ref_hk / hk_geo_means, FUN = '*')

    keep_rows <- fData$Gene[fData$Class %in% c('Endogenous', 'Housekeeping')]
    norm_counts <- counts_mat[keep_rows, ]

  } else {
    # Build NanoStringNorm input (requires Code.Class, Name, Accession + count columns)
    nsn_input = data.frame(
      Code.Class = fData$Class,
      Name = fData$Gene,
      Accession = fData$Accession,
      raw,
      check.names = FALSE
    )

    # If removing failed HK genes, reclassify them so NanoStringNorm excludes them
    # from the housekeeping normalization step
    if (remove_failed_HK == TRUE && length(bad_hk) > 0) {
      nsn_input$Code.Class[nsn_input$Name %in% bad_hk] <- 'Endogenous'
    }

    # Run NanoStringNorm with nSolver-equivalent settings:
    # geo.mean of positive controls -> background correction -> geo.mean of housekeeping genes
    bg_method <- if (background_correction) 'mean.2sd' else 'none'
    nsn_result = NanoStringNorm(
      x = nsn_input,
      CodeCount = 'geo.mean',
      Background = bg_method,
      SampleContent = 'housekeeping.geo.mean',
      OtherNorm = 'none',
      verbose = FALSE
    )

    norm_data = nsn_result$normalized.data
    keep_classes = c('Endogenous', 'Housekeeping')
    norm_counts = as.matrix(norm_data[norm_data$Code.Class %in% keep_classes,
                                       !colnames(norm_data) %in% c('Code.Class', 'Name', 'Accession')])
    rownames(norm_counts) = norm_data$Name[norm_data$Code.Class %in% keep_classes]
  }

  return(list(normalized_counts = norm_counts, QCdata = QCdata))
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
#' @export
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
#' Computes a geometric mean expression score (geoscore) for a defined gene set
#' and generates a per-sample boxplot comparing patient and control groups.
#' Genes absent from the data are reported and included as NA rows in the output.
#'
#' @param main_data A gene x sample matrix of normalized counts with lane IDs as
#' column names.
#' @param target_genes Character vector of gene names defining the gene set.
#' @param sample_data A sample sheet data frame with columns Sample, Group, and Lane.
#' @param group_name Character string used as the gene set label in the plot title.
#' Defaults to an empty string.
#' @param save_loc Path to the directory where the summary boxplot will be saved.
#' Defaults to the working directory.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{boxplot}: A ggplot2 boxplot of per-sample expression for the gene set.
#'   \item \code{results_df}: A data frame with the geoscore as the first row followed
#'     by per-gene normalized expression. Genes not found in the data appear as NA rows.
#' }
#'
#' @author Yohan Lefol
#' @export
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

  #Calculate geoscore using log-mean (robust to zeros and large values)
  geo_score <- exp(colMeans(log(pmax(sub_dta, 1))))

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
#' plots_and_table_wrapper
#'
#' @description
#' Generates a per-patient clinical report combining a quantile plot and a gene
#' expression summary table. For each patient, the geoscore is plotted against
#' control-derived quantile bands, and a table of per-gene normalized counts is
#' shown alongside. Multi-timepoint patients are handled automatically. Plots are
#' both printed to the active device (for RMarkdown rendering) and saved to disk.
#'
#' @param res_list The list returned by \code{analysis_with_geneset}, containing
#' \code{results_df} with geoscores and per-gene expression.
#' @param samp_sheet A sample sheet data frame with columns Sample, Group, and Lane.
#' @param norm_counts A gene x sample matrix of normalized counts with lane IDs as
#' column names.
#' @param target_genes Character vector of genes to display in the summary table.
#' @param save_loc Path to the directory where per-patient plots will be saved.
#'
#' @return
#' Called for its side effects. Prints one combined plot per patient to the active
#' device and saves each as a PNG file named by patient ID in \code{save_loc}.
#' Also prints control quantile values to the console.
#'
#' @author Yohan Lefol
#' @export
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
    ggsave(filename = paste0(save_loc,target_patient,'.png'),plot = combined_plot,dpi=600,bg='white')

    #Page break for markdown
    if(target_patient != all_patients[length(all_patients)]){#If it isn't the last patient
      cat('\\pagebreak')
    }
  }
}


#' @title
#' create_quantile_plot
#'
#' @description
#' Creates a ggplot2 visualization showing a patient's geoscore trajectory plotted
#' against control-derived quantile reference bands (0-25, 25-50, 50-75, 75-100
#' percentiles). Controls are shown as a filled background; patient samples are
#' shown as labeled points, connected by a line when multiple timepoints exist.
#'
#' @param plot_dta A data frame containing columns \code{geoscore_ISG}, \code{Group},
#' \code{Sample}, and \code{timepoint}. Should include all controls and the target
#' patient's rows.
#' @param quantiles A named numeric vector of five quantile values (0\%, 25\%, 50\%,
#' 75\%, 100\%) derived from the control group, as returned by \code{quantile()}.
#'
#' @return
#' A ggplot2 object. Print or pass to \code{grid.arrange()} for rendering.
#'
#' @author Yohan Lefol
#' @export
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


#' @title
#' plot_PCA
#'
#' @description
#' Creates a PCA plot from a normalized count matrix. Counts are log2-transformed
#' (pseudocount of 1 added), and PCA is computed on the top most-variable genes.
#' Samples are labeled and colored by group.
#'
#' @param normalized_counts A gene x sample matrix of normalized counts with sample
#' IDs as column names.
#' @param QCdata A data frame of QC metrics as returned by
#' \code{QC_plus_nSolver_processing} or \code{QC_plus_RUV_processing}, with
#' sample IDs as row names and a \code{Group} column.
#' @param n_genes Integer. Number of top most-variable genes to use for PCA.
#' Defaults to 500; if fewer genes are available all genes are used.
#'
#' @return
#' A ggplot2 object showing PC1 vs PC2 with percentage variance explained on each
#' axis and samples labeled by sample ID.
#'
#' @author Yohan Lefol
#' @export
plot_PCA <- function(normalized_counts, QCdata, n_genes = 500) {
  log_counts <- log2(normalized_counts + 1)

  # Select top n_genes most variable genes
  gene_vars <- apply(log_counts, 1, var)
  n_use <- min(n_genes, nrow(log_counts))
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_use]
  log_counts_sub <- log_counts[top_genes, ]

  # PCA on transposed matrix (samples as rows)
  pca_result <- prcomp(t(log_counts_sub), scale. = FALSE, center = TRUE)

  # Variance explained
  var_explained <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)

  # Build plot data frame
  pca_df <- data.frame(
    PC1   = pca_result$x[, 1],
    PC2   = pca_result$x[, 2],
    Sample = rownames(pca_result$x),
    Group  = QCdata$Group[match(rownames(pca_result$x), rownames(QCdata))]
  )

  ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3) +
    scale_color_manual(values = c('Kontroll' = '#00bfc4', 'Pasient' = '#f8766d')) +
    labs(
      x     = paste0('PC1 (', var_explained[1], '%)'),
      y     = paste0('PC2 (', var_explained[2], '%)'),
      title = paste0('PCA - top ', n_use, ' most variable genes')
    ) +
    theme_minimal()
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
##' \dontrun{
##' rcc <- system.file("extdata", "RCC", "20140604_C1-unstim_C1-unstim_01.RCC", package="NanoStringQCPro")
##' rcc.ls <- readRcc(rcc)
##' }
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





