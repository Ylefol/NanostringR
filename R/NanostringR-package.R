#' NanostringR: Quality Control and Normalization for NanoString Gene Expression Data
#'
#' @description
#' This package consolidates two unmaintained packages --- NanoStringNorm
#' (Waggott et al., 2012) and NanoStringQCPro (Bioconductor) --- with custom
#' functions for clinical NanoString analysis workflows. The source code of those
#' packages is included verbatim; all credit for those contributions belongs to
#' their original authors.
#'
#' @keywords internal
"_PACKAGE"

## Exports for functions from NanoStringNorm (Waggott et al., 2012).
## These files are not authored by the package maintainer and cannot carry
## @export tags directly; exports are declared here instead.
#' @rawNamespace export(NanoStringNorm)

## Exports for functions from NanoStringQCPro (Bioconductor).
#' @rawNamespace export(readRcc)

## Exports for functions from nanostring_RUV_functions.R (NanoStringQCPro).
#' @rawNamespace export(RUV_total)
#' @rawNamespace export(imagingQC)
#' @rawNamespace export(bindingDensityQC)
#' @rawNamespace export(limitOfDetectionQC)
#' @rawNamespace export(positiveLinQC)
#' @rawNamespace export(makeRLEplot)
NULL
