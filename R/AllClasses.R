#' BSseqExtra Class
#'
#' This class added extra data to BSseq class. Change the assays by storing
#'     M/Cov to any value1/2
#'
#' @name BSseqExtra-class
#' @docType class
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom methods setClass
#' @keywords classes
#' @seealso BSseqExtra class inherits RangedSummarizedExperiment class,
#'     other slots see \linkS4class{RangedSummarizedExperiment}
#' @exportClass BSseqExtra
setClass("BSseqExtra", contains = "RangedSummarizedExperiment")
