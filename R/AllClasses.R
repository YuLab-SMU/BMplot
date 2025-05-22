#' bmData Class
#'
#' This class added extra data to BSseq class. Change the assays by storing
#'     M/Cov to any value1/2
#'
#' @name bmData-class
#' @docType class
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom methods setClass
#' @keywords classes
#' @seealso bmData class inherits RangedSummarizedExperiment class,
#'     other slots see \linkS4class{RangedSummarizedExperiment}
#' @exportClass bmData
setClass("bmData", contains = "RangedSummarizedExperiment")
