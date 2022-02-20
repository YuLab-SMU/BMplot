#' Class BMplot
#' This class is a subclass from BSseq class, storing the information of any
#' value of interest. Acting as a reference role for visualizing base modification
#'
#' @name BMplot-class
#' @aliases BMplot
#'
#' @docType class
#' @slot valuenames the name of each value stored in BMplot object
#' @importClassesFrom bsseq BSseq
#' @import methods
#' @keywords classes
#' @exportClass BMplot
.BMplot <- setClass("BMplot",
                    slots = representation(valuenames = c("character")),
                    contains = "BSseq")
