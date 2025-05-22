#' getBaseModificationDf methods generics
#'
#'
#' @docType methods
#' @name getBaseModificationDf
#' @rdname getBaseModificationDf-methods
#' @importFrom methods setGeneric
#' @export
setGeneric("getBaseModificationDf",
           function(region,
                    input,
                    BSgenome,
                    base = NULL,
                    motif = NULL,
                    position_bias = NULL,
                    ...){

             standardGeneric("getBaseModificationDf")

           })

#' makebmDataFromData method generics
#'
#'
#' @docType methods
#' @name makebmDataFromData
#' @rdname makebmDataFromData-methods
#' @importFrom methods setGeneric
#' @export
setGeneric("makebmDataFromData", function(data,
                                          sampleNames=NULL){
  standardGeneric("makebmDataFromData")
})
