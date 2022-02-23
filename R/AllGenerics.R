##' getBaseModificationDf methods generics
##'
##'
##' @docType methods
##' @name getBaseModificationDf
##' @rdname getBaseModificationDf-methods
##' @importFrom methods setGeneric
##' @export
setGeneric("getBaseModificationDf",
           function(region,
                    input,
                    BSgenome,
                    strand = TRUE,
                    base = NULL,
                    motif = NULL,
                    position_bias = NULL,
                    ...){

             standardGeneric("getBaseModificationDf")

           })
