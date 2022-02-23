##' getBaseModificationDf method for \linkS4class{BSseq}
##'
##' @docType methods
##' @rdname getBaseModificationDf-methods
##' @title getBaseModificationDf method
##' @param region base modification region in the form of dataframe, having columns of "chr","start" and "end"
##' @param input the input data stored in BSseq objects or BSseqExtra objects
##' @param BSgenome genome reference
##' @param strand distinguish strand information or not
##' @param base one of A/T/G/C/U
##' @param motif the motif(e.g C:CG/CH, A:GAGG/AGG) of the base modification
##' @param position_bias 1-base bias. e.g position_bias = 1("C" in "CHH"), position_bias = 2("A" in "GAGG")
##' @param cover_depth take the depth of cover into account or not
##' @param ... other parameters
##' @aliases getBaseModificationDf, BSseq-methods
##' @return dataframe
##' @importFrom methods setMethod
##' @exportMethod getBaseModificationDf
setMethod("getBaseModificationDf",signature(input = "BSseq"),
          function(region,
                   input,
                   BSgenome,
                   strand = TRUE,
                   base = NULL,
                   motif = NULL,
                   position_bias = NULL,
                   cover_depth=TRUE,
                   ...){

            getBaseModificationDf.BSseq(region = region,
                                        input = input,
                                        BSgenome = BSgenome,
                                        strand = strand,
                                        base = base,
                                        motif = motif,
                                        position_bias = position_bias,
                                        cover_depth = cover_depth)

          })
