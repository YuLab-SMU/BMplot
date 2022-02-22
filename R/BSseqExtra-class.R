#' Constructor for BSseqExtra objects
#'
#' This is constructor fo BSseqExtra objects.
#'
#' @param value1 the first value to be stored, a matrix-like object
#' @param value2 the second value to be stored, a matrix-like object
#' @param pos A vector of locations
#' @param chr A vector of chromosomes
#' @param gr An object of type \linkS4class{GRanges}
#' @param sampleNames A vector of sample names
#' @param valueNames the name of value1 or value2 or both. The order maps to the value.
#' @param ... other parameters from \code{\link[bsseq]{BSseq}}
#' @importFrom bsseq BSseq
#' @importFrom bsseq pData
#' @importFrom SummarizedExperiment assays<-
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges granges
#' @importFrom methods new
#' @importFrom S4Vectors SimpleList
#' @export
BSseqExtra <- function(value1 = NULL, value2 = NULL,
                       pos = NULL, chr = NULL, gr = NULL,
                       sampleNames = NULL, valueNames = NULL,
                       ...){

  ## check value1 and value2
  if(is.null(value1) && is.null(value2)){

    ## value1 and value2 can not be NULL simultaneously
    stop("Need an input value...")

  }

  ## check the coordination of valueNames and value1/2
  if(is.null(value1) && !is.null(value2)){

    if (length(valueNames) != 1) {
      stop("ValueNames do not match the value...")
    }

    value1 <- value2
    valueNames <- rep(valueNames,2)

  }


  if(is.null(value2) && !is.null(value1)){

    if (length(valueNames) != 1) {
      stop("ValueNames do not match the value...")
    }

    value2<- value1
    valueNames <- rep(valueNames,2)

  }

  if(!is.null(value1) && !is.null(value2)){

    if(length(valueNames) != 2){
      stop("ValueNames do not match the value...")
    }

  }

  ## This parameter check comes from BSseq()
  ## https://github.com/hansenlab/bsseq/blob/master/R/BSseq-class.R
  if (is.null(gr)) {
    if (is.null(pos) || is.null(chr)) {
      stop("Need 'pos' and 'chr' if 'gr' not supplied.")
    }
    gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1L))
  }
  if (!is(gr, "GRanges")) {
    stop("'gr' needs to be a GRanges.")
  }

  ## deal with the duplicated location
  if(any(duplicated(gr))){

    ## We removed duplicated locations and keep only one value
    warning("There are duplicated locis, which will be removed and keep only one value...")
    gr <- gr[!duplicated(gr)]
    value1 <- value1[!duplicated(gr)]
    value2 <- value2[!duplicated(gr)]

    ## Users can use other methods to deal with duplicated value.
    ## rowsums() maybe a good method.
    ## An example is placed in BSseq-class.R, following codes are from
    ## https://github.com/hansenlab/bsseq/blob/master/R/BSseq-class.R
    ##
    ##
    # # Collapse duplicate loci --------------------------------------------------
    #
    # is_duplicated <- duplicated(gr)
    # if (any(is_duplicated)) {
    #   warning("Detected duplicate loci. Collapsing counts in 'M' and 'Cov' ",
    #           "at these positions.")
    #   if (!is.null(coef) || !is.null(se.coef)) {
    #     stop("Cannot collapse when 'coef' or 'se.coef' are non-NULL.")
    #   }
    #   loci <- gr[!is_duplicated]
    #   ol <- findOverlaps(gr, loci, type = "equal")
    #   M <- rowsum(x = M, group = subjectHits(ol), reorder = FALSE)
    #   rownames(M) <- NULL
    #   Cov <- rowsum(x = Cov, group = subjectHits(ol), reorder = FALSE)
    #   rownames(Cov) <- NULL
    # } else {
    #   loci <- gr
    # }

  }


  ## In order to use BSseq() as an constructor, we make vitual M and Cov
  ## to pass the BSseq() check
  M <- matrix(c(rep(1,length(gr))),nrow = length(gr))
  Cov <- matrix(c(rep(3,length(gr))),nrow = length(gr))

  tmp_bsseq <- BSseq(M = M, Cov = Cov,
                     gr = gr, sampleNames = NULL,
                     ...)

  ## Now we extract the assays from BSseq object and
  ## substitute it with the value we input.
  if(identical(value1,value2)){

    command <- paste0("assays(tmp_bsseq) <- SimpleList(",
                      valueNames[1],
                      "=value1)")

    eval(parse(text = command))

  }else{

    command <- paste0("assays(tmp_bsseq) <- SimpleList(",
                      valueNames[1],
                      "= value1, ",
                      valueNames[2],
                      "= value2)")

    eval(parse(text = command))

  }

  new_se <- SummarizedExperiment(assays = assays(tmp_bsseq),
                                 rowRanges = granges(tmp_bsseq),
                                 colData = pData(tmp_bsseq))

  return(new("BSseqExtra", new_se))

}
