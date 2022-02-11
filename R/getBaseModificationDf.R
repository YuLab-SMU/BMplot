##' get the information of base modification
##'
##'
##' This function retrieve the information of each base.Then organized it to dataframe.
##' Base modification includes modification in DNA and RNA.
##'
##' @param region base modification region in the form of dataframe, having columns of "chr","start" and "end"
##' @param BSseq BSseq objects
##' @param BSgenome genome reference
##' @param strand distinguish strand information or not
##' @param cover_depth take the depth of cover into account or not
##' @param base one of A/T/G/C/U
##' @param motif the motif(e.g C:CG/CH, A:GAGG/AGG) of the base modification
##' @param position_bias 1-base bias. e.g position_bias = 1("C" in "CHH"), position_bias = 2("A" in "GAGG")
##' @return dataframe
##' @export
getBaseModificationDf <- function(region,
                                  BSseq,
                                  BSgenome,
                                  strand = TRUE,
                                  cover_depth=TRUE,
                                  base = NULL,
                                  motif = NULL,
                                  position_bias = NULL){

  list_flag <- FALSE
  if(nrow(region) != 1){
    list_flag <- TRUE
  }

  ## check and assign the base
  if(is.null(base)){
    warning("No base was specified, \"C\"(cytosine) will be detected...")
    base <- "C"

  }else{

    if(!all(base %in% c("A","T","G","C","U"))){
      stop("please input legal base...")
    }
  }

  ## assign the motif
  if(is.null(motif)){

    ## assign the motif for cytosine
    if(base == "C"){
      motif <- c("CG","CHH","CHG")
      warning("No motif was specified, \"CG\",\"CHH\",\"CHG\" would be detected by default...")
    }

    ## assign the motif for adenine
    if(base == "A"){
      motif <- c("AGG","AGAA","TAGG")
      warning("No motif was specified, \"AGG\",\"AGAA\",\"TAGG\" would be detected by default...")
    }

    ## there is no default motif for Thymine, Uridine and Guanine
    if(base == "T" || base == "G" || base == "U"){
      stop("There is no default motif for Thymine, Uridine and Guanine.Please specify the motif...")
    }
  }

  ## we substitute the U with T in order to fit in the RNA situation
  motif <- gsub("U","T",motif)

  ## assign the position bias
  if(is.null(position_bias)){

    ## create the position bias by motif
    position_bias <- rep(1,length(motif))

  }else{

    if(length(position_bias) != length(motif)){
      stop("The length of position bias and motif are not equal...")
    }

  }

  ## check the position motif
  for (i in seq_len(length(motif))) {

    tmp <- unlist(strsplit(motif[i],split = "")[[1]])
    if(tmp[position_bias[i]] != base){
      stop("The ", position_bias[i]," position of the ",motif[i]," is not the correct base(",base,") to be detected.",
           "please cheak the position bias...")
    }
  }

  ## name the position bias
  names(position_bias) <- motif

  if(list_flag){

    df <- list()
    for (i in seq_len(nrow(region))) {
      df[[i]] <- getBaseModificationDf.internal(region = region[i,],
                                                BSseq = BSseq,
                                                BSgenome = BSgenome,
                                                strand = strand,
                                                cover_depth = cover_depth,
                                                motif = motif,
                                                base = base,
                                                position_bias = position_bias)
    }

  }else{
    df <- getBaseModificationDf.internal(region = region,
                                         BSseq = BSseq,
                                         BSgenome = BSgenome,
                                         strand = strand,
                                         cover_depth = cover_depth,
                                         motif = motif,
                                         base = base,
                                         position_bias = position_bias)
  }

  return(df)
}


##' @importFrom tidyr gather
##' @importFrom tidyselect all_of
##' @importFrom GenomicRanges mcols
##' @importFrom GenomicRanges start
getBaseModificationDf.internal <- function(region,
                                           BSseq,
                                           BSgenome,
                                           strand,
                                           cover_depth,
                                           motif,
                                           base,
                                           position_bias){


  ## load the reference of BSgenome
  BSgenome <- loadBSgenome(BSgenome)

  ## make the methylation from bsseq object
  methylation_reference <- make_Methylation_reference(BSseq,cover_depth)

  ## extract methylation information for region object
  dmR_melth <- detect_strand_and_motif(region = region,
                                       motif = motif,
                                       BSgenome = BSgenome,
                                       strand = strand,
                                       methylation_reference = methylation_reference,
                                       BSseq = BSseq,
                                       base = base,
                                       position_bias = position_bias)

  ## make the dataframe for plotting
  results <- as.data.frame(mcols(dmR_melth))
  results$coordinate <- start(dmR_melth)
  results$strand <- as.character(strand(dmR_melth))
  coordinate <- c(start(dmR_melth)[1]:start(dmR_melth)[length(start(dmR_melth))])


  ## fill in the non-mentioned site
  coordinate <- data.frame(coordinate = coordinate)
  results <- merge(results,coordinate,by = "coordinate",all=T)
  results$motif[is.na(results$motif)] <- "none"
  results[is.na(results)] <- 0

  ## global binding for value
  value <- NULL

  if (cover_depth) {
    depth_content <- names(results)[grep("depth",names(results))]
    methylation_content <- names(results)[grep("methylation",names(results))]
    content <- c(depth_content,methylation_content)
    df <- gather(results, sample, value,all_of(content))
    df$type[grep("depth",df$sample)] <- "depth"
    df$type[grep("methylation",df$sample)] <- "methylation"
    df$sample <- gsub("_depth","",df$sample)
    df$sample <- gsub("_methylation","",df$sample)
  }else{
    content <- names(results)[grep("methylation",names(results))]
    df <- gather(results, sample, value,all_of(content))
    df$sample <- gsub("_methylation","",df$sample)
  }

  df$strand[df$strand == 0] <- "*"

  ## assign attributes to df
  attr(df,"cover_depth") <- cover_depth
  attr(df,"chromosome") <- paste0("Chr",gsub("chr","",region$chr,ignore.case = T))
  return(df)
}

