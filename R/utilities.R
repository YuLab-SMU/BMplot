make_Methylation_reference <- function(input,cover_depth){

  ## make the methylation reference from bsseq object
  methylation_reference <- input@rowRanges

  ## add the methylation situation to methylation_reference
  for (i in input@colData@rownames) {
    methylation_M <- input@assays@data@listData[["M"]][,i]

    ## fill the 0 in Cov with 1
    methylation_cov <- input@assays@data@listData[["Cov"]][,i]

    if (cover_depth) {
      command <- paste0("mcols(methylation_reference)$",
                        i,'_depth <- methylation_cov')

      eval(parse(text = command))
    }


    methylation_cov[input@assays@data@listData[["Cov"]][,i] == 0] <- 1

    ## the situation of methylation is calculated by M/cov
    methylation <- round(methylation_M/methylation_cov,3)

    command <- paste0("mcols(methylation_reference)$",
                      i,'_methylation <- methylation')

    eval(parse(text = command))

  }

  return(methylation_reference)

}


##' @importFrom SummarizedExperiment assayNames
make_reference <- function(input){

  ## make the  reference from bmData object
  reference <- input@rowRanges

  ## extract the valueNames
  aName <- assayNames(input)

  n0 <- length(aName)

  if(n0 == 1){

    for (i in input@colData@rownames) {

      value <- input@assays@data@listData[[aName]][,i]

      command <- paste0("mcols(reference)$",i, "_",aName,'<- value')

      eval(parse(text = command))

    }

  }else{

    for (i in input@colData@rownames) {

      aName1 <- aName[1]
      aName2 <- aName[2]

      value1 <- input@assays@data@listData[[aName1]][,i]
      value2 <- input@assays@data@listData[[aName2]][,i]

      command <- paste0("mcols(reference)$",i, "_",aName1,'<- value1')

      eval(parse(text = command))

      command <- paste0("mcols(reference)$",i, "_",aName2,'<- value2')

      eval(parse(text = command))

    }

  }

  return(reference)

}

##' @importFrom methods is
loadBSgenome <- function(BSgenome){

  if(!is(BSgenome,"BSgenome")){
    stop(">> input must be BSgenome object...")
  }

  if(is.null(BSgenome)){
    stop(">> please specify BSgenome object...")
  }else{

    ## get the object from BSgenome
    BSgenome_name <- attr(BSgenome,"pkgname")
    text <- paste0("package:",BSgenome_name)
    object <- ls(text)[grep("^[^BSgenome]",ls(text))]
    command <- paste0(BSgenome_name,"::",object)

    ## execute the command "pkg::object"
    BSgenome <- eval(parse(text = command))
  }

  return(BSgenome)
}



##' @importFrom GenomicRanges GRanges
##' @importFrom GenomicRanges findOverlaps
##' @importFrom GenomicRanges strand
##' @importFrom GenomicRanges strand<-
##' @importFrom GenomicRanges mcols
##' @importFrom GenomicRanges mcols<-
##' @importFrom GenomicRanges start
##' @importFrom Biostrings DNAStringSet
##' @importFrom IRanges IRanges
detect_strand_and_motif <- function(region,
                                    motif,
                                    BSgenome,
                                    methylation_reference,
                                    input,
                                    base,
                                    position_bias){

  ## detect the strand and motif information
  dmRregion <- GRanges(seqnames = region$chr,
                       ranges = IRanges(start = region$start,
                                        end = region$end))

  hits <- findOverlaps(methylation_reference,
                       dmRregion,
                       type = "within")

  dmR_melth <- methylation_reference[hits@from]


  ## make the chromosome
  region$chr <- gsub("chr","",region$chr,ignore.case = T)
  region$chr <- as.numeric(region$chr)

  ## detect the positive and negative strand
  positive_index <- rep(FALSE,length(dmR_melth@ranges@start))

  for(i in motif){

    melth_sequence <- DNAStringSet(BSgenome[[region$chr]],
                                   start = dmR_melth@ranges@start+position_bias[i]-1,
                                   width = 1)

    tmp_positive_index <- vapply(melth_sequence,
                                 function(x) grepl(base,as.character(x))==1,
                                 FUN.VALUE = logical(1))

    positive_index <- positive_index | tmp_positive_index
  }

  negetive_index <- !positive_index

  ## deal with the positive strand
  dmR_melth_positive <- dmR_melth[positive_index]
  strand(dmR_melth_positive) <- "+"

  for(i in motif){

    ## create the regex for the motif
    regex <- create_regex_patterns_positive(i)

    ## count the length of the motif
    width <- length(strsplit(i,split = "")[[1]])

    ## create the mapping sequence
    sequence <- DNAStringSet(BSgenome[[region$chr]],
                             start = dmR_melth@ranges@start[positive_index]-position_bias[i]+1,
                             width = width)

    ## get the index of the specific motif
    index <- vapply(sequence,
                    function(x) grepl(regex,as.character(x))==1,
                    FUN.VALUE = logical(1))

    ## fill in the motif columns
    mcols(dmR_melth_positive)[["motif"]][index] <- i
  }


  ## deal with the negetive strand
  dmR_melth_negetive <- dmR_melth[negetive_index]
  strand(dmR_melth_negetive) <- "-"

  for(i in motif){

    ## create the regex for the motif
    regex <- create_regex_patterns_negative(i)

    ## count the length of the motif
    width <- length(strsplit(i,split = "")[[1]])

    ## create the mapping sequence
    sequence <- DNAStringSet(BSgenome[[region$chr]],
                             end = dmR_melth@ranges@start[negetive_index]+position_bias[i]-1,
                             width = width)

    ## get the index of the specific motif
    index <- vapply(sequence,
                    function(x) grepl(regex,as.character(x))==1,
                    FUN.VALUE = logical(1))

    ## fill in the motif columns
    mcols(dmR_melth_negetive)[["motif"]][index] <- i

  }


  ## combine and reorder the results from positive and negative strand
  dmR_melth <- c(dmR_melth_positive,dmR_melth_negetive)
  dmR_melth <- dmR_melth[order(start(dmR_melth)),]

  ## filtered the unidentified melthylation sites
  dmR_melth <- dmR_melth[which(!is.na(mcols(dmR_melth)[["motif"]]))]

  return(dmR_melth)
}


##' create regex patterns in positive strand
##'
##' @param motif the motif(e.g C:CG/CH, A:GAGG/AGG) of the base modification
create_regex_patterns_positive <- function(motif){

  ## split the motif
  motif_list <- as.list(strsplit(motif,split = "")[[1]])

  ## check the character in the motif
  legal_character <- c("R","Y","M","K",
                       "S","W","H","B",
                       "V","D","N","A",
                       "T","G","C")

  if(!all(unlist(motif_list) %in% legal_character)){
    stop("please input legal motif...")
  }

  ## make the connection between degenerate bases and normal bases
  converted_motif <- lapply(motif_list, function(x){

    if(x == "R"){
      return("[AG]")
    }

    if(x == "Y"){
      return("[CT]")
    }

    if(x == "M"){
      return("[AC]")
    }

    if(x == "K"){
      return("[GT]")
    }

    if(x == "S"){
      return("[GC]")
    }

    if(x == "W"){
      return("[AT]")
    }

    if(x == "H"){
      return("[ATC]")
    }

    if(x == "B"){
      return("[GTC]")
    }

    if(x == "V"){
      return("[GAC]")
    }

    if(x == "D"){
      return("[GAT]")
    }

    if(x == "N"){
      return("[ATGC]")
    }

    if(x == "A"){
      return("A")
    }

    if(x == "T"){
      return("T")
    }

    if(x == "G"){
      return("G")
    }

    if(x == "C"){
      return("C")
    }

  })


  ## create the regex pattern
  regex <- paste0(unlist(converted_motif),collapse = "")

  return(regex)
}


##' create regex patterns in negative strand
##'
##' @param motif the motif(e.g C:CG/CH, A:GAGG/AGG) of the base modification
create_regex_patterns_negative <- function(motif){

  ## split the motif
  motif_list <- rev(as.list(strsplit(motif,split = "")[[1]]))

  ## check the character in the motif
  legal_character <- c("R","Y","M","K",
                       "S","W","H","B",
                       "V","D","N","A",
                       "T","G","C")

  if(!all(unlist(motif_list) %in% legal_character)){
    stop("please input legal motif...")
  }

  ## make the connection between degenerate bases and normal bases
  converted_motif <- lapply(motif_list, function(x){

    ## R->A/G->T/C
    if(x == "R"){
      return("[TC]")
    }

    ## Y->C/T->G/A
    if(x == "Y"){
      return("[GA]")
    }

    ## M->A/C->T/G
    if(x == "M"){
      return("[TG]")
    }

    ## K->G/T->C/A
    if(x == "K"){
      return("[CA]")
    }

    ## S->G/C->C/G
    if(x == "S"){
      return("[CG]")
    }

    ## W->A/T->T/A
    if(x == "W"){
      return("[TA]")
    }

    ## H->A/T/C->T/A/G
    if(x == "H"){
      return("[TAG]")
    }

    ## B->G/T/C->C/A/G
    if(x == "B"){
      return("[CAG]")
    }

    ## V->G/A/C->C/T/G
    if(x == "V"){
      return("[CTG]")
    }

    ## D->G/A/T->C/T/A
    if(x == "D"){
      return("[CTA]")
    }

    ## N->A/T/C/G->T/A/G/C
    if(x == "N"){
      return("[TAGC]")
    }

    ## A->A->T
    if(x == "A"){
      return("T")
    }

    ## T->T->A
    if(x == "T"){
      return("A")
    }

    ## G->G->C
    if(x == "G"){
      return("C")
    }

    ## C->C->G
    if(x == "C"){
      return("G")
    }

  })

  ## create the regex pattern
  regex <- paste0(unlist(converted_motif),collapse = "")

  return(regex)
}


##' @importFrom DSS makeBSseqData
##'
##' @export
DSS::makeBSseqData


##' @importFrom bsseq BSseq
##'
##' @export
bsseq::BSseq


.check_valueNames <- function(valueNames, n0){

  if(is.null(valueNames)){
    valueNames <- paste0("value",1:n0)
  }

  ## check the coordination of valueNames and value1/2
  if (length(valueNames) != n0) {
    stop("ValueNames do not match the value...")
  }

  return(valueNames)

}


.check_and_make_sampleNames <- function(data,sampleNames){

  n0 <- length(data)
  if(!is.null(names(data)) && !any(is.na(names(data)))){
    return(names(data))
  }

  if(is.null(sampleNames)){
    sampleNames <- paste("sample", 1:n0, sep="")
    return(sampleNames)
  }

  if(length(data) != length(sampleNames)){
    stop("sampleNames should have equal length with data...")
  }

  return(sampleNames)

}


.check_variables_names <- function(data){

  variables_names <- colnames(data[[1]])

  tmp <- unlist(lapply(data, function(x){
    if(!identical(variables_names,colnames(x))){
      return(TRUE)
    }

    return(FALSE)
  }))

  if(any(tmp)){
    return(FALSE)
  }else{
    return(TRUE)
  }

}


##' @importFrom SummarizedExperiment assays
##'
##' @export
SummarizedExperiment::assays

##' @importFrom SummarizedExperiment assay
##'
##' @export
SummarizedExperiment::assay
