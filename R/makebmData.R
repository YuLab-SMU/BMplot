##' makebmDataFromData method for \code{CompressedGRangesList} objects
##'
##' @name makebmDataFromData
##' @docType methods
##' @rdname makebmDataFromData-methods
##' @aliases makebmDataFromData,CompressedGRangesList-method
##' @param data lists object
##' @param sampleNames the name of each samples
##' @importFrom methods setMethod
##' @importFrom GenomicRanges seqnames
##' @importFrom GenomicRanges start
##' @importFrom dplyr mutate
##' @exportMethod makebmDataFromData
setMethod("makebmDataFromData", signature(data="CompressedGRangesList"),
          function(data,
                   sampleNames=NULL){

            data_list <- lapply(data, function(x){

              chr <- as.character(seqnames(x))
              pos <- start(x)
              valueNames <- names(mcols(x))

              n0 <- length(valueNames)

              df <- data.frame(chr=chr,pos=pos)

              mcol <- as.data.frame(mcols(x))
              df <- data.frame(df,mcol)

              colnames(df) <- c("chr","pos",valueNames)

            })

            makebmDataFromData.internal(data = data_list, sampleNames=sampleNames)

          })


##' makebmDataFromData method for \code{GRanges} objects
##'
##' @name makebmDataFromData
##' @docType methods
##' @rdname makebmDataFromData-methods
##' @aliases makebmDataFromData,GRanges-method
##' @param data lists object
##' @param sampleNames the name of each samples
##' @importFrom methods setMethod
##' @importFrom GenomicRanges seqnames
##' @importFrom GenomicRanges start
##' @importFrom dplyr mutate
##' @exportMethod makebmDataFromData
setMethod("makebmDataFromData", signature(data="GRanges"),
          function(data,
                   sampleNames=NULL){

            chr <- as.character(seqnames(data))
            pos <- start(data)
            valueNames <- names(mcols(data))

            n0 <- length(valueNames)

            df <- data.frame(chr=chr,pos=pos)

            mcol <- as.data.frame(mcols(data))
            df <- data.frame(df,mcol)

            colnames(df) <- c("chr","pos",valueNames)

            data_list <- list(df)

            makebmDataFromData.internal(data = data_list, sampleNames=sampleNames)


          })

##' makebmDataFromData method for \code{list} objects
##'
##' @name makebmDataFromData
##' @docType methods
##' @rdname makebmDataFromData-methods
##' @aliases makebmDataFromData,list-method
##' @param data lists object
##' @param sampleNames the name of each samples
##' @details The objects in \code{data} must have specific forms. Colunms should be
##'    features, which should be organized in the order of "chr", "pos", "value1",
##'    "value2(optional)". chr stands for chromosome. pos stands for position on
##'    chromosome, also known as coordinates. value1/2 stands for the value on each base.
##'    The colnames can be any character but must be in the order. Rows stands for each
##'    observation.
##' @importFrom methods setMethod
##' @exportMethod makebmDataFromData
setMethod("makebmDataFromData", signature(data="list"),
          function(data,
                   sampleNames=NULL){

            if(any(unlist(lapply(list, function(x) is.null(colnames(x)))))){

              stop("pls input colnames to each objects...")

            }

            makebmDataFromData.internal(data = data, sampleNames=sampleNames)

          })


##' makebmDataFromData method for \code{data.frame} objects
##'
##' @name makebmDataFromData
##' @docType methods
##' @rdname makebmDataFromData-methods
##' @aliases makebmDataFromData,data.frame-method
##' @param data lists object
##' @param sampleNames the name of each samples
##' @details The objects in \code{data} must have specific forms. Colunms should be
##'    features, which should be organized in the order of "chr", "pos", "value1",
##'    "value2(optional)". chr stands for chromosome. pos stands for position on
##'    chromosome, also known as coordinates. value1/2 stands for the value on each base.
##'    The colnames can be any character but must be in the order. Rows stands for each
##'    observation.
##' @importFrom methods setMethod
##' @exportMethod makebmDataFromData
setMethod("makebmDataFromData", signature(data="data.frame"),
          function(data,
                   sampleNames=NULL){

            if(is.null(colnames(data))){
              stop("pls input colnames...")
            }

            data_list <- list(data)

            makebmDataFromData.internal(data = data_list,
                                        sampleNames = sampleNames)

          })


##' make dmData object from data
##'
##' @title makebmDataFromData.internal
##' @rdname makebmDataFromData
##' @param data lists object
##' @param sampleNames the name of each samples
##' @details This internal function was inspired by DSS::makeBSseqData.
##'
##'    The objects in \code{data} must have specific forms. Colunms should be
##'    features, which should be organized in the order of "chr", "pos", "value1",
##'    "value2(optional)". chr stands for chromosome. pos stands for position on
##'    chromosome, also known as coordinates. value1/2 stands for the value on each base.
##'    The colnames can be any character but must be in the order. Rows stands for each
##'    observation.
##' @return dmData object
makebmDataFromData.internal <- function(data,
                                        sampleNames=NULL){

  n0 <- length(colnames(data[[1]]))-2

  ## check the sampleNames
  sampleNames <- .check_and_make_sampleNames(data = data, sampleNames = sampleNames)

  ## check variables names in data
  if(!.check_variables_names(data)){
    stop("pls input the same column names in each object in data...")
  }


  if(n0 == 1){
    bmData <- make_bmData_from_value1(data = data, sampleNames = sampleNames)
  }else{
    bmData <- make_bmData_from_value1_and_value2(data = data, sampleNames = sampleNames)
  }

  return(bmData)

}


##' @importFrom dplyr mutate
##' @importFrom dplyr arrange
make_bmData_from_value1 <- function(data, sampleNames){

  n0 <- length(data)
  allDat <- data.frame(data[[1]][,c(1:2)])
  valueNames <- colnames(data[[1]])[3]

  ## merge data
  for(i in 1:n0){
    allDat <- data.frame(allDat, data[[i]][,3])
  }

  colnames(allDat) <- c("chr","pos",sampleNames)

  ## order the allDat
  chr <- pos <- sampleNames <- NULL

  allDat.ordered <- arrange(allDat,chr,pos)

  value1 <- as.matrix(allDat.ordered[,sampleNames])


  bmData <- bmData(chr = allDat.ordered$chr,
                   pos = allDat.ordered$pos,
                   value1 = value1,
                   sampleNames = sampleNames,
                   valueNames = valueNames)

  return(bmData)

}

##' @importFrom dplyr mutate
##' @importFrom dplyr arrange
make_bmData_from_value1_and_value2 <- function(data, sampleNames){

  n0 <- length(data)
  allDat_value1 <- data.frame(data[[1]][,c(1:2)])
  allDat_value2 <- data.frame(data[[1]][,c(1:2)])
  value1_name <- colnames(data[[1]])[3]
  value2_name <- colnames(data[[1]])[4]

  ## merge data
  for(i in 1:n0){
    allDat_value1 <- data.frame(allDat_value1, data[[i]][,3])
    allDat_value2 <- data.frame(allDat_value2, data[[i]][,4])
  }

  colnames(allDat_value1) <- c("chr","pos",sampleNames)
  colnames(allDat_value2) <- c("chr","pos",sampleNames)

  ## order the allDat
  chr <- pos <- sampleNames <- NULL

  allDat_value1.ordered <- arrange(allDat_value1,chr,pos)
  allDat_value2.ordered <- arrange(allDat_value2,chr,pos)

  value1 <- as.matrix(allDat_value1.ordered[,sampleNames])
  value2 <- as.matrix(allDat_value2.ordered[,sampleNames])

  bmData <- bmData(chr = allDat_value1.ordered$chr,
                   pos = allDat_value1.ordered$pos,
                   value1 = value1,
                   value2 = value2,
                   sampleNames = sampleNames,
                   valueNames = c(value1_name,value2_name))

  return(bmData)

}


##‘ make bmData from files
##'
##' This function makes bmData object from files. Users can input
##' the name of a file or a file folder.
##'
##' @name makebmDataFromFiles
##' @param name the name of files or file folder
##' @param sampleNames the name for each file
##' @param variablesNames the names of the first two columns will be assigned c("chr","pos"),
##'     the names of the following columns will be assigned by variablesNames
##' @importFrom utils file_test
##' @details bed files and txt files are supported. Bed files can
##'    only contain no more than two metadata, as it stands for value1/2. Txt files
##'    should organize the columns as chr, pos, value1, value2(optional).
##' @export
makebmDataFromFiles <- function(name,
                                sampleNames = NULL,
                                variablesNames = NULL){

  ## check the name of file
  if(!is.character(name)) stop("pls input character as file name...")

  if(file_test("-d", name)){

    ## deal with file folder
    data <- makebmDataFromFiles.folder(name = name, variablesNames = variablesNames)

  }else{

    ## deal with file
    data <- makebmDataFromFiles.file(name = name, variablesNames = variablesNames)
  }

  result <- makebmDataFromData(data = data,
                               sampleNames = sampleNames)

}


##' @importFrom ChIPseeker readPeakFile
##' @importFrom GenomicRanges mcols
##' @importFrom GenomicRanges GRangesList
##' @importFrom data.table fread
##' @importFrom utils getFromNamespace
makebmDataFromFiles.folder <- function(name, variablesNames){

  ## check the file type
  file_type <- gsub(".*\\.","",list.files(name))

  if(file_type[!duplicated(file_type)] != 1){

    stop("pls input files with the same type in the folder...")

  }

  if(is.null(variablesNames)){
    cat(">> no variable name is assigned,default names will be assigned",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }

  isBedFile <- getFromNamespace("isBedFile","ChIPseeker")

  if(isBedFile(file_type)){

    data <- lapply(list.files(name),function(x){

      cat(">> reading",x, format(Sys.time(), "%Y-%m-%d %X"), "\n")

      tmp=readPeakFile(file.path(paste0(name,"/"),x))

      if(is.null(variablesNames)){
        n0 <- length(names(mcols(tmp)))
        variablesNames <- paste0("value",1:n0)
      }

      names(mcols(tmp)) <- variablesNames

      return(tmp)
    })

    data_list <- GRangesList(data)

    return(data_list)
  }

  data_list <- lapply(list.files(name),function(x){

    cat(">> reading",x, format(Sys.time(), "%Y-%m-%d %X"), "\n")
    tmp=fread(file.path(paste0(name,"/"),x))

    if(is.null(variablesNames)){
      n0 <- ncol(tmp)-2
      variablesNames <- paste0("value",1:n0)
    }

    colnames(tmp)=c('chr', 'pos' ,variablesNames)
    return(tmp)
  })

  return(data_list)

}


##' @importFrom ChIPseeker readPeakFile
##' @importFrom GenomicRanges mcols
##' @importFrom data.table fread
makebmDataFromFiles.file <- function(name, variablesNames){

  if(is.null(variablesNames)){
    cat(">> no variable name is assigned,default names will be assigned",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }


  isBedFile <- getFromNamespace("isBedFile","ChIPseeker")


  if(isBedFile(name)){

    cat(">> reading",name, format(Sys.time(), "%Y-%m-%d %X"), "\n")

    data=readPeakFile(file.path("./",name))

    if(is.null(variablesNames)){
      n0 <- length(names(mcols(data)))
      variablesNames <- paste0("value",1:n0)
    }

    names(mcols(data)) <- variablesNames

    return(data)
  }

  cat(">> reading",name, format(Sys.time(), "%Y-%m-%d %X"), "\n")
  data=fread(file.path("./",name))

  if(is.null(variablesNames)){
    n0 <- ncol(data)-2
    variablesNames <- paste0("value",1:n0)
  }

  colnames(data)=c('chr', 'pos' ,variablesNames)
  return(data)

}
