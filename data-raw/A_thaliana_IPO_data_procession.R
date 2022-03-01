# IPO and IPD ratio are important indicators in single-molecule
# real-time sequencing. We can distinguish modified bases with it.
#
#
# [1] Flusberg B A, Webster D R, Lee J H, et al. Direct
# detection of DNA methylation during single-molecule,
# real-time sequencing[J]. Nature methods, 2010, 7(6): 461-465.

# For the sake of convenience, we simulated IPD and IPD ratio from GSE62206
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62206
# These simulated data have no biological meanings at all and is designed
# for demonstration of visualization only.


library(data.table)
library(stringr)
library(tidyverse)
library(DSS)
## use the data from GSE62206
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62206
## download the row data
## and select two samples
allDat <- lapply(list.files('./GSE62206/',pattern='*txt.gz')[c(1:2)],function(f){
  # f="GSM1522201_WT_seedlings_BS_seq_genome-wide_CX_report.txt.gz"
  # f="GSM1522202_ddm1_seedlings_BS_seq_genome-wide_CX_report.txt.gz"
  print(f)
  tmp=fread(file.path('GSE62206/',f))
  tmp$all <- as.numeric(tmp$V5)+as.numeric(tmp$V4)
  tmp$V4 <- as.numeric(tmp$V4)
  tmp <- tmp[c(20000000:20040000),c(1,2,8,4)]
  colnames(tmp)=c('chr', 'pos' ,'N' ,'X')
  return(tmp)
})


## here we had read the GSE62206 and organized it in to a list.
## Now we assign random data to IPO and IPO ratio.
set.seed(706)
for (i in seq_len(length(allDat))) {

  IPD <- round(runif(length(allDat[[i]]$chr),5,20))
  IPD_ratio<- signif(runif(length(allDat[[i]]$chr),0.5,3),2)
  allDat[[i]]$N <- IPD
  allDat[[i]]$X <- IPD_ratio
  colnames(allDat[[i]])[3:4] <- c("IPD","IPD_ratio")

}

## Now we output the allDat into compressed txt file and bed file
## in order to test the BMplot::makebmData()
write.table(allDat[[1]],file = "GSM1522201_WT_seedlings_simulated.txt",
            row.names = F)
system("gzip ./GSM1522201_WT_seedlings_simulated.txt")

write.table(allDat[[2]],file = "GSM1522202_ddm1_seedlings_simulated.txt",
            row.names = F)
system("gzip ./GSM1522202_ddm1_seedlings_simulated.txt")

## And put the two files above together in folder "GSE62206_simulated"

## Now started to make bed files
allDatbed <- allDat
for (i in seq_len(length(allDat))) {

  allDatbed[[i]] <- mutate(allDat[[i]],end=allDat[[i]]$pos) %>%
    select(chr,pos,end,everything())

}

write.table(allDatbed[[1]],file = "GSM1522201_WT_seedlings_simulated.bed",
            row.names = F,col.names = F, sep = "\t")
system("gzip ./GSM1522201_WT_seedlings_simulated.bed")

write.table(allDatbed[[2]],file = "GSM1522202_ddm1_seedlings_simulated.bed",
            row.names = F,col.names = F, sep = "\t")
system("gzip ./GSM1522202_ddm1_seedlings_simulated.bed")

## And put the two bed files into folder "GSE62206_simulated_bed"

## now make the bmData object
A_thaliana_bmData <- makebmDataFromFiles(name = "GSE62206_simulated",
                                         variablesNames = c("IPD","IPD_ratio"),
                                         sampleNames = c("GSM1522201_WT",
                                                         "GSM1522202_ddm1"))
