## the following codes were inspired by
## https://cloud.tencent.com/developer/article/1605975
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
  tmp <- tmp[,c(1,2,8,4)]
  colnames(tmp)=c('chr', 'pos' ,'N' ,'X')
  return(tmp)
})

sn=gsub('.txt.gz','',list.files('GSE62206/',pattern='*.txt.gz')[c(1:2)])
sn=gsub('GSM.*?_','',sn)
sn=gsub('_BS.*','',sn)

A_thaliana_BSobj <- makeBSseqData(allDat,sn)[c(20000000:20040000)]

dmlTest <- DMLtest(A_thaliana_BSobj, group1=c("WT_seedlings"),
                   group2=c("ddm1_seedlings"),smoothing = T)

dmls <- callDML(dmlTest, p.threshold=0.01)
A_thaliana_dmR <- callDMR(dmlTest, p.threshold=0.01)

usethis::use_data(A_thaliana_BSobj,A_thaliana_dmR)
