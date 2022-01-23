## the following codes come from
## https://cloud.tencent.com/developer/article/1605975
library(data.table)
library(stringr)
library(tidyverse)
library(DSS)

## use the data from GSE52140
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52140
## download the row data
## https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE52140&format=file
## we only use two samples with replicates.
allDat <- lapply(list.files('./GSE52140/',pattern='*cpgs.txt.gz'),function(f){
  # f="GSM1084238_A0R_d0_rep1.cpgs.txt.gz"
  # f="GSM1084239_A3R_d0_rep1.cpgs.txt.gz"
  # f="GSM1251236_A0R_d0_rep2.cpgs.txt.gz"
  # f="GSM1251237_A3R_d0_rep2.cpgs.txt.gz"
  print(f)
  tmp=fread(file.path('GSE52140/',f))
  chr=as.character(tmp$CHR)
  pos=as.character(tmp$POS)
  newTmp=separate(tmp,col =3,into = c("methy", "unmethy"), sep = "/")
  newTmp$methy <- as.numeric(newTmp$methy)
  newTmp$all=as.numeric(newTmp$methy)+as.numeric(newTmp$unmethy)
  newTmp=as.data.frame(newTmp[,c(1,2,5,3)])
  colnames(newTmp)=c('chr', 'pos' ,'N' ,'X')
  return(newTmp)
})

sn=gsub('.cpgs.txt.gz','',list.files('GSE52140/',pattern='*cpgs.txt.gz'))
sn=gsub('GSM.*?_','',sn)


Human_BSobj <- makeBSseqData(allDat,sn)
Human_BSobj <- Human_BSobj[c(1:10000,1000000:1010000,2000000:2010000)]

dmlTest <- DMLtest(Human_BSobj, group1=c("A0R_d0_rep1","A0R_d0_rep2"),
                   group2=c("A3R_d0_rep1","A3R_d0_rep2"),smoothing=T)

dmls <- callDML(dmlTest, p.threshold=0.001)
Human_dmR <- callDMR(dmlTest, p.threshold=0.01)

usethis::use_data(Human_dmR,Human_BSobj)

