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
  colnames(tmp)=c('chr', 'pos' ,'IPO' ,'IPO_ratio')
  tmp$IPO_ratio <- signif(tmp$IPO_ratio/tmp$IPO,2)
  tmp$IPO <- round(tmp$IPO*0.1)
  return(tmp)
})


simulated_IPO <- allDat
use_data(simulated_IPO,overwrite = T)

