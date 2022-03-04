## this script show the procession of making simulated bed file
## the data was drived from ChIPseeker,by using getSampleFiles()
## the file stored in inst/extdata/
library(ChIPseeker)
library(GenomicRanges)

# make bmData from grange object
files <- getSampleFiles()

set.seed(706)

gr1 <- readPeakFile(files[[4]])
mcols(gr1)[[1]] <- signif(runif(length(gr1)),2)

gr2 <- readPeakFile(files[[5]])
mcols(gr2)[[1]] <- signif(runif(length(gr2)),2)

df1 <- data.frame(chr=as.character(seqnames(gr1)),
                  start=start(gr1), end=start(gr1),
                  V4=mcols(gr1)[[1]], V5=mcols(gr1)[[2]])

write.table(df1,file = "GSM1295076_CBX6_BF_simulated.bed",
            row.names = F, col.names = F, quote = F, sep = "\t")

system("gzip ./GSM1295076_CBX6_BF_simulated.bed")


df2 <- data.frame(chr=as.character(seqnames(gr2)),
                  start=start(gr2), end=start(gr2),
                  V4=mcols(gr2)[[1]], V5=mcols(gr2)[[2]])

write.table(df2,file = "GSM1295077_CBX7_BF_simulated.bed",
            row.names = F, col.names = F, quote = F, sep = "\t")

system("gzip ./GSM1295077_CBX7_BF_simulated.bed")

## put the two bed file into simulated_bed_file/

## make txt file
df1_1 <- data.frame(chr=as.character(seqnames(gr1)),
                    start=start(gr1),
                    V4=mcols(gr1)[[1]], V5=mcols(gr1)[[2]])

write.table(df1_1,file = "GSM1295076_CBX6_BF_simulated.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

system("gzip ./GSM1295076_CBX6_BF_simulated.txt")


df2_2 <- data.frame(chr=as.character(seqnames(gr2)),
                    start=start(gr2),
                    V4=mcols(gr2)[[1]], V5=mcols(gr2)[[2]])

write.table(df2_2,file = "GSM1295077_CBX7_BF_simulated.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

system("gzip ./GSM1295077_CBX7_BF_simulated.txt")


## put the two txt file into simulated_txt_file
