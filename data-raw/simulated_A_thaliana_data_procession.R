## There is no much sequencing data about adenine modification
## So this script simulated some data about adenine modification
## This simulation is only about visualization rather than biological meaning
library(BSgenome.Athaliana.TAIR.TAIR9)
BSgenome <- BSgenome.Athaliana.TAIR.TAIR9::Athaliana

sequence <- DNAStringSet(BSgenome[[5]],
                         start = 10000,
                         width = 2000)

A_postition <- which(strsplit(as.character(sequence),split="")[[1]] == "A") + 10000

normal <- data.frame(chr=rep("Chr5",400))
normal$pos <- sample(A_postition,400)
normal <- normal[order(normal$pos),]
normal$N <- sample(c(10:50),400,replace = TRUE)
normal$X <- round(sample(seq(0,1,by = 0.005),400,replace = TRUE)*normal$N)

test <- data.frame(chr=rep("Chr5",400))
test$pos <- sample(A_postition,400)
test <- test[order(test$pos),]
test$N <- sample(c(10:50),400,replace = TRUE)
test$X <- round(sample(seq(0,1,by = 0.005),400,replace = TRUE)*test$N)

allDat <- list(normal,test)
sn <- c("normal","test")
library(DSS)
simulated_BSobj <- makeBSseqData(allDat,sn)
simulated_dmR <- data.frame(chr="Chr5",start=10000,end=12000)

usethis::use_data(simulated_BSobj,simulated_dmR,overwrite = T)
