##' Arabidopsis thaliana BSseq object created by DSS
##'
##' A data set contains methylation information
"A_thaliana_BSobj"

##' Different methylation region created by DSS
##'
##' A data set contains information of different methylation region
##' @format A data frame with 1 row and 9 variables
##' \describe{
##'   \item{chr}{chromosome, the chromosome information of dmR}
##'   \item{start}{the start site of dmR}
##'   \item{end}{the end site of dmR}
##'   \item{length}{the length of dmR}
##'   \item{nCG}{Number of CpG sites contained in the DMR}
##'   \item{meanMethy1}{Average methylation levels in two conditions}
##'   \item{meanMethy2}{Average methylation levels in two conditions}
##'   \item{diff.Methy}{The difference in the methylation levels between two conditions}
##'   \item{areaStat}{	The sum of the test statistics of all CpG sites within the DMR}
##' }
"A_thaliana_dmR"


##' Human BSseq object created by DSS
##'
##' A data set contains methylation information
"Human_BSobj"

##' Different methylation region created by DSS
##'
##' A data set contains information of different methylation region
##' @format A data frame with 29 row and 9 variables
##' \describe{
##'   \item{chr}{chromosome, the chromosome information of dmR}
##'   \item{start}{the start site of dmR}
##'   \item{end}{the end site of dmR}
##'   \item{length}{the length of dmR}
##'   \item{nCG}{Number of CpG sites contained in the DMR}
##'   \item{meanMethy1}{Average methylation levels in two conditions}
##'   \item{meanMethy2}{Average methylation levels in two conditions}
##'   \item{diff.Methy}{The difference in the methylation levels between two conditions}
##'   \item{areaStat}{	The sum of the test statistics of all CpG sites within the DMR}
##' }
"Human_dmR"


##' Simulated arabidopsis thaliana BSseq object created by DSS
##'
##' A simulated data set contains methylation information.
##' These data were created by finding the adenine positions in a specific region,
##' which was supposed to be the dmR, and assign the adenine sites with random values.
##' At last created BSseq object by DSS.
"simulated_BSobj"

##' Different methylation region created by simulation
##'
##' A data set contains information of different methylation region
##' @format A data frame with 1 row and 3 variables
##' \describe{
##'   \item{chr}{chromosome, the chromosome information of dmR}
##'   \item{start}{the start site of dmR}
##'   \item{end}{the end site of dmR}
##' }
"simulated_dmR"

##' Simulated IPO data of A_thaliana
##'
##' A data set contains simulated IPO information. This data set is simulated
##'    from GSE62206. We use two samples from GSE62206 and replace the previous value
##'    in GSE62206 with random IPO and IPO ratio.
"simulated_IPO"
