library(BMplot)

test_that("getBaseModificationDf() return a data.frame", {

  skip_if_not_installed("BSgenome.Athaliana.TAIR.TAIR9")

  library(BSgenome.Athaliana.TAIR.TAIR9)

  expect_true(is(getBaseModificationDf(region = A_thaliana_dmR[1,],
                                       BSseq = A_thaliana_BSobj,
                                       BSgenome = BSgenome.Athaliana.TAIR.TAIR9::Athaliana,
                                       motif = c("CHH","CHG"),
                                       base = "C"), "data.frame"))
})
