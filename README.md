<!-- README.md is generated from README.Rmd. Please edit that file -->

# BMplot: Base Modification Visulization

<!-- badges: start -->
<!-- badges: end -->

This package implements functions to visualize the profile of base
modification in both DNA and RNA.

## Introduction

Epigenetics plays an important role in Biology in the way of changing
genetic information without changing nucleic acids. Histone modification
and Methylation are two vital mechanisms in Epigenetics.

Many methods are used to detect the methylation in Epigenentics with the
combination of sequencing and chemical methods. `Me‑DIP` and `oxBS‑seq`
for `5mC`. `TAB-seq` and `hme‑DIP` for `5hmC`. `6mA‑DIP` and
`6mA‑RE-seq` for `6mA`. These methods can be divided into two
categories, including `Antibody enrichment techniques` for base
fragments and `Single-base resolution techniques` for single base.
BMplot is designed to visualize the `Single-base resolution techniques`.

With the data created by the methods mentioned above, BMplot can map the
data into the sequence data stored in the set of BSgenome packages and
visualize it with great concern to the `strand` and `motif` information.

## Workflow

The BMplot provide two functions to visualize the base modification.
`getBaseModificationDf()` gets the information of base modification in
each base and organizes it into a data frame.
`plotBaseModificationProf()` can return a ggplot object to visualize the
base modification through the data frame created by
`getBaseModificationDf()`.

DNA Methylation is one of the base modifications. It is well worthy to
investigate the different methylation region(dmR) in DNA methylayion.

### Detect dmR and Create BSseq object

Any method can be used to detect dmR and create BSseq object. We use
`DSS` package to achieve this goal. For better experience, We made some
example data for users to go through BMplot.

### Get the information of base modification

We extracted the information from BSseq object and mapped to the
sequence data stored in `BSgenome.Hsapiens.UCSC.hg19`. Then we extracted
the motif and strand information and organized it in the form of data
frame. Cytosine was the interesting base to detect and the
`CG, CHH and CHG` were the motifs to investigate.

``` r
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome_human <- BSgenome.Hsapiens.UCSC.hg19
human_df <- getBaseModificationDf(region = Human_dmR[1,],
                                  BSgenome = BSgenome_human,
                                  BSseq = Human_BSobj,
                                  base = "C",
                                  motif = c("CG"))
```

The data frame have six columns. `coordinate` stands for the the
coordinate of each base in the region. `motif` stands for the motif of
the methylation. `strand` stands for the strand information of the base
modification. `sample` stands for the sample name of the base
modification sample. `value` stands for the coverage depth of each base
or the methylation value of each base. `type` decides the category
should be depth value or methylation value.

### Visualize the base modification

``` r
plotBaseModificationProf(human_df)
```

![](./vignettes/figures/readme.png)

## Installation

You can install the development version of BMplot from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("YuLab-SMU/BMplot")
```
