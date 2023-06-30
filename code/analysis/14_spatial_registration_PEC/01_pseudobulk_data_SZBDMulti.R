library("zellkonverter")
library("SingleCellExperiment")
library("spatialLIBD")
library("jaffelab")
library("here")
library("sessioninfo")

#### load dataset  ####
args <- commandArgs(trailingOnly = TRUE)
# print(args)
version <- args[1]
input_file <- args[2]
message("input = ",version, " ", input_file)

dataset <- ss(input_file, "_")
message("\n#### Running: ", dataset, " ####")

## for V5 data
filepath <- here("raw-data", "psychENCODE", version, input_file)

stopifnot(file.exists(filepath))

message(Sys.time(), " - Reading data from: ", filepath)
sce <- readH5AD(file = filepath)

message("\nSCE Dimesions:")
dim(sce)

print(colnames(colData(sce)))
message("Cell Types:")
## must be syntactically valid
colData(sce)$cellType <- as.factor(make.names(colData(sce)$subclass))
table(sce$cellType)

# identify annotation/cluster labels
rowData(sce)$gene_name <- rownames(sce) # save gene name as column of rowData
rownames(sce) <- rowData(sce)$featureid # have to make row names of object the ensembl id instead of gene names

## Logcounts
# default “X” contain the log-normalized counts
message(Sys.time(), " revert to counts")

## check for all 0s (just first 100 cols for mem)
stopifnot(any(assays(sce)$X[, 1:100] != 0))

counts(sce) <- assays(sce)$X # change to normalized counts
# counts(sce)[counts(sce) != 0] <- (2^counts(sce)[counts(sce) != 0])-1 # Replace just non-zero values
counts(sce)@x <- 2^(counts(sce)@x) - 1 ## remove log2(counts + 1)

#### Pseudobulk ####
message(Sys.time(), " Pseudobulk")
sce_pseudo <- registration_pseudobulk(sce,
    var_registration = "cellType",
    # var_sample_id = "sampleID",
    var_sample_id = "individualID", ## for SZDBMulti-Seq
    covars = NULL
)

message("\nSCE Pseudobulk Dimesions:")
dim(sce_pseudo)

## Save results
saveRDS(sce_pseudo,
    file = here(
        "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
        paste0("pseudobulk_", dataset, ".rds")
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.2 Patched (2023-01-31 r83743)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-02-01
# pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version   date (UTC) lib source
# AnnotationDbi            1.60.0    2022-11-01 [2] Bioconductor
# AnnotationHub            3.6.0     2022-11-01 [2] Bioconductor
# assertthat               0.2.1     2019-03-21 [2] CRAN (R 4.2.1)
# attempt                  0.3.1     2020-05-03 [2] CRAN (R 4.2.1)
# basilisk                 1.10.2    2022-11-08 [2] Bioconductor
# basilisk.utils           1.10.0    2022-11-01 [2] Bioconductor
# beachmat                 2.14.0    2022-11-01 [2] Bioconductor
# beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
# benchmarkme              1.0.8     2022-06-12 [2] CRAN (R 4.2.1)
# benchmarkmeData          1.0.4     2020-04-23 [2] CRAN (R 4.2.1)
# Biobase                * 2.58.0    2022-11-01 [2] Bioconductor
# BiocFileCache            2.6.0     2022-11-01 [2] Bioconductor
# BiocGenerics           * 0.44.0    2022-11-01 [2] Bioconductor
# BiocIO                   1.8.0     2022-11-01 [2] Bioconductor
# BiocManager              1.30.19   2022-10-25 [2] CRAN (R 4.2.2)
# BiocNeighbors            1.16.0    2022-11-01 [2] Bioconductor
# BiocParallel             1.32.5    2022-12-23 [2] Bioconductor
# BiocSingular             1.14.0    2022-11-01 [2] Bioconductor
# BiocVersion              3.16.0    2022-04-26 [2] Bioconductor
# Biostrings               2.66.0    2022-11-01 [2] Bioconductor
# bit                      4.0.5     2022-11-15 [2] CRAN (R 4.2.2)
# bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
# bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# blob                     1.2.3     2022-04-10 [2] CRAN (R 4.2.1)
# bslib                    0.4.2     2022-12-16 [2] CRAN (R 4.2.2)
# cachem                   1.0.6     2021-08-19 [2] CRAN (R 4.2.1)
# cli                      3.6.0     2023-01-09 [2] CRAN (R 4.2.2)
# codetools                0.2-18    2020-11-04 [3] CRAN (R 4.2.2)
# colorspace               2.1-0     2023-01-23 [2] CRAN (R 4.2.2)
# config                   0.3.1     2020-12-17 [2] CRAN (R 4.2.1)
# cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
# crayon                   1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
# curl                     5.0.0     2023-01-12 [2] CRAN (R 4.2.2)
# data.table               1.14.6    2022-11-16 [2] CRAN (R 4.2.2)
# DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.2.1)
# dbplyr                   2.3.0     2023-01-16 [2] CRAN (R 4.2.2)
# DelayedArray             0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats       1.20.0    2022-11-01 [2] Bioconductor
# desc                     1.4.2     2022-09-08 [2] CRAN (R 4.2.1)
# digest                   0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
# dir.expiry               1.6.0     2022-11-01 [2] Bioconductor
# doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
# dotCall64                1.0-2     2022-10-03 [2] CRAN (R 4.2.1)
# dplyr                    1.1.0     2023-01-29 [2] CRAN (R 4.2.2)
# dqrng                    0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# DropletUtils             1.18.1    2022-11-22 [2] Bioconductor
# DT                       0.27      2023-01-17 [2] CRAN (R 4.2.2)
# edgeR                    3.40.2    2023-01-19 [2] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.2.1)
# ExperimentHub            2.6.0     2022-11-01 [2] Bioconductor
# fansi                    1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# fastmap                  1.1.0     2021-01-25 [2] CRAN (R 4.2.1)
# fields                   14.1      2022-08-12 [2] CRAN (R 4.2.1)
# filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
# foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
# fs                       1.6.0     2023-01-23 [2] CRAN (R 4.2.2)
# gargle                   1.3.0     2023-01-30 [2] CRAN (R 4.2.2)
# generics                 0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb           * 1.34.7    2023-01-24 [2] Bioconductor
# GenomeInfoDbData         1.2.9     2022-09-29 [2] Bioconductor
# GenomicAlignments        1.34.0    2022-11-01 [2] Bioconductor
# GenomicRanges          * 1.50.2    2022-12-16 [2] Bioconductor
# ggbeeswarm               0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
# ggplot2                  3.4.0     2022-11-04 [2] CRAN (R 4.2.2)
# ggrepel                  0.9.2     2022-11-06 [2] CRAN (R 4.2.2)
# glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# golem                    0.3.5     2022-10-18 [2] CRAN (R 4.2.1)
# googledrive              2.0.0     2021-07-08 [2] CRAN (R 4.2.1)
# gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
# gtable                   0.3.1     2022-09-01 [2] CRAN (R 4.2.1)
# HDF5Array                1.26.0    2022-11-01 [2] Bioconductor
# here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# htmltools                0.5.4     2022-12-07 [2] CRAN (R 4.2.2)
# htmlwidgets              1.6.1     2023-01-07 [1] CRAN (R 4.2.2)
# httpuv                   1.6.8     2023-01-12 [2] CRAN (R 4.2.2)
# httr                     1.4.4     2022-08-17 [2] CRAN (R 4.2.1)
# interactiveDisplayBase   1.36.0    2022-11-01 [2] Bioconductor
# IRanges                * 2.32.0    2022-11-01 [2] Bioconductor
# irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
# iterators                1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
# jaffelab               * 0.99.32   2022-12-07 [1] Github (LieberInstitute/jaffelab@7b7afe3)
# jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.2.1)
# jsonlite                 1.8.4     2022-12-06 [2] CRAN (R 4.2.2)
# KEGGREST                 1.38.0    2022-11-01 [2] Bioconductor
# knitr                    1.42      2023-01-25 [2] CRAN (R 4.2.2)
# later                    1.3.0     2021-08-18 [2] CRAN (R 4.2.1)
# lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.2.2)
# lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
# lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                    3.54.1    2023-01-26 [2] Bioconductor
# locfit                   1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# magick                   2.7.3     2021-08-18 [2] CRAN (R 4.2.1)
# magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# maps                     3.4.1     2022-10-30 [2] CRAN (R 4.2.2)
# MASS                     7.3-58.2  2023-01-23 [3] CRAN (R 4.2.2)
# Matrix                   1.5-3     2022-11-11 [2] CRAN (R 4.2.2)
# MatrixGenerics         * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats            * 0.63.0    2022-11-18 [2] CRAN (R 4.2.2)
# memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
# mime                     0.12      2021-09-28 [2] CRAN (R 4.2.1)
# munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
# nlme                     3.1-162   2023-01-31 [2] CRAN (R 4.2.2)
# paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.2.1)
# pillar                   1.8.1     2022-08-19 [2] CRAN (R 4.2.1)
# pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# pkgload                  1.3.2     2022-11-16 [2] CRAN (R 4.2.2)
# plotly                   4.10.1    2022-11-07 [2] CRAN (R 4.2.2)
# png                      0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
# promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.2.1)
# purrr                    1.0.1     2023-01-10 [2] CRAN (R 4.2.2)
# R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
# R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
# R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
# R6                       2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# rafalib                * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
# rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
# RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
# Rcpp                     1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                    1.98-1.10 2023-01-27 [2] CRAN (R 4.2.2)
# rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.2.1)
# restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
# reticulate               1.28      2023-01-27 [2] CRAN (R 4.2.2)
# rhdf5                    2.42.0    2022-11-01 [2] Bioconductor
# rhdf5filters             1.10.0    2022-11-01 [2] Bioconductor
# Rhdf5lib                 1.20.0    2022-11-01 [2] Bioconductor
# rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                    1.0.6     2022-09-24 [2] CRAN (R 4.2.1)
# roxygen2                 7.2.3     2022-12-08 [2] CRAN (R 4.2.2)
# rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# Rsamtools                2.14.0    2022-11-01 [2] Bioconductor
# RSQLite                  2.2.20    2022-12-22 [2] CRAN (R 4.2.2)
# rstudioapi               0.14      2022-08-22 [2] CRAN (R 4.2.1)
# rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
# rtracklayer              1.58.0    2022-11-01 [2] Bioconductor
# S4Vectors              * 0.36.1    2022-12-05 [2] Bioconductor
# sass                     0.4.5     2023-01-24 [2] CRAN (R 4.2.2)
# ScaledMatrix             1.6.0     2022-11-01 [2] Bioconductor
# scales                   1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
# scater                   1.26.1    2022-11-13 [2] Bioconductor
# scuttle                  1.8.3     2022-12-14 [1] Bioconductor
# segmented                1.6-2     2022-12-09 [1] CRAN (R 4.2.2)
# sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# shiny                    1.7.4     2022-12-15 [1] CRAN (R 4.2.2)
# shinyWidgets             0.7.6     2023-01-08 [2] CRAN (R 4.2.2)
# SingleCellExperiment   * 1.20.0    2022-11-01 [2] Bioconductor
# spam                     2.9-1     2022-08-07 [2] CRAN (R 4.2.1)
# sparseMatrixStats        1.10.0    2022-11-01 [2] Bioconductor
# SpatialExperiment      * 1.8.0     2022-11-01 [2] Bioconductor
# spatialLIBD            * 1.11.4    2022-12-16 [1] Github (LieberInstitute/spatialLIBD@951dded)
# statmod                  1.5.0     2023-01-06 [2] CRAN (R 4.2.2)
# stringi                  1.7.12    2023-01-11 [2] CRAN (R 4.2.2)
# stringr                  1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
# SummarizedExperiment   * 1.28.0    2022-11-01 [2] Bioconductor
# tibble                   3.1.8     2022-07-22 [2] CRAN (R 4.2.1)
# tidyr                    1.3.0     2023-01-24 [2] CRAN (R 4.2.2)
# tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# usethis                  2.1.6     2022-05-25 [2] CRAN (R 4.2.1)
# utf8                     1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                    0.5.2     2023-01-23 [2] CRAN (R 4.2.2)
# vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
# viridis                  0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
# viridisLite              0.4.1     2022-08-22 [2] CRAN (R 4.2.1)
# withr                    2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
# xfun                     0.37      2023-01-31 [2] CRAN (R 4.2.2)
# XML                      3.99-0.13 2022-12-04 [2] CRAN (R 4.2.2)
# xml2                     1.3.3     2021-11-30 [2] CRAN (R 4.2.1)
# xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.2.1)
# XVector                  0.38.0    2022-11-01 [2] Bioconductor
# yaml                     2.3.7     2023-01-23 [2] CRAN (R 4.2.2)
# zellkonverter          * 1.8.0     2022-11-01 [1] Bioconductor
# zlibbioc                 1.44.0    2022-11-01 [2] Bioconductor
#
# [1] /users/lhuuki/R/4.2.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
#
