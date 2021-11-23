## This script requires R 4.1
# module load conda_R/4.1.x

## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")

## vis
library("spatialLIBD")


## Define some info for the samples
sample_info <- data.frame(
  sample_id = c(
    "DLPFC_Br2743_ant_manual_alignment",
    "DLPFC_Br2743_mid_manual_alignment_extra_reads",
    "DLPFC_Br2743_post_manual_alignment",
    "DLPFC_Br3942_ant_manual_alignment",
    "DLPFC_Br3942_mid_manual_alignment",
    "DLPFC_Br3942_post_manual_alignment",
    "DLPFC_Br6423_ant_manual_alignment_extra_reads",
    "DLPFC_Br6423_mid_manual_alignment",
    "DLPFC_Br6423_post_extra_reads",
    "DLPFC_Br8492_ant_manual_alignment",
    "DLPFC_Br8492_mid_manual_alignment_extra_reads",
    "DLPFC_Br8492_post_manual_alignment",
    "Round2/DLPFC_Br2720_ant_manual_alignment",
    "Round2/DLPFC_Br2720_mid_manual_alignment",
    "Round2/DLPFC_Br2720_post_extra_reads",
    "Round2/DLPFC_Br6432_ant_manual_alignment",
    "Round2/DLPFC_Br6432_mid_manual_alignment",
    "Round2/DLPFC_Br6432_post_manual_alignment",
    "Round3/DLPFC_Br6471_ant_manual_alignment_all",
    "Round3/DLPFC_Br6471_mid_manual_alignment_all",
    "Round3/DLPFC_Br6471_post_manual_alignment_all",
    "Round3/DLPFC_Br6522_ant_manual_alignment_all",
    "Round3/DLPFC_Br6522_mid_manual_alignment_all",
    "Round3/DLPFC_Br6522_post_manual_alignment_all",
    "Round3/DLPFC_Br8325_ant_manual_alignment_all",
    "Round3/DLPFC_Br8325_mid_manual_alignment_all",
    "Round3/DLPFC_Br8325_post_manual_alignment_all",
    "Round3/DLPFC_Br8667_ant_extra_reads",
    "Round3/DLPFC_Br8667_mid_manual_alignment_all",
    "Round3/DLPFC_Br8667_post_manual_alignment_all",
    "Round4/DLPFC_Br2720_ant_2",
    "Round4/DLPFC_Br6432_ant_2",
    "Round4/DLPFC_Br8325_mid_2"
  ),
  subjects = c(rep(
    c("Br2743", "Br3942", "Br6423", "Br8492", "Br2720","Br6432","Br6471","Br6522","Br8325","Br8667"),
    each = 3
  ),"Br2720","Br6432","Br8325"),
  regions = c(rep(
    c("anterior", "middle", "posterior"),
    10
  ),"anterior","anterior","middle"),
  sex = c(rep(
    c("M", "M", "M", "F","F","M","M","M","F","F"),
    each = 3
  ),"F","M","F"),
  age = c(rep(
    c(61.54, 47.53, 51.73, 53.40,48.22,48.88,55.46,33.39,57.62,37.33),
    each = 3
  ),48.22,48.88,57.62),
  diagnosis = "Control"
)
sample_info$sample_path <- file.path(
  here::here("processed-data", "NextSeq"),
  sample_info$sample_id,
  "outs"
)
stopifnot(all(file.exists(sample_info$sample_path)))

#clean up sample_id
sample_info$sample_id <- gsub("_all|_extra_reads|DLPFC_|_manual_alignment", "", basename(sample_info$sample_id))

## Build SPE object
Sys.time()
spe_wrapper <- read10xVisiumWrapper(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres", "hires", "detected", "aligned"),
  load = TRUE,
  reference_gtf = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
)
Sys.time()

# 2021-11-23 15:05:05 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2021-11-23 15:12:23 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2021-11-23 15:13:09 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2021-11-23 15:13:15 rtracklayer::import: reading the reference GTF file
# 2021-11-23 15:14:03 adding gene information to the SPE object
# 2021-11-23 15:14:04 adding information used by spatialLIBD

## Add the experimental information
### Equivalent to:
## colData(spe_wrapper)$subject <- ...
spe_wrapper$subject <- sample_info$subjects[match(spe_wrapper$sample_id, sample_info$sample_id)]
spe_wrapper$region <- sample_info$regions[match(spe_wrapper$sample_id, sample_info$sample_id)]
spe_wrapper$sex <- sample_info$sex[match(spe_wrapper$sample_id, sample_info$sample_id)]
spe_wrapper$age <- sample_info$age[match(spe_wrapper$sample_id, sample_info$sample_id)]
spe_wrapper$diagnosis <- sample_info$diagnosis[match(spe_wrapper$sample_id, sample_info$sample_id)]

## Remove genes with no data
no_expr <- which(rowSums(counts(spe_wrapper)) == 0)
length(no_expr)
# [1] 7468
length(no_expr) / nrow(spe_wrapper) * 100
# [1] 20.40381

#remove genes that are not expressed in any spots
spe_wrapper <- spe_wrapper[-no_expr, ]

spe_wrapper <- spe_wrapper[, which(spatialData(spe_wrapper)$in_tissue=="TRUE")]

## Remove spots without counts
spe_wrapper <- spe_wrapper[, -which(colSums(counts(spe_wrapper)) == 0)]

## Read in cell counts and segmentation results
segmentations_list <- lapply(sample_info$sample_id, function(sampleid) {
  current<-sample_info$sample_path[sample_info$sample_id==sampleid]
  file <- file.path(current, "spatial", "tissue_spot_counts.csv")
  if(!file.exists(file)) return(NULL)
  x <- read.csv(file)
  x$key <- paste0(x$barcode, "_", sampleid)
  return(x)
})
## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <- Reduce(function(...) merge(..., all = TRUE), segmentations_list[lengths(segmentations_list) > 0])

## Add the information
segmentation_match <- match(spe_wrapper$key, segmentations$key)
segmentation_info <- segmentations[segmentation_match, - which(colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")), drop=FALSE]
colData(spe_wrapper) <- cbind(colData(spe_wrapper), segmentation_info)

load(here::here("processed-data","rdata", "spe", "spe_final.Rdata"),verbose = TRUE)
spe$cell_count <- NULL

#fixing the key since spe$Barcode doesn't exist anymore
spe$key <- paste0(colnames(spe), '_', spe$sample_id)

#checks that all spots are in the same order
stopifnot(identical(spe_wrapper$key,spe$key))

#fix rowRanges
rowData(spe)$gene_search <- paste0(rowData(spe)$gene_name, "; ", rowData(spe)$gene_id)

#checks that all the genes are in the same order
stopifnot(identical(rowRanges(spe_wrapper)$gene_search,rowRanges(spe)$gene_search))


##adding missing elements
logcounts(spe_wrapper) <- logcounts(spe)
rowRanges(spe_wrapper)$is.HVG <- rowRanges(spe)$is.HVG
stopifnot(identical(rowRanges(spe_wrapper),rowRanges(spe)))
reducedDims(spe_wrapper) <- c(reducedDims(spe),reducedDims(spe_wrapper))

spe$BayesSpace_cluster.init <- spe$cluster.init
spe$cluster.init <-NULL
spe$BayesSpace <- spe$spatial.cluster
spe$spatial.cluster <-NULL

colData(spe_wrapper) <- cbind(colData(spe_wrapper),colData(spe)[,!colnames(colData(spe))%in%colnames(colData(spe_wrapper))])
spe<-spe_wrapper
save(spe, file = here::here("processed-data", "rdata", "spe", "spe_merged_final.Rdata"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


# ─ Session info  ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# hash: fuel pump, flag: Botswana, person walking: dark skin tone
# 
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2021-11-23
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version  date (UTC) lib source
# AnnotationDbi            1.56.2   2021-11-09 [2] Bioconductor
# AnnotationHub            3.2.0    2021-10-26 [2] Bioconductor
# assertthat               0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# attempt                  0.3.1    2020-05-03 [1] CRAN (R 4.1.1)
# beachmat                 2.10.0   2021-10-26 [2] Bioconductor
# beeswarm                 0.4.0    2021-06-01 [1] CRAN (R 4.1.1)
# benchmarkme              1.0.7    2021-03-21 [1] CRAN (R 4.1.1)
# benchmarkmeData          1.0.4    2020-04-23 [1] CRAN (R 4.1.1)
# Biobase                * 2.54.0   2021-10-26 [2] Bioconductor
# BiocFileCache            2.2.0    2021-10-26 [2] Bioconductor
# BiocGenerics           * 0.40.0   2021-10-26 [2] Bioconductor
# BiocIO                   1.4.0    2021-10-26 [2] Bioconductor
# BiocManager              1.30.16  2021-06-15 [2] CRAN (R 4.1.2)
# BiocNeighbors            1.11.0   2021-05-19 [1] Bioconductor
# BiocParallel             1.28.1   2021-11-18 [2] Bioconductor
# BiocSingular             1.9.1    2021-06-08 [1] Bioconductor
# BiocVersion              3.14.0   2021-05-19 [2] Bioconductor
# Biostrings               2.62.0   2021-10-26 [2] Bioconductor
# bit                      4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
# bit64                    4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
# bitops                   1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# blob                     1.2.2    2021-07-23 [2] CRAN (R 4.1.0)
# bslib                    0.3.1    2021-10-06 [2] CRAN (R 4.1.2)
# cachem                   1.0.6    2021-08-19 [2] CRAN (R 4.1.2)
# callr                    3.7.0    2021-04-20 [2] CRAN (R 4.1.0)
# cli                      3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
# codetools                0.2-18   2020-11-04 [3] CRAN (R 4.1.2)
# colorout               * 1.2-2    2021-09-09 [1] Github (jalvesaq/colorout@79931fd)
# colorspace               2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
# config                   0.3.1    2020-12-17 [1] CRAN (R 4.1.1)
# cowplot                  1.1.1    2020-12-30 [1] CRAN (R 4.1.1)
# crayon                   1.4.2    2021-10-29 [2] CRAN (R 4.1.2)
# curl                     4.3.2    2021-06-23 [2] CRAN (R 4.1.0)
# data.table               1.14.2   2021-09-27 [2] CRAN (R 4.1.2)
# DBI                      1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
# dbplyr                   2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray             0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats       1.16.0   2021-10-26 [2] Bioconductor
# desc                     1.4.0    2021-09-28 [2] CRAN (R 4.1.2)
# digest                   0.6.28   2021-09-23 [2] CRAN (R 4.1.2)
# dockerfiler              0.1.4    2021-09-03 [1] CRAN (R 4.1.2)
# doParallel               1.0.16   2020-10-16 [2] CRAN (R 4.1.0)
# dotCall64                1.0-1    2021-02-11 [2] CRAN (R 4.1.0)
# dplyr                    1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
# dqrng                    0.3.0    2021-05-01 [1] CRAN (R 4.1.1)
# DropletUtils             1.13.2   2021-08-12 [1] Bioconductor
# DT                       0.20     2021-11-15 [2] CRAN (R 4.1.2)
# edgeR                    3.36.0   2021-10-26 [2] Bioconductor
# ellipsis                 0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# ExperimentHub            2.2.0    2021-10-26 [2] Bioconductor
# fansi                    0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
# fastmap                  1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
# fields                   13.3     2021-10-30 [2] CRAN (R 4.1.2)
# filelock                 1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
# foreach                  1.5.1    2020-10-15 [2] CRAN (R 4.1.0)
# fs                       1.5.0    2020-07-31 [2] CRAN (R 4.1.0)
# generics                 0.1.1    2021-10-25 [2] CRAN (R 4.1.2)
# GenomeInfoDb           * 1.30.0   2021-10-26 [2] Bioconductor
# GenomeInfoDbData         1.2.7    2021-11-01 [2] Bioconductor
# GenomicAlignments        1.30.0   2021-10-26 [2] Bioconductor
# GenomicRanges          * 1.46.1   2021-11-18 [2] Bioconductor
# ggbeeswarm               0.6.0    2017-08-07 [1] CRAN (R 4.1.1)
# ggplot2                  3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# ggrepel                  0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# glue                     1.5.0    2021-11-07 [2] CRAN (R 4.1.2)
# golem                    0.3.1    2021-04-17 [1] CRAN (R 4.1.1)
# gridExtra                2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                   0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# HDF5Array                1.22.1   2021-11-14 [2] Bioconductor
# here                   * 1.0.1    2020-12-13 [1] CRAN (R 4.1.1)
# htmltools                0.5.2    2021-08-25 [2] CRAN (R 4.1.2)
# htmlwidgets              1.5.4    2021-09-08 [2] CRAN (R 4.1.2)
# httpuv                   1.6.3    2021-09-09 [2] CRAN (R 4.1.2)
# httr                     1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# interactiveDisplayBase   1.32.0   2021-10-26 [2] Bioconductor
# IRanges                * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                    2.3.3    2019-02-05 [2] CRAN (R 4.1.0)
# iterators                1.0.13   2020-10-15 [2] CRAN (R 4.1.0)
# jquerylib                0.1.4    2021-04-26 [2] CRAN (R 4.1.0)
# jsonlite                 1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
# KEGGREST                 1.34.0   2021-10-26 [2] Bioconductor
# knitr                    1.36     2021-09-29 [2] CRAN (R 4.1.2)
# later                    1.3.0    2021-08-18 [2] CRAN (R 4.1.2)
# lattice                  0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lazyeval                 0.2.2    2019-03-15 [2] CRAN (R 4.1.0)
# lifecycle                1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                    3.50.0   2021-10-26 [2] Bioconductor
# locfit                   1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
# magick                   2.7.3    2021-08-18 [2] CRAN (R 4.1.2)
# magrittr                 2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# maps                     3.4.0    2021-09-25 [2] CRAN (R 4.1.2)
# Matrix                   1.3-4    2021-06-01 [3] CRAN (R 4.1.2)
# MatrixGenerics         * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats            * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# memoise                  2.0.0    2021-01-26 [2] CRAN (R 4.1.0)
# mime                     0.12     2021-09-28 [2] CRAN (R 4.1.2)
# munsell                  0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                   1.6.4    2021-10-18 [2] CRAN (R 4.1.2)
# pkgbuild                 1.2.0    2020-12-15 [2] CRAN (R 4.1.0)
# pkgconfig                2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# pkgload                  1.2.3    2021-10-13 [2] CRAN (R 4.1.2)
# plotly                   4.10.0   2021-10-09 [2] CRAN (R 4.1.2)
# png                      0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
# Polychrome               1.3.1    2021-07-16 [1] CRAN (R 4.1.1)
# prettyunits              1.1.1    2020-01-24 [2] CRAN (R 4.1.0)
# processx                 3.5.2    2021-04-30 [2] CRAN (R 4.1.0)
# promises                 1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)
# ps                       1.6.0    2021-02-28 [2] CRAN (R 4.1.0)
# purrr                    0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3              1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
# R.oo                     1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
# R.utils                  2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
# R6                       2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rappdirs                 0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
# RColorBrewer             1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                     1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
# RCurl                    1.98-1.5 2021-09-17 [2] CRAN (R 4.1.2)
# remotes                  2.4.1    2021-09-29 [2] CRAN (R 4.1.2)
# restfulr                 0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
# rhdf5                    2.38.0   2021-10-26 [2] Bioconductor
# rhdf5filters             1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib                 1.16.0   2021-10-26 [2] Bioconductor
# rjson                    0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
# rlang                    0.4.12   2021-10-18 [2] CRAN (R 4.1.2)
# roxygen2                 7.1.2    2021-09-08 [2] CRAN (R 4.1.2)
# rprojroot                2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# Rsamtools                2.10.0   2021-10-26 [2] Bioconductor
# RSQLite                  2.2.8    2021-08-21 [2] CRAN (R 4.1.2)
# rstudioapi               0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rsvd                     1.0.5    2021-04-16 [1] CRAN (R 4.1.1)
# rtracklayer              1.54.0   2021-10-26 [2] Bioconductor
# S4Vectors              * 0.32.3   2021-11-21 [2] Bioconductor
# sass                     0.4.0    2021-05-12 [2] CRAN (R 4.1.0)
# ScaledMatrix             1.1.0    2021-05-19 [1] Bioconductor
# scales                   1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scater                   1.21.3   2021-08-01 [1] Bioconductor
# scatterplot3d            0.3-41   2018-03-14 [1] CRAN (R 4.1.1)
# scuttle                  1.3.1    2021-08-05 [1] Bioconductor
# sessioninfo            * 1.2.1    2021-11-02 [2] CRAN (R 4.1.2)
# shiny                    1.7.1    2021-10-02 [2] CRAN (R 4.1.2)
# shinyWidgets             0.6.2    2021-09-17 [1] CRAN (R 4.1.2)
# SingleCellExperiment   * 1.16.0   2021-10-26 [2] Bioconductor
# spam                     2.7-0    2021-06-25 [2] CRAN (R 4.1.0)
# sparseMatrixStats        1.6.0    2021-10-26 [2] Bioconductor
# SpatialExperiment      * 1.3.4    2021-08-24 [1] Bioconductor
# spatialLIBD            * 1.7.3    2021-11-23 [1] Github (LieberInstitute/spatialLIBD@2516ba6)
# stringi                  1.7.5    2021-10-04 [2] CRAN (R 4.1.2)
# stringr                  1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment   * 1.24.0   2021-10-26 [2] Bioconductor
# testthat                 3.1.0    2021-10-04 [2] CRAN (R 4.1.2)
# tibble                   3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyr                    1.1.4    2021-09-27 [2] CRAN (R 4.1.2)
# tidyselect               1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# usethis                  2.1.3    2021-10-27 [2] CRAN (R 4.1.2)
# utf8                     1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                    0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# vipor                    0.4.5    2017-03-22 [1] CRAN (R 4.1.1)
# viridis                  0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite              0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                    2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
# xfun                     0.28     2021-11-04 [2] CRAN (R 4.1.2)
# XML                      3.99-0.8 2021-09-17 [2] CRAN (R 4.1.2)
# xml2                     1.3.2    2020-04-23 [2] CRAN (R 4.1.0)
# xtable                   1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
# XVector                  0.34.0   2021-10-26 [2] Bioconductor
# yaml                     2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
# zlibbioc                 1.40.0   2021-10-26 [2] Bioconductor
# 
# [1] /users/aspangle/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
