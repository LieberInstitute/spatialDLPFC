library('basilisk')
library('scRNAseq')
library('SingleCellExperiment')
library('zellkonverter')
library('data.table')
library('tidyverse')
library("DeconvoBuddies")
library("here")
library("jaffelab")

#  Path to write the python AnnData object
dir.create(file.path(here::here(), "tangram_libd/out/"), showWarnings = FALSE)
visium_out = file.path(here::here(), "tangram_libd/out/visium_dlpfc.h5ad")
sc_out = file.path(here::here(), "tangram_libd/out/sce_dlpfc.h5ad")

#  An example SingleCellExperiment object
load(file.path(here::here(), "tangram_libd/data/SCE_DLPFC_tran-etal.rda"))
load(file.path(here::here(), "tangram_libd/data/sce_combined.rda"))

# read in brain sample names
brain_samples = readLines(file.path(here::here(), "tangram_libd/data/brain_samples.txt"))

# filtering to only brain samples of interest in the spatial data
sce <- sce[, sce$sample_name %in% brain_samples]

rna.sce <- sce.dlpfc
spatial.seq <- sce

# prepairing for expression_cutoff
seed <- 20210324

# Find expression cutoff for each dataset and filter out genes below that
# cutoff
counts_rna = as.matrix(assays(rna.sce)$counts)
cutoff_rna <- max(expression_cutoff(counts_rna, seed = seed))
rna.sce = rna.sce[rowMeans(counts_rna) > cutoff_rna, ]

counts_spat = as.matrix(assays(spatial.seq)$counts)
cutoff_spat <- max(expression_cutoff(counts_spat, seed = seed))
spatial.seq = spatial.seq[rowMeans(counts_spat) > cutoff_spat, ]

# Getting overlaps between spatial and scRNAseq data
overlaps <-
    rna.sce[rowData(rna.sce)$Symbol %in% rowData(spatial.seq)$gene_name,]

mean_ratio <-
    map(
        "cell_type",
        ~ get_mean_ratio2(
            overlaps,
            cellType_col = .x,
            assay_name = "counts",
            add_symbol = TRUE
        )
    )
markers_1vAll <-
    map(
        "cell_type",
        ~ findMarkers_1vAll(
            overlaps,
            cellType_col = .x,
            assay_name = "counts",
            add_symbol = TRUE
        )
    )

marker_stats <- map2(mean_ratio, markers_1vAll,
                     ~ left_join(.x, .y, by = c("gene", "cellType.target", "Symbol")))

tangram_markers <-
    marker_stats %>% as.data.table() %>% group_by("cellType") %>% slice_head(n = 1000) %>% ungroup() %>% select("Symbol")
tangram_markers

write.csv(
    tangram_markers,
    file = file.path(here::here(), "tangram_libd/data/marker_stats.csv"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

pdf(file.path(here::here(), "tangram_libd/out/marker_stats.pdf"))
#### Plot ####
ratio_plot <- map(
    marker_stats,
    ~ .x %>%
        ggplot(aes(ratio, std.logFC)) +
        geom_point(size = 0.5) +
        facet_wrap( ~ cellType.target, scales = "free_x") +
        labs(x = "mean(target logcount)/mean(highest non-target logcount)") +
        theme_bw() +
        NULL
)

ratio_plot

dev.off()

rm(sce, sce.dlpfc)

###############################################################################
#  The main code we'll use in general to convert SCE R objects to AnnData
#  python objects, as a preprocessing step to running tangram on our own data
###############################################################################

write_anndata = function(sce, out_path) {
    invisible(
        basiliskRun(
            fun = function(sce, filename) {
                library('zellkonverter')
                library('reticulate')
                
                # Convert SCE to AnnData:
                adata <- SCE2AnnData(sce)
                
                #  Write AnnData object to disk
                adata$write(filename = filename)
                
                return()
            },
            env = zellkonverter:::anndata_env,
            sce = sce,
            filename = out_path
        )
    )
}

write_anndata(rna.sce, sc_out)
write_anndata(spatial.seq, visium_out)

# sessioninfo::session_info()
# ─ Session info ───────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.0.4 RC (2021-02-08 r79975)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2021-05-06
#
# ─ Packages ───────────────────────────────────────────────────────────────────
# package                * version  date       lib
# AnnotationDbi            1.52.0   2020-10-27 [2]
# AnnotationFilter         1.14.0   2020-10-27 [2]
# AnnotationHub            2.22.1   2021-04-16 [2]
# askpass                  1.1      2019-01-13 [2]
# assertthat               0.2.1    2019-03-21 [2]
# backports                1.2.1    2020-12-09 [2]
# basilisk               * 1.2.1    2020-12-16 [1]
# basilisk.utils           1.2.2    2021-01-27 [1]
# beachmat                 2.6.4    2020-12-20 [2]
# Biobase                * 2.50.0   2020-10-27 [2]
# BiocFileCache            1.14.0   2020-10-27 [2]
# BiocGenerics           * 0.36.1   2021-04-16 [2]
# BiocManager              1.30.12  2021-03-28 [2]
# BiocNeighbors            1.8.2    2020-12-07 [1]
# BiocParallel             1.24.1   2020-11-06 [2]
# BiocSingular             1.6.0    2020-10-27 [1]
# BiocVersion              3.12.0   2020-04-27 [2]
# biomaRt                  2.46.3   2021-02-09 [2]
# Biostrings               2.58.0   2020-10-27 [2]
# bit                      4.0.4    2020-08-04 [2]
# bit64                    4.0.5    2020-08-30 [2]
# bitops                   1.0-7    2021-04-24 [2]
# blob                     1.2.1    2020-01-20 [2]
# bluster                  1.0.0    2020-10-27 [1]
# broom                    0.7.6    2021-04-05 [2]
# cachem                   1.0.4    2021-02-13 [2]
# cellranger               1.1.0    2016-07-27 [2]
# cli                      2.5.0    2021-04-26 [2]
# colorspace               2.0-0    2020-11-11 [2]
# crayon                   1.4.1    2021-02-08 [2]
# curl                     4.3      2019-12-02 [2]
# data.table             * 1.14.0   2021-02-21 [2]
# DBI                      1.1.1    2021-01-15 [2]
# dbplyr                   2.1.1    2021-04-06 [2]
# DeconvoBuddies         * 0.99.0   2021-03-05 [1]
# DelayedArray             0.16.3   2021-03-24 [2]
# DelayedMatrixStats       1.12.3   2021-02-03 [2]
# digest                   0.6.27   2020-10-24 [2]
# dplyr                  * 1.0.5    2021-03-05 [1]
# dqrng                    0.3.0    2021-05-01 [1]
# edgeR                    3.32.1   2021-01-14 [2]
# ellipsis                 0.3.2    2021-04-29 [2]
# ensembldb                2.14.1   2021-04-19 [2]
# ExperimentHub            1.16.1   2021-04-16 [2]
# fansi                    0.4.2    2021-01-15 [2]
# farver                   2.1.0    2021-02-28 [2]
# fastmap                  1.1.0    2021-01-25 [2]
# filelock                 1.0.2    2018-10-05 [1]
# forcats                * 0.5.1    2021-01-27 [2]
# fs                       1.5.0    2020-07-31 [2]
# generics                 0.1.0    2020-10-31 [2]
# GenomeInfoDb           * 1.26.7   2021-04-08 [2]
# GenomeInfoDbData         1.2.4    2020-11-30 [2]
# GenomicAlignments        1.26.0   2020-10-27 [2]
# GenomicFeatures          1.42.3   2021-04-01 [2]
# GenomicRanges          * 1.42.0   2020-10-27 [2]
# ggplot2                * 3.3.3    2020-12-30 [2]
# glue                     1.4.2    2020-08-27 [2]
# googledrive              1.0.1    2020-05-05 [1]
# gtable                   0.3.0    2019-03-25 [2]
# haven                    2.4.1    2021-04-23 [2]
# here                   * 1.0.1    2020-12-13 [1]
# hms                      1.0.0    2021-01-13 [2]
# htmltools                0.5.1.1  2021-01-22 [2]
# httpuv                   1.6.0    2021-04-23 [2]
# httr                     1.4.2    2020-07-20 [2]
# igraph                   1.2.6    2020-10-06 [2]
# interactiveDisplayBase   1.28.0   2020-10-27 [2]
# IRanges                * 2.24.1   2020-12-12 [2]
# irlba                    2.3.3    2019-02-05 [2]
# jaffelab               * 0.99.30  2021-02-16 [1]
# jsonlite                 1.7.2    2020-12-09 [2]
# labeling                 0.4.2    2020-10-20 [2]
# later                    1.2.0    2021-04-23 [2]
# lattice                  0.20-41  2020-04-02 [3]
# lazyeval                 0.2.2    2019-03-15 [2]
# lifecycle                1.0.0    2021-02-15 [2]
# limma                    3.46.0   2020-10-27 [2]
# locfit                   1.5-9.4  2020-03-25 [2]
# lubridate                1.7.9.2  2020-11-13 [1]
# magrittr                 2.0.1    2020-11-17 [2]
# Matrix                   1.3-2    2021-01-06 [3]
# MatrixGenerics         * 1.2.1    2021-01-30 [2]
# matrixStats            * 0.58.0   2021-01-29 [2]
# memoise                  2.0.0    2021-01-26 [2]
# mime                     0.10     2021-02-13 [2]
# modelr                   0.1.8    2020-05-19 [1]
# munsell                  0.5.0    2018-06-12 [2]
# openssl                  1.4.3    2020-09-18 [2]
# pillar                   1.6.0    2021-04-13 [2]
# pkgconfig                2.0.3    2019-09-22 [2]
# png                      0.1-7    2013-12-03 [2]
# prettyunits              1.1.1    2020-01-24 [2]
# progress                 1.2.2    2019-05-16 [2]
# promises                 1.2.0.1  2021-02-11 [2]
# ProtGenerics             1.22.0   2020-10-27 [2]
# ps                       1.6.0    2021-02-28 [2]
# purrr                  * 0.3.4    2020-04-17 [2]
# R6                       2.5.0    2020-10-28 [2]
# rafalib                * 1.0.0    2015-08-09 [1]
# rappdirs                 0.3.3    2021-01-31 [2]
# RColorBrewer             1.1-2    2014-12-07 [2]
# Rcpp                     1.0.6    2021-01-15 [2]
# RCurl                    1.98-1.3 2021-03-16 [2]
# readr                  * 1.4.0    2020-10-05 [2]
# readxl                   1.3.1    2019-03-13 [2]
# reprex                   2.0.0    2021-04-02 [1]
# reticulate             * 1.20     2021-05-03 [1]
# rlang                    0.4.10   2020-12-30 [2]
# rprojroot                2.0.2    2020-11-15 [2]
# Rsamtools                2.6.0    2020-10-27 [2]
# RSQLite                  2.2.7    2021-04-22 [2]
# rstudioapi               0.13     2020-11-12 [2]
# rsvd                     1.0.5    2021-04-16 [1]
# rtracklayer              1.50.0   2020-10-27 [2]
# rvest                    1.0.0    2021-03-09 [2]
# S4Vectors              * 0.28.1   2020-12-09 [2]
# scales                   1.1.1    2020-05-11 [2]
# scran                    1.18.5   2021-02-04 [1]
# scRNAseq               * 2.4.0    2020-11-09 [1]
# scuttle                  1.0.4    2020-12-17 [1]
# segmented                1.3-4    2021-04-22 [1]
# sessioninfo              1.1.1    2018-11-05 [2]
# shiny                    1.6.0    2021-01-25 [2]
# SingleCellExperiment   * 1.12.0   2020-10-27 [1]
# sparseMatrixStats        1.2.1    2021-02-02 [2]
# statmod                  1.4.35   2020-10-19 [2]
# stringi                  1.5.3    2020-09-09 [2]
# stringr                * 1.4.0    2019-02-10 [2]
# SummarizedExperiment   * 1.20.0   2020-10-27 [2]
# tibble                 * 3.1.1    2021-04-18 [2]
# tidyr                  * 1.1.3    2021-03-03 [2]
# tidyselect               1.1.1    2021-04-30 [2]
# tidyverse              * 1.3.0    2019-11-21 [1]
# utf8                     1.2.1    2021-03-12 [2]
# vctrs                    0.3.8    2021-04-29 [2]
# withr                    2.4.2    2021-04-18 [2]
# XML                      3.99-0.6 2021-03-16 [2]
# xml2                     1.3.2    2020-04-23 [2]
# xtable                   1.8-4    2019-04-21 [2]
# XVector                  0.30.0   2020-10-27 [2]
# yaml                     2.2.1    2020-02-01 [2]
# zellkonverter          * 1.0.2    2021-01-28 [1]
# zlibbioc                 1.36.0   2020-10-27 [2]
# source
# Bioconductor
# Bioconductor
# Bioconductor
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# Bioconductor
# Bioconductor
# Bioconductor
# Bioconductor
# Bioconductor
# Bioconductor
# CRAN (R 4.0.4)
# Bioconductor
# Bioconductor
# Bioconductor
# Bioconductor
# Bioconductor
# Bioconductor
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# Bioconductor
# CRAN (R 4.0.4)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# Github (lahuuki/DeconvoBuddies@73e9e66)
# Bioconductor
# Bioconductor
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.4)
# Bioconductor
# CRAN (R 4.0.4)
# Bioconductor
# Bioconductor
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# Bioconductor
# Bioconductor
# Bioconductor
# Bioconductor
# Bioconductor
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# Bioconductor
# Bioconductor
# CRAN (R 4.0.3)
# Github (LieberInstitute/jaffelab@42637ff)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# Bioconductor
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# Bioconductor
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.2)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# Bioconductor
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# Bioconductor
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.4)
# Bioconductor
# CRAN (R 4.0.4)
# Bioconductor
# CRAN (R 4.0.3)
# Bioconductor
# Bioconductor
# Bioconductor
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# Bioconductor
# Bioconductor
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# Bioconductor
# CRAN (R 4.0.4)
# CRAN (R 4.0.4)
# CRAN (R 4.0.4)
# CRAN (R 4.0.2)
# CRAN (R 4.0.4)
# CRAN (R 4.0.4)
# CRAN (R 4.0.4)
# CRAN (R 4.0.4)
# CRAN (R 4.0.3)
# CRAN (R 4.0.3)
# Bioconductor
# CRAN (R 4.0.3)
# Bioconductor
# Bioconductor
#
# [1] /users/aseyedia/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
