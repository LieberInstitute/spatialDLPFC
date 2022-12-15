library("spatialLIBD")
library("here")
library("sessioninfo")

## Load check functions
source(here("code", "analysis", "01_build_spe", "check_pxl_functions.R"), echo = TRUE, max.deparse.length = 500)

# $ ls -lht
# total 11G
# -rw-rw---- 1 lcollado lieber_lcolladotor 540M Dec 14 20:08 spe_subset_for_spatialLIBD.rds
# -rw-rw---- 1 lcollado lieber_lcolladotor 1.4G Dec 14 20:05 spe_filtered_final_with_clusters_and_deconvolution_results.rds
# -rw-rw---- 1 aspangle lieber_lcolladotor 1.4G May 12  2022 spe_filtered_final_with_clusters.Rdata
# -rw-rw---- 1 aspangle lieber_lcolladotor 1.4G May  6  2022 spe_filtered_final.Rdata
# -rw-rw---- 1 aspangle lieber_lcolladotor 141K May  3  2022 top.hvgs_all.Rdata
# -rw-rw---- 1 aspangle lieber_lcolladotor 890M May  3  2022 spe_final.Rdata
# -rw-rw---- 1 aspangle lieber_lcolladotor 909M May  3  2022 spe_raw_final.Rdata
# -rw-rw---- 1 aspangle lieber_lcolladotor 1.4G Mar  4  2022 spe_filtered_final_old.Rdata
# -rw-rw---- 1 aspangle lieber_lcolladotor 141K Mar  3  2022 top.hvgs_all_old.Rdata
# -rw-rw---- 1 aspangle lieber_lcolladotor 887M Mar  3  2022 spe_final_filtered.Rdata
# -rw-rw---- 1 aspangle lieber_lcolladotor 890M Mar  3  2022 spe_final_old.Rdata
# -rw-rw---- 1 aspangle lieber_lcolladotor 909M Mar  3  2022 spe_raw_final_old.Rdata

## Re-read the data
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
        "DLPFC_Br2720_ant_2",
        "DLPFC_Br2720_mid_manual_alignment",
        "DLPFC_Br2720_post_extra_reads",
        "DLPFC_Br6432_ant_2",
        "DLPFC_Br6432_mid_manual_alignment",
        "DLPFC_Br6432_post_manual_alignment",
        "DLPFC_Br6471_ant_manual_alignment_all",
        "DLPFC_Br6471_mid_manual_alignment_all",
        "DLPFC_Br6471_post_manual_alignment_all",
        "DLPFC_Br6522_ant_manual_alignment_all",
        "DLPFC_Br6522_mid_manual_alignment_all",
        "DLPFC_Br6522_post_manual_alignment_all",
        "DLPFC_Br8325_ant_manual_alignment_all",
        "DLPFC_Br8325_mid_2",
        "DLPFC_Br8325_post_manual_alignment_all",
        "DLPFC_Br8667_ant_extra_reads",
        "DLPFC_Br8667_mid_manual_alignment_all",
        "DLPFC_Br8667_post_manual_alignment_all"
    ),
    subjects = c(rep(
        c("Br2743", "Br3942", "Br6423", "Br8492", "Br2720", "Br6432", "Br6471", "Br6522", "Br8325", "Br8667"),
        each = 3
    )),
    regions = c(rep(
        c("anterior", "middle", "posterior"),
        10
    )),
    sex = c(rep(
        c("M", "M", "M", "F", "F", "M", "M", "M", "F", "F"),
        each = 3
    )),
    age = c(rep(
        c(61.54, 47.53, 51.73, 53.40, 48.22, 48.88, 55.46, 33.39, 57.62, 37.33),
        each = 3
    )),
    diagnosis = "Control"
)
sample_info$sample_path <- file.path(
    here::here("processed-data", "rerun_spaceranger"),
    sample_info$sample_id,
    "outs"
)
stopifnot(all(file.exists(sample_info$sample_path)))

# clean up sample_id
sample_info$sample_id <- gsub("_all|_extra_reads|DLPFC_|_manual_alignment", "", basename(sample_info$sample_id))


## Build spe_reloaded object
Sys.time()
spe_reloaded <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    reference_gtf = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
)
Sys.time()

# [1] "2022-03-01 11:22:49 EST"

## Add the experimental information
spe_reloaded$key <- paste0(colnames(spe_reloaded), "_", spe_reloaded$sample_id)
spe_reloaded$subject <- sample_info$subjects[match(spe_reloaded$sample_id, sample_info$sample_id)]
spe_reloaded$region <- sample_info$regions[match(spe_reloaded$sample_id, sample_info$sample_id)]
spe_reloaded$sex <- sample_info$sex[match(spe_reloaded$sample_id, sample_info$sample_id)]
spe_reloaded$age <- sample_info$age[match(spe_reloaded$sample_id, sample_info$sample_id)]
spe_reloaded$diagnosis <- sample_info$diagnosis[match(spe_reloaded$sample_id, sample_info$sample_id)]
spe_reloaded$sample_id_complete <- spe_reloaded$sample_id
spe_reloaded$sample_id <- gsub("_2", "", spe_reloaded$sample_id)

## Add key
spe_reloaded <- add_key(spe_reloaded)

## Save for later use
Sys.time()
save(spe_reloaded, file = here::here("processed-data", "rdata", "spe", "01_build_spe", "spe_reloaded.Rdata"))
Sys.time()

## Read the latest version that I recently made
# spe_latest <- readRDS(
#     here(
#         "processed-data",
#         "rdata",
#         "spe",
#         "01_build_spe",
#         "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
#     )
# )


## Read the latest object most people are using
Sys.time()
load(
    here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters.Rdata"
    ),
    verbose = TRUE
)
Sys.time()
spe_filtered_final_with_clusters <- spe
rm(spe)

## Read the one SPE prior to it, which is the filtered version
Sys.time()
load(
    here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final.Rdata"
    ),
    verbose = TRUE
)
Sys.time()
spe_filtered_final <- spe
rm(spe)

## Read the version prior to filtering low library size spots
Sys.time()
load(
    here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_final.Rdata"
    ),
    verbose = TRUE
)
Sys.time()
spe_final <- spe
rm(spe)

## Read the raw version before filtering
Sys.time()
load(
    here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_raw_final.Rdata"
    ),
    verbose = TRUE
)
Sys.time()
spe_raw_recent <- spe_raw
rm(spe_raw)


## Read the raw version before filtering called "old"
Sys.time()
load(
    here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_raw_final_old.Rdata"
    ),
    verbose = TRUE
)
Sys.time()
spe_raw_old <- spe_raw
rm(spe_raw)

check_pxl(spe_filtered_final_with_clusters, spe_reloaded)
# $array_row
#
#   TRUE
# 113927
#
# $array_col
#
#   TRUE
# 113927
#
# $pxl_row_in_fullres_diff
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   -2095    -992     450    1258    2275    9472
#
# $pxl_col_in_fullres_diff
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   -4797   -3453   -2621   -2097    -842    3361
#
# $swap_row_col_diff
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  -19546   -3184    2033    2104    7393   27001
#
# $swap_col_row_diff
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  -22734   -8127   -2872   -2943    2243   17337
#
# $pxl_row_in_fullres_tab
#
#  FALSE   TRUE
# 113922      5
#
# $pxl_row_in_fullretab
#
#  FALSE   TRUE
# 113922      5
check_pxl(spe_filtered_final_with_clusters, spe_filtered_final)
# [1] "All good!"
check_pxl(spe_filtered_final_with_clusters, spe_final)
# [1] "All good!"
check_pxl(spe_filtered_final_with_clusters, spe_raw_recent)
# [1] "All good!"
check_pxl(spe_filtered_final_with_clusters, spe_raw_old)
# Subsetting to those spots that match in spe_source
#
#  FALSE   TRUE
# 103136  10791
# [1] "All good!"

## This is the same as check_pxl(spe_filtered_final_with_clusters, spe_reloaded)
# check_pxl(spe_filtered_final, spe_reloaded)


check_limits(spe_filtered_final_with_clusters, spe_reloaded)
#        y_min y_max x_min x_max
# source   452   468   124   388
# target   162   501    89   453
check_limits(spe_filtered_final_with_clusters, spe_raw_recent)
#        y_min y_max x_min x_max
# source   452   468   124   388
# target   135   468    59   416

pdf(
    here(
        "plots",
        "01_build_spe",
        "check_pxl_spe_filtered_with_clusters_vs_spe_reloaded.pdf"
    ),
    width = 14
)
check_plot(spe_filtered_final_with_clusters, spe_reloaded)
check_plot(spe_filtered_final_with_clusters, spe_reloaded, auto_crop = TRUE)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.2 Patched (2022-12-14 r83465)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2022-12-15
#  pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version   date (UTC) lib source
#  AnnotationDbi            1.60.0    2022-11-01 [2] Bioconductor
#  AnnotationHub            3.6.0     2022-11-01 [2] Bioconductor
#  assertthat               0.2.1     2019-03-21 [2] CRAN (R 4.2.1)
#  attempt                  0.3.1     2020-05-03 [2] CRAN (R 4.2.1)
#  beachmat                 2.14.0    2022-11-01 [2] Bioconductor
#  beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
#  benchmarkme              1.0.8     2022-06-12 [2] CRAN (R 4.2.1)
#  benchmarkmeData          1.0.4     2020-04-23 [2] CRAN (R 4.2.1)
#  Biobase                * 2.58.0    2022-11-01 [1] Bioconductor
#  BiocFileCache            2.6.0     2022-11-01 [1] Bioconductor
#  BiocGenerics           * 0.44.0    2022-11-01 [2] Bioconductor
#  BiocIO                   1.8.0     2022-11-01 [2] Bioconductor
#  BiocManager              1.30.19   2022-10-25 [2] CRAN (R 4.2.2)
#  BiocNeighbors            1.16.0    2022-11-01 [2] Bioconductor
#  BiocParallel             1.32.4    2022-12-01 [1] Bioconductor
#  BiocSingular             1.14.0    2022-11-01 [2] Bioconductor
#  BiocVersion              3.16.0    2022-04-26 [2] Bioconductor
#  Biostrings               2.66.0    2022-11-01 [2] Bioconductor
#  bit                      4.0.5     2022-11-15 [2] CRAN (R 4.2.2)
#  bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
#  bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
#  blob                     1.2.3     2022-04-10 [2] CRAN (R 4.2.1)
#  bslib                    0.4.1     2022-11-02 [2] CRAN (R 4.2.2)
#  cachem                   1.0.6     2021-08-19 [2] CRAN (R 4.2.1)
#  cli                      3.4.1     2022-09-23 [2] CRAN (R 4.2.1)
#  codetools                0.2-18    2020-11-04 [3] CRAN (R 4.2.2)
#  colorout                 1.2-2     2022-11-02 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace               2.0-3     2022-02-21 [2] CRAN (R 4.2.1)
#  config                   0.3.1     2020-12-17 [2] CRAN (R 4.2.1)
#  cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
#  crayon                   1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
#  curl                     4.3.3     2022-10-06 [2] CRAN (R 4.2.1)
#  data.table               1.14.6    2022-11-16 [2] CRAN (R 4.2.2)
#  DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.2.1)
#  dbplyr                   2.2.1     2022-06-27 [2] CRAN (R 4.2.1)
#  DelayedArray             0.24.0    2022-11-01 [2] Bioconductor
#  DelayedMatrixStats       1.20.0    2022-11-01 [1] Bioconductor
#  desc                     1.4.2     2022-09-08 [2] CRAN (R 4.2.1)
#  digest                   0.6.31    2022-12-11 [2] CRAN (R 4.2.2)
#  doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
#  dotCall64                1.0-2     2022-10-03 [2] CRAN (R 4.2.1)
#  dplyr                    1.0.10    2022-09-01 [2] CRAN (R 4.2.1)
#  dqrng                    0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
#  DropletUtils             1.18.1    2022-11-22 [2] Bioconductor
#  DT                       0.26      2022-10-19 [2] CRAN (R 4.2.1)
#  edgeR                    3.40.0    2022-11-01 [2] Bioconductor
#  ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.2.1)
#  ExperimentHub            2.6.0     2022-11-01 [2] Bioconductor
#  fansi                    1.0.3     2022-03-24 [2] CRAN (R 4.2.1)
#  farver                   2.1.1     2022-07-06 [2] CRAN (R 4.2.1)
#  fastmap                  1.1.0     2021-01-25 [2] CRAN (R 4.2.1)
#  fields                   14.1      2022-08-12 [2] CRAN (R 4.2.1)
#  filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
#  foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
#  fs                       1.5.2     2021-12-08 [2] CRAN (R 4.2.1)
#  generics                 0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
#  GenomeInfoDb           * 1.34.4    2022-12-01 [1] Bioconductor
#  GenomeInfoDbData         1.2.9     2022-09-29 [2] Bioconductor
#  GenomicAlignments        1.34.0    2022-11-01 [2] Bioconductor
#  GenomicRanges          * 1.50.1    2022-11-06 [2] Bioconductor
#  ggbeeswarm               0.6.0     2017-08-07 [2] CRAN (R 4.2.1)
#  ggplot2                  3.4.0     2022-11-04 [2] CRAN (R 4.2.2)
#  ggrepel                  0.9.2     2022-11-06 [2] CRAN (R 4.2.2)
#  glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
#  golem                    0.3.5     2022-10-18 [2] CRAN (R 4.2.1)
#  gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
#  gtable                   0.3.1     2022-09-01 [2] CRAN (R 4.2.1)
#  HDF5Array                1.26.0    2022-11-01 [2] Bioconductor
#  here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
#  htmltools                0.5.4     2022-12-07 [2] CRAN (R 4.2.2)
#  htmlwidgets              1.5.4     2021-09-08 [2] CRAN (R 4.2.1)
#  httpuv                   1.6.7     2022-12-14 [2] CRAN (R 4.2.2)
#  httr                     1.4.4     2022-08-17 [2] CRAN (R 4.2.1)
#  interactiveDisplayBase   1.36.0    2022-11-01 [2] Bioconductor
#  IRanges                * 2.32.0    2022-11-01 [2] Bioconductor
#  irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
#  iterators                1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
#  jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.2.1)
#  jsonlite                 1.8.4     2022-12-06 [2] CRAN (R 4.2.2)
#  KEGGREST                 1.38.0    2022-11-01 [2] Bioconductor
#  knitr                    1.41      2022-11-18 [2] CRAN (R 4.2.2)
#  labeling                 0.4.2     2020-10-20 [2] CRAN (R 4.2.1)
#  later                    1.3.0     2021-08-18 [2] CRAN (R 4.2.1)
#  lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.2.2)
#  lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
#  lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
#  limma                    3.54.0    2022-11-01 [1] Bioconductor
#  locfit                   1.5-9.6   2022-07-11 [2] CRAN (R 4.2.1)
#  magick                   2.7.3     2021-08-18 [2] CRAN (R 4.2.1)
#  magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
#  maps                     3.4.1     2022-10-30 [2] CRAN (R 4.2.2)
#  Matrix                   1.5-3     2022-11-11 [2] CRAN (R 4.2.2)
#  MatrixGenerics         * 1.10.0    2022-11-01 [1] Bioconductor
#  matrixStats            * 0.63.0    2022-11-18 [2] CRAN (R 4.2.2)
#  memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
#  mime                     0.12      2021-09-28 [2] CRAN (R 4.2.1)
#  munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
#  paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.2.1)
#  pillar                   1.8.1     2022-08-19 [2] CRAN (R 4.2.1)
#  pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
#  pkgload                  1.3.2     2022-11-16 [2] CRAN (R 4.2.2)
#  plotly                   4.10.1    2022-11-07 [2] CRAN (R 4.2.2)
#  png                      0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
#  promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.2.1)
#  purrr                    0.3.5     2022-10-06 [2] CRAN (R 4.2.1)
#  R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
#  R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
#  R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
#  R6                       2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
#  rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
#  RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
#  Rcpp                     1.0.9     2022-07-08 [2] CRAN (R 4.2.1)
#  RCurl                    1.98-1.9  2022-10-03 [2] CRAN (R 4.2.1)
#  rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.2.1)
#  restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
#  rhdf5                    2.42.0    2022-11-01 [2] Bioconductor
#  rhdf5filters             1.10.0    2022-11-01 [2] Bioconductor
#  Rhdf5lib                 1.20.0    2022-11-01 [2] Bioconductor
#  rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
#  rlang                    1.0.6     2022-09-24 [2] CRAN (R 4.2.1)
#  rmote                    0.3.4     2022-11-02 [1] Github (cloudyr/rmote@fbce611)
#  roxygen2                 7.2.3     2022-12-08 [2] CRAN (R 4.2.2)
#  rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
#  Rsamtools                2.14.0    2022-11-01 [2] Bioconductor
#  RSQLite                  2.2.19    2022-11-24 [2] CRAN (R 4.2.2)
#  rstudioapi               0.14      2022-08-22 [2] CRAN (R 4.2.1)
#  rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
#  rtracklayer              1.58.0    2022-11-01 [2] Bioconductor
#  S4Vectors              * 0.36.1    2022-12-05 [1] Bioconductor
#  sass                     0.4.4     2022-11-24 [2] CRAN (R 4.2.2)
#  ScaledMatrix             1.6.0     2022-11-01 [1] Bioconductor
#  scales                   1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
#  scater                   1.26.1    2022-11-13 [2] Bioconductor
#  scuttle                  1.8.2     2022-12-07 [2] Bioconductor
#  servr                    0.25      2022-11-04 [1] CRAN (R 4.2.2)
#  sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
#  shiny                    1.7.3     2022-10-25 [2] CRAN (R 4.2.2)
#  shinyWidgets             0.7.5     2022-11-17 [2] CRAN (R 4.2.2)
#  SingleCellExperiment   * 1.20.0    2022-11-01 [1] Bioconductor
#  spam                     2.9-1     2022-08-07 [2] CRAN (R 4.2.1)
#  sparseMatrixStats        1.10.0    2022-11-01 [2] Bioconductor
#  SpatialExperiment      * 1.8.0     2022-11-01 [2] Bioconductor
#  spatialLIBD            * 1.11.3    2022-12-15 [1] Github (LieberInstitute/spatialLIBD@362bfb6)
#  statmod                  1.4.37    2022-08-12 [2] CRAN (R 4.2.1)
#  stringi                  1.7.8     2022-07-11 [2] CRAN (R 4.2.1)
#  stringr                  1.5.0     2022-12-02 [2] CRAN (R 4.2.2)
#  SummarizedExperiment   * 1.28.0    2022-11-01 [2] Bioconductor
#  tibble                   3.1.8     2022-07-22 [2] CRAN (R 4.2.1)
#  tidyr                    1.2.1     2022-09-08 [2] CRAN (R 4.2.1)
#  tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
#  usethis                  2.1.6     2022-05-25 [2] CRAN (R 4.2.1)
#  utf8                     1.2.2     2021-07-24 [2] CRAN (R 4.2.1)
#  vctrs                    0.5.1     2022-11-16 [2] CRAN (R 4.2.2)
#  vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
#  viridis                  0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
#  viridisLite              0.4.1     2022-08-22 [2] CRAN (R 4.2.1)
#  withr                    2.5.0     2022-03-03 [2] CRAN (R 4.2.1)
#  xfun                     0.35      2022-11-16 [2] CRAN (R 4.2.2)
#  XML                      3.99-0.13 2022-12-04 [2] CRAN (R 4.2.2)
#  xml2                     1.3.3     2021-11-30 [2] CRAN (R 4.2.1)
#  xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.2.1)
#  XVector                  0.38.0    2022-11-01 [2] Bioconductor
#  yaml                     2.3.6     2022-10-18 [2] CRAN (R 4.2.1)
#  zlibbioc                 1.44.0    2022-11-01 [1] Bioconductor
#
#  [1] /users/lcollado/R/4.2.x
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
