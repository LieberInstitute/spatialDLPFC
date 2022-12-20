library("here")
library("sessioninfo")
library("spatialLIBD")
library("Polychrome")

# load spe object and clusters
load(
    file = here::here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters.Rdata"
    ),
    verbose = TRUE
)

## Set k
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## Set colors
mycolors <- Polychrome::palette36.colors(k)

## Note that some k's have missing clusters, like k = 28 is missing cluster 18 and 21
k_present <-
    sort(unique(colData(spe)[[paste0("bayesSpace_harmony_", k)]]))
mycolors <- mycolors[k_present]
names(mycolors) <- k_present

## Choose variables of interest
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
sample_ids <- unique(colData(spe)$sample_id)

pdf(file = here::here(
    "plots",
    "03_BayesSpace",
    paste0("polychrome_vis_clus_bayesSpace_harmony_square", k, ".pdf")
))
for (i in seq_along(sample_ids)) {
    p <- vis_clus(
        spe = spe,
        clustervar = bayesSpace_name,
        sampleid = sample_ids[i],
        colors = mycolors
    )
    print(p)
}
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.2 Patched (2022-12-14 r83479)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2022-12-20
#  pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/pandoc
#
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
#  bslib                    0.4.2     2022-12-16 [2] CRAN (R 4.2.2)
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
#  edgeR                    3.40.1    2022-12-14 [2] Bioconductor
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
#  GenomicRanges          * 1.50.2    2022-12-16 [2] Bioconductor
#  ggbeeswarm               0.7.1     2022-12-16 [2] CRAN (R 4.2.2)
#  ggplot2                  3.4.0     2022-11-04 [2] CRAN (R 4.2.2)
#  ggrepel                  0.9.2     2022-11-06 [2] CRAN (R 4.2.2)
#  glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
#  golem                    0.3.5     2022-10-18 [2] CRAN (R 4.2.1)
#  gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
#  gtable                   0.3.1     2022-09-01 [2] CRAN (R 4.2.1)
#  HDF5Array                1.26.0    2022-11-01 [2] Bioconductor
#  here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
#  htmltools                0.5.4     2022-12-07 [2] CRAN (R 4.2.2)
#  htmlwidgets              1.6.0     2022-12-15 [2] CRAN (R 4.2.2)
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
#  Polychrome             * 1.5.1     2022-05-03 [1] CRAN (R 4.2.2)
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
#  scatterplot3d            0.3-42    2022-09-08 [1] CRAN (R 4.2.2)
#  scuttle                  1.8.3     2022-12-14 [2] Bioconductor
#  servr                    0.25      2022-11-04 [1] CRAN (R 4.2.2)
#  sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
#  shiny                    1.7.4     2022-12-15 [2] CRAN (R 4.2.2)
#  shinyWidgets             0.7.5     2022-11-17 [2] CRAN (R 4.2.2)
#  SingleCellExperiment   * 1.20.0    2022-11-01 [1] Bioconductor
#  spam                     2.9-1     2022-08-07 [2] CRAN (R 4.2.1)
#  sparseMatrixStats        1.10.0    2022-11-01 [2] Bioconductor
#  SpatialExperiment      * 1.8.0     2022-11-01 [2] Bioconductor
#  spatialLIBD            * 1.11.4    2022-12-17 [1] Github (LieberInstitute/spatialLIBD@1aecde8)
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
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
