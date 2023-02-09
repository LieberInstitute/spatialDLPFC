library("spatialLIBD")
library("here")
library("jaffelab")
library("sessioninfo")
library("tidyverse")

#   Paths for manual annotation of layers and wrinkles, which was done for
#   3 of the 30 samples
anno_samples <- c("Br6522_ant", "Br6522_mid", "Br8667_post")

anno_wrinkle_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF",
    "spatialLIBD_ManualAnnotation_{sample_id}_wrinkle.csv"
)

anno_layers_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF",
    "spatialLIBD_ManualAnnotation_{sample_id}_layers.csv"
)

## Takes about 4 min
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

## Check the size in GB
Sys.time()
lobstr::obj_size(spe) ## takes about 2 min
Sys.time()
# 6.92 GB

## Fix some variables
vars <- colnames(colData(spe))
colnames(colData(spe))[grep("^X10x", vars)] <-
    gsub("X10x", "SpaceRanger_10x", vars[grep("^X10x", vars)])
colnames(colData(spe))[grep("^bayes", vars)] <-
    gsub("bayes", "Bayes", vars[grep("^bayes", vars)])
colnames(colData(spe))[grep("^region$", vars)] <-
    gsub("region", "position", vars[grep("^region$", vars)])

## Fix numerical columns
vars <- colnames(colData(spe))
colnames(colData(spe))[grep("SpaceRanger_10x_kmeans_", vars)] <-
    paste0("SpaceRanger_10x_kmeans_", sprintf("%02d", as.integer(ss(
        ss(vars[grep("SpaceRanger_10x_kmeans_", vars)], "_kmeans_", 2), "_"
    ))))
colnames(colData(spe))[grep("_[0-9]$", vars)] <-
    paste0(ss(vars[grep("_[0-9]$", vars)], "_[0-9]", 1), "_", sprintf("%02d", as.integer(ss(vars[grep("_[0-9]$", vars)], "[y|a]_", 2))))

## Import spot deconvolution results
for (resolution in c("broad", "layer")) {
    for (deconvo in c("01-tangram", "03-cell2location", "04-spotlight")) {
        message(Sys.time(), " importing: ", resolution, " for ", deconvo)
        spe <- cluster_import(
            spe,
            here(
                "processed-data",
                "spot_deconvo",
                deconvo,
                "nonIF",
                resolution,
                "raw_results"
            ),
            prefix = paste0(resolution, "_", gsub(".+-", "", deconvo), "_")
        )
    }
}

## re-sort columns after fixing the numerical ones
colData(spe) <- colData(spe)[, sort(colnames(colData(spe)))]

## Check the size in GB
lobstr::obj_size(spe)
# 6.96 GB

## Delete original VistoSeg counts (prior to updates in VistoSeg)
spe$VistoSeg_count_deprecated <- spe$count
spe$count <- NULL

## Import updated VistoSeg cell counts
nonIF_id_path <-
    here("processed-data", "spot_deconvo", "nonIF_ID_table.csv")
nonIF_counts_path <- here(
    "processed-data",
    "rerun_spaceranger",
    "{sample_id}",
    "outs",
    "spatial",
    "tissue_spot_counts.csv"
)

## Import code from
## https://github.com/LieberInstitute/spatialDLPFC/blob/444ce5cb408b90aa59066379df7a0e30cf0d2447/code/spot_deconvo/05-shared_utilities/01-r_to_python.R#L130-L167
id_table <- read.csv(nonIF_id_path)

## Initialize the new variables
spe$VistoSeg_count <- NA
spe$VistoSeg_proportion <- NA

for (sample_id in unique(spe$sample_id)) {
    message(Sys.time(), " processing sample id ", sample_id)
    #   Correctly determine the path for the cell counts for this sample, then
    #   read in
    long_id <-
        id_table[match(sample_id, id_table$short_id), "long_id"]
    this_path <-
        sub("{sample_id}", long_id, nonIF_counts_path, fixed = TRUE)
    cell_counts <- read.csv(this_path)

    #   All spots in the object should have counts
    stopifnot(all(colnames(spe[, spe$sample_id == sample_id]) %in%
        cell_counts$barcode))

    #   Line up the rows of 'cell_counts' with the sample-subsetted SPE object
    cell_counts <- cell_counts[match(
        colnames(spe[, spe$sample_id == sample_id]),
        cell_counts$barcode
    ), ]

    #   Add this sample's counts to the SPE object
    spe$VistoSeg_count[spe$sample_id == sample_id] <-
        cell_counts$Nmask_dark_blue

    ## Also add the percent of the spot covered, which can be useful
    ## for detecting neuropil spots
    spe$VistoSeg_proportion[spe$sample_id == sample_id] <-
        cell_counts$Pmask_dark_blue
}

#   Ensure counts were read in for all spots in the object
if (any(is.na(spe$VistoSeg_count))) {
    stop("Did not find cell counts for all non-IF spots.")
}
stopifnot(all(!is.na(spe$VistoSeg_proportion)))

#   Collect manual annotation of layers and wrinkles for the 3 annotated
#   samples, forming a list of tibbles
anno_list <- list()
for (sample_id in anno_samples) {
    #   Adjust paths for this sample
    this_anno_wrinkle_path <- sub(
        "\\{sample_id\\}", sample_id, anno_wrinkle_path
    )
    this_anno_layers_path <- sub("\\{sample_id\\}", sample_id, anno_layers_path)
    
    #   Read in wrinkle annotation and use unique and informative colnames
    anno_wrinkle <- this_anno_wrinkle_path |>
        read.csv() |>
        as_tibble() |>
        rename(wrinkle_type = ManualAnnotation, barcode = spot_name)
    
    #   Read in layer annotation and use unique and informative colnames
    anno_layers <- this_anno_layers_path |>
        read.csv() |>
        as_tibble() |>
        rename(manual_layer_label = ManualAnnotation, barcode = spot_name)
    
    anno_list[[sample_id]] <- anno_layers |>
        full_join(anno_wrinkle, by = c("sample_id", "barcode"))
}

#   Combine into a single tibble with all 3 samples
anno <- do.call(rbind, anno_list) |>
    as_tibble()

stopifnot(all(anno$barcode %in% colnames(spe)))

added_coldata <- tibble(
    "barcode" = colnames(spe), "sample_id" = spe$sample_id
) |>
    left_join(anno)

#   Add the wrinkle and layer annotation to the SPE's colData
colData(spe) <- cbind(
    colData(spe),
    added_coldata |> select(c("manual_layer_label", "wrinkle_type"))
)

#   Colnames started out sorted, so we'll mantain that organization
colData(spe) <- colData(spe)[, sort(colnames(colData(spe)))]

#   Interactively double-check that the annotations merged as expected

# table(is.na(spe$manual_layer_label))
# FALSE   TRUE 
# 11991 101936

# table(is.na(spe$wrinkle_type))
# FALSE   TRUE 
# 2094 111833

# vis_clus(spe, sampleid = anno_samples[1], clustervar = "wrinkle_type") +
#     #   Not sure why these are needed
#     scale_color_discrete() +
#     scale_fill_discrete()
# vis_clus(spe, sampleid = anno_samples[1], clustervar = "manual_layer_label")

## Check the size in GB
Sys.time()
lobstr::obj_size(spe) ## takes about 2 min
Sys.time()
# 6.96 GB

## Save for later use, takes about 11 min
Sys.time()
saveRDS(
    spe,
    file = here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
    )
)
Sys.time()


## Start removing pieces for the shiny apps
imgData(spe) <- imgData(spe)[imgData(spe)$image_id == "lowres", ]
## Check the size in GB
lobstr::obj_size(spe)
# 4.59 GB

## Drop the counts which take quite a bit of space
counts(spe) <- NULL
## Check the size in GB
lobstr::obj_size(spe)
# 2.45 GB


## Save for use in the spatialLIBD shiny apps, takes about 3 min
Sys.time()
saveRDS(
    spe,
    file = here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_subset_for_spatialLIBD.rds"
    )
)
Sys.time()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.2 Patched (2022-11-23 r83388)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2022-11-30
#  pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/pandoc
#
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────
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
#  BiocParallel             1.32.3    2022-11-23 [1] Bioconductor
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
#  digest                   0.6.30    2022-10-18 [2] CRAN (R 4.2.1)
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
#  fastmap                  1.1.0     2021-01-25 [2] CRAN (R 4.2.1)
#  fields                   14.1      2022-08-12 [2] CRAN (R 4.2.1)
#  filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
#  foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
#  fs                       1.5.2     2021-12-08 [2] CRAN (R 4.2.1)
#  gargle                   1.2.1     2022-09-08 [2] CRAN (R 4.2.1)
#  generics                 0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
#  GenomeInfoDb           * 1.34.3    2022-11-10 [1] Bioconductor
#  GenomeInfoDbData         1.2.9     2022-09-29 [2] Bioconductor
#  GenomicAlignments        1.34.0    2022-11-01 [2] Bioconductor
#  GenomicRanges          * 1.50.1    2022-11-06 [2] Bioconductor
#  ggbeeswarm               0.6.0     2017-08-07 [2] CRAN (R 4.2.1)
#  ggplot2                  3.4.0     2022-11-04 [2] CRAN (R 4.2.2)
#  ggrepel                  0.9.2     2022-11-06 [2] CRAN (R 4.2.2)
#  glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
#  golem                    0.3.5     2022-10-18 [2] CRAN (R 4.2.1)
#  googledrive              2.0.0     2021-07-08 [2] CRAN (R 4.2.1)
#  gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
#  gtable                   0.3.1     2022-09-01 [2] CRAN (R 4.2.1)
#  HDF5Array                1.26.0    2022-11-01 [2] Bioconductor
#  here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
#  htmltools                0.5.3     2022-07-18 [2] CRAN (R 4.2.1)
#  htmlwidgets              1.5.4     2021-09-08 [2] CRAN (R 4.2.1)
#  httpuv                   1.6.6     2022-09-08 [2] CRAN (R 4.2.1)
#  httr                     1.4.4     2022-08-17 [2] CRAN (R 4.2.1)
#  interactiveDisplayBase   1.36.0    2022-11-01 [2] Bioconductor
#  IRanges                * 2.32.0    2022-11-01 [2] Bioconductor
#  irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.2.1)
#  iterators                1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
#  jaffelab               * 0.99.32   2022-11-02 [1] Github (LieberInstitute/jaffelab@7b7afe3)
#  jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.2.1)
#  jsonlite                 1.8.3     2022-10-21 [2] CRAN (R 4.2.2)
#  KEGGREST                 1.38.0    2022-11-01 [2] Bioconductor
#  knitr                    1.41      2022-11-18 [2] CRAN (R 4.2.2)
#  later                    1.3.0     2021-08-18 [2] CRAN (R 4.2.1)
#  lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.2.2)
#  lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
#  lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
#  limma                    3.54.0    2022-11-01 [1] Bioconductor
#  lobstr                   1.1.2     2022-06-22 [2] CRAN (R 4.2.1)
#  locfit                   1.5-9.6   2022-07-11 [2] CRAN (R 4.2.1)
#  magick                   2.7.3     2021-08-18 [2] CRAN (R 4.2.1)
#  magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
#  maps                     3.4.1     2022-10-30 [2] CRAN (R 4.2.2)
#  MASS                     7.3-58.1  2022-08-03 [3] CRAN (R 4.2.2)
#  Matrix                   1.5-3     2022-11-11 [2] CRAN (R 4.2.2)
#  MatrixGenerics         * 1.10.0    2022-11-01 [1] Bioconductor
#  matrixStats            * 0.63.0    2022-11-18 [2] CRAN (R 4.2.2)
#  memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
#  mime                     0.12      2021-09-28 [2] CRAN (R 4.2.1)
#  munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
#  nlme                     3.1-160   2022-10-10 [2] CRAN (R 4.2.1)
#  paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.2.1)
#  pillar                   1.8.1     2022-08-19 [2] CRAN (R 4.2.1)
#  pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
#  pkgload                  1.3.2     2022-11-16 [2] CRAN (R 4.2.2)
#  plotly                   4.10.1    2022-11-07 [2] CRAN (R 4.2.2)
#  png                      0.1-8     2022-11-29 [2] CRAN (R 4.2.2)
#  prettyunits              1.1.1     2020-01-24 [2] CRAN (R 4.2.1)
#  promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.2.1)
#  purrr                    0.3.5     2022-10-06 [2] CRAN (R 4.2.1)
#  R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
#  R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
#  R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
#  R6                       2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
#  rafalib                * 1.0.0     2015-08-09 [1] CRAN (R 4.2.2)
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
#  roxygen2                 7.2.2     2022-11-11 [2] CRAN (R 4.2.2)
#  rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
#  Rsamtools                2.14.0    2022-11-01 [2] Bioconductor
#  RSQLite                  2.2.19    2022-11-24 [2] CRAN (R 4.2.2)
#  rstudioapi               0.14      2022-08-22 [2] CRAN (R 4.2.1)
#  rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
#  rtracklayer              1.58.0    2022-11-01 [2] Bioconductor
#  S4Vectors              * 0.36.0    2022-11-01 [1] Bioconductor
#  sass                     0.4.4     2022-11-24 [2] CRAN (R 4.2.2)
#  ScaledMatrix             1.6.0     2022-11-01 [1] Bioconductor
#  scales                   1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
#  scater                   1.26.1    2022-11-13 [2] Bioconductor
#  scuttle                  1.8.1     2022-11-20 [2] Bioconductor
#  segmented                1.6-1     2022-11-08 [1] CRAN (R 4.2.2)
#  servr                    0.25      2022-11-04 [1] CRAN (R 4.2.2)
#  sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
#  shiny                    1.7.3     2022-10-25 [2] CRAN (R 4.2.2)
#  shinyWidgets             0.7.5     2022-11-17 [2] CRAN (R 4.2.2)
#  SingleCellExperiment   * 1.20.0    2022-11-01 [1] Bioconductor
#  spam                     2.9-1     2022-08-07 [2] CRAN (R 4.2.1)
#  sparseMatrixStats        1.10.0    2022-11-01 [2] Bioconductor
#  SpatialExperiment      * 1.8.0     2022-11-01 [2] Bioconductor
#  spatialLIBD            * 1.11.1    2022-11-30 [1] Github (LieberInstitute/spatialLIBD@e16f4e0)
#  statmod                  1.4.37    2022-08-12 [2] CRAN (R 4.2.1)
#  stringi                  1.7.8     2022-07-11 [2] CRAN (R 4.2.1)
#  stringr                  1.4.1     2022-08-20 [2] CRAN (R 4.2.1)
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
#  xfun                     0.35      2022-11-16 [2] CRAN (R 4.2.2)
#  XML                      3.99-0.12 2022-10-28 [2] CRAN (R 4.2.2)
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
# ─────────────────────────────────────────────────────────────────────────────────────────────────
