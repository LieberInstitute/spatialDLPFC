library("SpatialExperiment")
library("scater")
library("here")
library("sessioninfo")
library("ggplot2")
library("Polychrome")

## Plot directory
dir_plots <- here::here(
    "plots",
    "09_position_differential_expression"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## load sce_pseudo data
sce_pseudo <-
    readRDS(
        here(
            "processed-data",
            "rdata",
            "spe",
            "07_layer_differential_expression",
            paste0("sce_pseudo_BayesSpace_k09.rds")
        )
    )

## Define variables to use
vars <- c(
    "age",
    "sample_id",
    "BayesSpace",
    "subject",
    "sex",
    "position"
)

## Obtain percent of variance explained at the gene level
## using scater::getVarianceExplained()
vars <- getVarianceExplained(sce_pseudo,
    variables = vars
)

## Now visualize the percent of variance explained across all genes
pdf(
    file = file.path(dir_plots, "sce_pseudo_gene_explanatory_vars_k09_large.pdf"),
    width = 18,
    height = 5
)
plotExplanatoryVariables(vars) + theme_classic(base_size = 30)
dev.off()

## Load Sp09 DE results
load(here("code", "deploy_app_k09", "sig_genes_subset_k09.Rdata"),
    verbose = TRUE
)
sig_domain <- sig_genes
rm(sig_genes)



## Load position DE results
load(
    here(
        "code",
        "deploy_app_k09_position",
        "sig_genes_subset_k09_position.Rdata"
    ),
    verbose = TRUE
)
sig_position <- sig_genes
rm(sig_genes)


## Combine results
enriched <- as.data.frame(rbind(
    cbind(subset(sig_domain, model_type == "enrichment"), "Analysis" = "Sp09"),
    cbind(subset(sig_position, model_type == "enrichment"), "Analysis" = "position")
))
enriched$Analysis <-
    factor(enriched$Analysis, levels = c("Sp09", "position"))

## Specify colors
colors <- palette36.colors(9 + 3)
names(colors) <- unique(enriched$test)

pdf(
    file.path(dir_plots, "density_histogram_enriched_tstats.pdf"),
    width = 18,
    height = 5
)
## Plot densities of t-statistics
ggplot(enriched, aes(x = stat, fill = test)) +
    geom_density() +
    facet_grid(~ Analysis + test, margins = "test") +
    xlab("Enrichment t-statistic") +
    scale_color_manual(values = colors, name = "Test") +
    scale_fill_manual(values = colors, name = "Test") +
    theme_bw(base_size = 20) +
    coord_flip() +
    scale_y_continuous(breaks = c(0, 0.3))

## Plot histogram of t-statistics with FDR < 5%
ggplot(subset(enriched, fdr < 0.05), aes(x = stat, fill = test)) +
    geom_histogram() +
    facet_grid(~ Analysis + test, margins = "test") +
    xlab("Enrichment t-statistic (FDR <5%)") +
    scale_color_manual(values = colors, name = "Test") +
    scale_fill_manual(values = colors, name = "Test") +
    theme_bw(base_size = 20) +
    coord_flip() +
    scale_y_log10() +
    theme(axis.text.x = element_text(size = 8)) +
    ylab("Count (log10 scale)")
dev.off()


## Summarize into a table


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.2 (2022-10-31)
#  os       macOS Ventura 13.0.1
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/Mexico_City
#  date     2023-02-06
#  rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
#  pandoc   2.17.1.1 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  beachmat               2.14.0    2022-11-01 [1] Bioconductor
#  beeswarm               0.4.0     2021-06-01 [1] CRAN (R 4.2.0)
#  Biobase              * 2.58.0    2022-11-01 [1] Bioconductor
#  BiocGenerics         * 0.44.0    2022-11-01 [1] Bioconductor
#  BiocNeighbors          1.16.0    2022-11-01 [1] Bioconductor
#  BiocParallel           1.32.5    2022-12-25 [1] Bioconductor
#  BiocSingular           1.14.0    2022-11-01 [1] Bioconductor
#  biocthis               1.9.1     2022-11-01 [1] Github (lcolladotor/biocthis@af38c7c)
#  bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.2.0)
#  brio                   1.1.3     2021-11-30 [1] CRAN (R 4.2.0)
#  cachem                 1.0.6     2021-08-19 [1] CRAN (R 4.2.0)
#  callr                  3.7.3     2022-11-02 [1] CRAN (R 4.2.2)
#  cli                    3.6.0     2023-01-09 [1] CRAN (R 4.2.0)
#  codetools              0.2-19    2023-02-01 [1] CRAN (R 4.2.0)
#  colorout               1.2-2     2022-03-01 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.2.0)
#  cowplot                1.1.1     2020-12-30 [1] CRAN (R 4.2.0)
#  crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.2.0)
#  data.table             1.14.6    2022-11-16 [1] CRAN (R 4.2.0)
#  DelayedArray           0.24.0    2022-11-01 [1] Bioconductor
#  DelayedMatrixStats     1.20.0    2022-11-01 [1] Bioconductor
#  devtools             * 2.4.5     2022-10-11 [1] CRAN (R 4.2.0)
#  digest                 0.6.31    2022-12-11 [1] CRAN (R 4.2.0)
#  dplyr                  1.1.0     2023-01-29 [1] CRAN (R 4.2.0)
#  dqrng                  0.3.0     2021-05-01 [1] CRAN (R 4.2.0)
#  DropletUtils           1.18.1    2022-11-23 [1] Bioconductor
#  edgeR                  3.40.2    2023-01-22 [1] Bioconductor
#  ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.2.0)
#  fansi                  1.0.4     2023-01-22 [1] CRAN (R 4.2.0)
#  farver                 2.1.1     2022-07-06 [1] CRAN (R 4.2.1)
#  fastmap                1.1.0     2021-01-25 [1] CRAN (R 4.2.0)
#  fs                     1.6.0     2023-01-23 [1] CRAN (R 4.2.0)
#  generics               0.1.3     2022-07-05 [1] CRAN (R 4.2.0)
#  GenomeInfoDb         * 1.34.9    2023-02-02 [1] Bioconductor
#  GenomeInfoDbData       1.2.9     2022-11-02 [1] Bioconductor
#  GenomicRanges        * 1.50.2    2022-12-18 [1] Bioconductor
#  ggbeeswarm             0.7.1     2022-12-16 [1] CRAN (R 4.2.2)
#  ggplot2              * 3.4.0     2022-11-04 [1] CRAN (R 4.2.0)
#  ggrepel                0.9.3     2023-02-03 [1] CRAN (R 4.2.0)
#  glue                   1.6.2     2022-02-24 [1] CRAN (R 4.2.0)
#  gridExtra              2.3       2017-09-09 [1] CRAN (R 4.2.0)
#  gtable                 0.3.1     2022-09-01 [1] CRAN (R 4.2.1)
#  HDF5Array              1.26.0    2022-11-01 [1] Bioconductor
#  here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.2.0)
#  hms                    1.1.2     2022-08-19 [1] CRAN (R 4.2.0)
#  htmltools              0.5.4     2022-12-07 [1] CRAN (R 4.2.0)
#  htmlwidgets            1.6.1     2023-01-07 [1] CRAN (R 4.2.0)
#  httpuv                 1.6.8     2023-01-12 [1] CRAN (R 4.2.0)
#  IRanges              * 2.32.0    2022-11-01 [1] Bioconductor
#  irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.2.1)
#  labeling               0.4.2     2020-10-20 [1] CRAN (R 4.2.0)
#  later                  1.3.0     2021-08-18 [1] CRAN (R 4.2.0)
#  lattice                0.20-45   2021-09-22 [1] CRAN (R 4.2.2)
#  lifecycle              1.0.3     2022-10-07 [1] CRAN (R 4.2.1)
#  limma                  3.54.1    2023-01-26 [1] Bioconductor
#  locfit                 1.5-9.7   2023-01-02 [1] CRAN (R 4.2.0)
#  lubridate              1.9.1     2023-01-24 [1] CRAN (R 4.2.0)
#  magick                 2.7.3     2021-08-18 [1] CRAN (R 4.2.0)
#  magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.2.0)
#  Matrix                 1.5-3     2022-11-11 [1] CRAN (R 4.2.0)
#  MatrixGenerics       * 1.10.0    2022-11-01 [1] Bioconductor
#  matrixStats          * 0.63.0    2022-11-18 [1] CRAN (R 4.2.0)
#  memoise                2.0.1     2021-11-26 [1] CRAN (R 4.2.0)
#  mime                   0.12      2021-09-28 [1] CRAN (R 4.2.0)
#  miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.2.0)
#  munsell                0.5.0     2018-06-12 [1] CRAN (R 4.2.0)
#  pillar                 1.8.1     2022-08-19 [1] CRAN (R 4.2.0)
#  pkgbuild               1.4.0     2022-11-27 [1] CRAN (R 4.2.2)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.2.0)
#  pkgload                1.3.2     2022-11-16 [1] CRAN (R 4.2.2)
#  Polychrome           * 1.5.1     2022-05-03 [1] CRAN (R 4.2.0)
#  prettyunits            1.1.1     2020-01-24 [1] CRAN (R 4.2.0)
#  processx               3.8.0     2022-10-26 [1] CRAN (R 4.2.0)
#  profvis                0.3.7     2020-11-02 [1] CRAN (R 4.2.0)
#  promises               1.2.0.1   2021-02-11 [1] CRAN (R 4.2.0)
#  prompt                 1.0.1     2022-03-01 [1] Github (gaborcsardi/prompt@7ef0f2e)
#  ps                     1.7.2     2022-10-26 [1] CRAN (R 4.2.0)
#  purrr                  1.0.1     2023-01-10 [1] CRAN (R 4.2.0)
#  R.cache                0.16.0    2022-07-21 [1] CRAN (R 4.2.0)
#  R.methodsS3            1.8.2     2022-06-13 [1] CRAN (R 4.2.0)
#  R.oo                   1.25.0    2022-06-12 [1] CRAN (R 4.2.0)
#  R.utils                2.12.2    2022-11-11 [1] CRAN (R 4.2.0)
#  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.2.0)
#  Rcpp                   1.0.10    2023-01-22 [1] CRAN (R 4.2.0)
#  RCurl                  1.98-1.10 2023-01-27 [1] CRAN (R 4.2.0)
#  remotes                2.4.2     2021-11-30 [1] CRAN (R 4.2.0)
#  rhdf5                  2.42.0    2022-11-01 [1] Bioconductor
#  rhdf5filters           1.10.0    2022-11-01 [1] Bioconductor
#  Rhdf5lib               1.20.0    2022-11-01 [1] Bioconductor
#  rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.2.0)
#  rlang                  1.0.6     2022-09-24 [1] CRAN (R 4.2.0)
#  rprojroot              2.0.3     2022-04-02 [1] CRAN (R 4.2.0)
#  rsthemes               0.3.1     2022-03-01 [1] Github (gadenbuie/rsthemes@bbe73ca)
#  rstudioapi             0.14      2022-08-22 [1] CRAN (R 4.2.0)
#  rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.2.0)
#  S4Vectors            * 0.36.1    2022-12-07 [1] Bioconductor
#  ScaledMatrix           1.6.0     2022-11-01 [1] Bioconductor
#  scales                 1.2.1     2022-08-20 [1] CRAN (R 4.2.0)
#  scater               * 1.26.1    2022-11-13 [1] Bioconductor
#  scatterplot3d          0.3-42    2022-09-08 [1] CRAN (R 4.2.0)
#  scuttle              * 1.8.4     2023-01-22 [1] Bioconductor
#  sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.2.0)
#  shiny                  1.7.4     2022-12-15 [1] CRAN (R 4.2.2)
#  SingleCellExperiment * 1.20.0    2022-11-01 [1] Bioconductor
#  sparseMatrixStats      1.10.0    2022-11-01 [1] Bioconductor
#  SpatialExperiment    * 1.8.0     2022-11-01 [1] Bioconductor
#  stringi                1.7.12    2023-01-11 [1] CRAN (R 4.2.0)
#  stringr                1.5.0     2022-12-02 [1] CRAN (R 4.2.0)
#  styler                 1.9.0     2023-01-15 [1] CRAN (R 4.2.0)
#  SummarizedExperiment * 1.28.0    2022-11-01 [1] Bioconductor
#  suncalc                0.5.1     2022-09-29 [1] CRAN (R 4.2.0)
#  testthat             * 3.1.6     2022-12-09 [1] CRAN (R 4.2.0)
#  tibble                 3.1.8     2022-07-22 [1] CRAN (R 4.2.1)
#  tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.2.0)
#  timechange             0.2.0     2023-01-11 [1] CRAN (R 4.2.0)
#  urlchecker             1.0.1     2021-11-30 [1] CRAN (R 4.2.0)
#  usethis              * 2.1.6     2022-05-25 [1] CRAN (R 4.2.0)
#  utf8                   1.2.3     2023-01-31 [1] CRAN (R 4.2.0)
#  vctrs                  0.5.2     2023-01-23 [1] CRAN (R 4.2.0)
#  vipor                  0.4.5     2017-03-22 [1] CRAN (R 4.2.0)
#  viridis                0.6.2     2021-10-13 [1] CRAN (R 4.2.0)
#  viridisLite            0.4.1     2022-08-22 [1] CRAN (R 4.2.0)
#  withr                  2.5.0     2022-03-03 [1] CRAN (R 4.2.0)
#  xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.2.0)
#  XVector                0.38.0    2022-11-01 [1] Bioconductor
#  zlibbioc               1.44.0    2022-11-01 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
