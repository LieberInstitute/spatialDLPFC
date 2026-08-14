suppressPackageStartupMessages(library("spatialLIBD"))
library("here")
library("Banksy")
library("tidyr")
library("ggplot2")
library("ggpubr")
library("RColorBrewer")
library("ggrepel")
library("sessioninfo")

## Download HumanPilot project data
spe <- spatialLIBD::fetch_data(type = "spe")

bayes_res <- dir(pattern = "BayesSpace_version")
bayes_res
for(i in seq_along(bayes_res)) {
  spe <- cluster_import(spe, bayes_res[i], prefix = bayes_res[i])
}
colnames(colData(spe)) <- gsub("spatial.cluster", "", colnames(colData(spe)))


spe_sub <- spe[, !is.na(spe$layer_guess_reordered_short)]
colData(spe_sub) <- colData(spe_sub)[, c("layer_guess_reordered_short", bayes_res)]
colnames(colData(spe_sub)) <- paste0("clust_", colnames(colData(spe_sub)))

ari <- compareClusters(spe_sub, func = "ARI")
rownames(ari) <- gsub("clust_", "", rownames(ari))
colnames(ari) <- gsub("clust_", "", colnames(ari))
ari[-6, -6]
#                             layer_guess_reordered_short BayesSpace_version_1.14.0 BayesSpace_version_1.16.0 BayesSpace_version_1.17.0 BayesSpace_version_1.20.2
# layer_guess_reordered_short                       1.000                     0.396                     0.407                     0.407                     0.407
# BayesSpace_version_1.14.0                         0.396                     1.000                     0.856                     0.856                     0.856
# BayesSpace_version_1.16.0                         0.407                     0.856                     1.000                     1.000                     1.000
# BayesSpace_version_1.17.0                         0.407                     0.856                     1.000                     1.000                     1.000
# BayesSpace_version_1.20.2                         0.407                     0.856                     1.000                     1.000                     1.000

nmi <- compareClusters(spe_sub, func = "NMI")
rownames(nmi) <- gsub("clust_", "", rownames(nmi))
colnames(nmi) <- gsub("clust_", "", colnames(nmi))
nmi[-6, -6]
#                             layer_guess_reordered_short BayesSpace_version_1.14.0 BayesSpace_version_1.16.0 BayesSpace_version_1.17.0 BayesSpace_version_1.20.2
# layer_guess_reordered_short                       1.000                     0.537                     0.544                     0.544                     0.544
# BayesSpace_version_1.14.0                         0.537                     1.000                     0.842                     0.842                     0.842
# BayesSpace_version_1.16.0                         0.544                     0.842                     1.000                     1.000                     1.000
# BayesSpace_version_1.17.0                         0.544                     0.842                     1.000                     1.000                     1.000
# BayesSpace_version_1.20.2                         0.544                     0.842                     1.000                     1.000                     1.000

## Compute them by sample
ari_sample <- lapply(unique(spe_sub$sample_id), function(x) {
  spe_sub2 <- spe_sub[, spe_sub$sample_id == x]
  x <- compareClusters(spe_sub2, func = "ARI")
  rownames(x) <- gsub("clust_", "", rownames(x))
  colnames(x) <- gsub("clust_", "", colnames(x))
  x[-6, -6] ## Drop the "sample_id" comparison
})

nmi_sample <- lapply(unique(spe_sub$sample_id), function(x) {
  spe_sub2 <- spe_sub[, spe_sub$sample_id == x]
  x <- compareClusters(spe_sub2, func = "NMI")
  rownames(x) <- gsub("clust_", "", rownames(x))
  colnames(x) <- gsub("clust_", "", colnames(x))
  x[-6, -6] ## Drop the "sample_id" comparison
})
names(ari_sample) <- names(nmi_sample) <- unique(spe_sub$sample_id)

## Tidy into long format the information by sample
ari_sample_long <- lapply(names(ari_sample), function(x) {
  x2 <- ari_sample[[x]]
  x2 <- as.data.frame(x2)
  x2$Var1 <- rownames(x2)
  x2 <- tidyr::pivot_longer(x2, cols = -Var1, names_to = "Var2", values_to = "ARI")
  x2$sample_id <- x
  x2
})
ari_sample_long <- do.call(rbind, ari_sample_long)

nmi_sample_long <- lapply(names(nmi_sample), function(x) {
  x2 <- nmi_sample[[x]]
  x2 <- as.data.frame(x2)
  x2$Var1 <- rownames(x2)
  x2 <- tidyr::pivot_longer(x2, cols = -Var1, names_to = "Var2", values_to = "NMI")
  x2$sample_id <- x
  x2
})
nmi_sample_long <- do.call(rbind, nmi_sample_long)

comp_values <- merge(ari_sample_long, nmi_sample_long, by = c("Var1", "Var2", "sample_id"))
colnames(comp_values) <- c("Reference", "Target", "sample_id", "ARI", "NMI")
comp_values$Reference <- gsub("BayesSpace_version_", "", comp_values$Reference)
comp_values$Target <- gsub("BayesSpace_version_", "", comp_values$Target)

## Plot ARI and NMI values by sample
## Code adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/05_ARI/02_plot_ARI.R
p_ari <- ggpaired(
    subset(comp_values, Reference == "layer_guess_reordered_short" & ! Target %in% c("layer_guess_reordered_short", "1.17.0")),
    x = "Target",
    y = "ARI",
    id = "sample_id",
    xlab = "BayesSpace version (from BioC)",
    ylab = "Adjusted Rand Index",
    palette = brewer.pal("Paired", n = 10)[c(1:4, 10)],
    fill = "Target",
    line.color = "gray",
    line.size = 0.4,
    point.size = 3
) + font("xy.text", size = 20) +
    font("xlab", size = 30) +
    font("ylab", size = 30) + ggrepel::geom_text_repel(
        aes(label = sample_id),
        max.overlaps = Inf,
        fontface = "italic",
        size = 3,
        seed = 20221216,
        force = 30,
        max.time = 4,
        segment.color = "lightpink",
        color = "deeppink"
    ) + theme(legend.position = "none")

p_nmi <- ggpaired(
    subset(comp_values, Reference == "layer_guess_reordered_short" & ! Target %in% c("layer_guess_reordered_short", "1.17.0")),
    x = "Target",
    y = "NMI",
    id = "sample_id",
    xlab = "BayesSpace version (from BioC)",
    ylab = "Normalized Mutual Information",
    palette = brewer.pal("Paired", n = 10)[c(1:4, 10)],
    fill = "Target",
    line.color = "gray",
    line.size = 0.4,
    point.size = 3
) + font("xy.text", size = 20) +
    font("xlab", size = 30) +
    font("ylab", size = 30) + ggrepel::geom_text_repel(
        aes(label = sample_id),
        max.overlaps = Inf,
        fontface = "italic",
        size = 3,
        seed = 20221216,
        force = 30,
        max.time = 4,
        segment.color = "lightpink",
        color = "deeppink"
    ) + theme(legend.position = "none")

p_ari
p_nmi

## Repeat using version 1.14.0 as the reference
p_ari_1.14.0 <- ggpaired(
    subset(comp_values, Reference == "1.14.0" & ! Target %in% c("layer_guess_reordered_short", "1.17.0", "1.14.0")),
    x = "Target",
    y = "ARI",
    id = "sample_id",
    xlab = "BayesSpace version (from BioC)",
    ylab = "Adjusted Rand Index",
    palette = brewer.pal("Paired", n = 10)[c(1:4, 10)],
    fill = "Target",
    line.color = "gray",
    line.size = 0.4,
    point.size = 3
) + font("xy.text", size = 20) +
    font("xlab", size = 30) +
    font("ylab", size = 30) + ggrepel::geom_text_repel(
        aes(label = sample_id),
        max.overlaps = Inf,
        fontface = "italic",
        size = 3,
        seed = 20221216,
        force = 30,
        max.time = 4,
        segment.color = "lightpink",
        color = "deeppink"
    ) + theme(legend.position = "none")

p_nmi_1.14.0 <- ggpaired(
    subset(comp_values, Reference == "1.14.0" & ! Target %in% c("layer_guess_reordered_short", "1.17.0", "1.14.0")),
    x = "Target",
    y = "NMI",
    id = "sample_id",
    xlab = "BayesSpace version (from BioC)",
    ylab = "Normalized Mutual Information",
    palette = brewer.pal("Paired", n = 10)[c(1:4, 10)],
    fill = "Target",
    line.color = "gray",
    line.size = 0.4,
    point.size = 3
) + font("xy.text", size = 20) +
    font("xlab", size = 30) +
    font("ylab", size = 30) + ggrepel::geom_text_repel(
        aes(label = sample_id),
        max.overlaps = Inf,
        fontface = "italic",
        size = 3,
        seed = 20221216,
        force = 30,
        max.time = 4,
        segment.color = "lightpink",
        color = "deeppink"
    ) + theme(legend.position = "none")

p_ari_1.14.0
p_nmi_1.14.0

## Compare changes across versions in BayesSpace
tab <- table(
  "BayesSpace 1.14.0" = spe_sub$clust_BayesSpace_version_1.14.0,
  "BayesSpace1.16.0" = spe_sub$clust_BayesSpace_version_1.16.0
)
addmargins(tab)
#                  BayesSpace1.16.0
# BayesSpace 1.14.0     1     2     3     4     5     6     7   Sum
#               1   11045     0     0     9   371   103    25 11553
#               2      30  2222    50     4     6   165    27  2504
#               3       0   110  1895     1     0     0     7  2013
#               4       8     0     0  4779    37     5    64  4893
#               5     577     0     0   313 10399     1     9 11299
#               6     733    27     0     7     1 12257    67 13092
#               7       4    23     0    13     6    63  1866  1975
#               Sum 12397  2382  1945  5126 10820 12594  2065 47329

## Percent agreement between BayesSpace 1.14.0 and 1.16.0
round(sum(diag(tab)) / sum(tab) * 100, 2)
# [1] 93.94

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120) ## Makes it easier to read later
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.6.1 (2026-06-24)
#  os       macOS Tahoe 26.6
#  system   aarch64, darwin23
#  ui       Positron
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2026-08-13
#  pandoc   3.10.1 @ /opt/homebrew/bin/pandoc
#  quarto   1.10.18 @ /Applications/quarto/bin/quarto

# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version    date (UTC) lib source
#  abind                  1.4-8      2024-09-12 [1] CRAN (R 4.6.0)
#  AnnotationDbi          1.74.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  AnnotationHub          4.2.2      2026-06-30 [1] Bioconductor 3.23 (R 4.6.1)
#  aricode                1.1.0      2026-05-13 [1] CRAN (R 4.6.0)
#  attempt                0.3.1      2020-05-03 [1] CRAN (R 4.6.0)
#  backports              1.5.1      2026-04-03 [1] CRAN (R 4.6.0)
#  Banksy               * 1.8.1      2026-05-02 [1] https://bioc-release.r-universe.dev (R 4.6.0)
#  beachmat               2.28.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  beeswarm               0.4.0      2021-06-01 [1] CRAN (R 4.6.0)
#  benchmarkme            1.0.8      2022-06-12 [1] CRAN (R 4.6.0)
#  benchmarkmeData        2.0.0      2026-01-19 [1] CRAN (R 4.6.0)
#  Biobase              * 2.72.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  BiocFileCache          3.2.0      2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  BiocGenerics         * 0.58.1     2026-05-14 [1] https://bioc-release.r-universe.dev (R 4.6.0)
#  BiocIO                 1.22.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  BiocManager            1.30.27    2025-11-14 [1] CRAN (R 4.6.0)
#  BiocNeighbors          2.6.0      2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  BiocParallel           1.46.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  BiocSingular           1.28.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  BiocVersion            3.23.1     2025-10-30 [1] https://bioc.r-universe.dev (R 4.6.0)
#  Biostrings             2.80.1     2026-05-22 [1] https://bioc-release.r-universe.dev (R 4.6.0)
#  bit                    4.6.0      2025-03-06 [1] CRAN (R 4.6.0)
#  bit64                  4.8.2      2026-05-19 [1] CRAN (R 4.6.0)
#  bitops                 1.1-0      2026-07-30 [1] CRAN (R 4.6.1)
#  blob                   1.3.0      2026-01-14 [1] CRAN (R 4.6.0)
#  brio                   1.1.5      2024-04-24 [1] CRAN (R 4.6.0)
#  broom                  1.0.13     2026-05-14 [1] CRAN (R 4.6.0)
#  bslib                  0.12.0     2026-08-04 [1] CRAN (R 4.6.1)
#  cachem                 1.1.0      2024-05-16 [1] CRAN (R 4.6.0)
#  callr                  3.8.0      2026-06-05 [1] CRAN (R 4.6.0)
#  car                    3.1-5      2026-02-03 [1] CRAN (R 4.6.0)
#  carData                3.0-6      2026-01-30 [1] CRAN (R 4.6.0)
#  cigarillo              1.2.1      2026-07-13 [1] Bioconductor 3.23 (R 4.6.1)
#  circlize               0.4.18     2026-04-04 [1] CRAN (R 4.6.0)
#  cli                    3.6.6      2026-04-09 [1] CRAN (R 4.6.0)
#  clue                   0.3-68     2026-03-26 [1] CRAN (R 4.6.0)
#  cluster                2.1.8.3    2026-07-30 [1] CRAN (R 4.6.1)
#  codetools              0.2-20     2024-03-31 [1] CRAN (R 4.6.1)
#  colorout             * 1.3-3      2026-07-30 [1] Github (jalvesaq/colorout@64863bb)
#  colorspace             2.1-3      2026-07-12 [1] CRAN (R 4.6.1)
#  ComplexHeatmap         2.28.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  config                 0.3.2      2023-08-30 [1] CRAN (R 4.6.0)
#  coro                   1.1.0      2024-11-05 [1] CRAN (R 4.6.0)
#  cowplot                1.2.0      2025-07-07 [1] CRAN (R 4.6.0)
#  crayon                 1.5.3      2024-06-20 [1] CRAN (R 4.6.0)
#  curl                   7.1.0      2026-04-22 [1] CRAN (R 4.6.0)
#  data.table             1.18.4     2026-05-06 [1] CRAN (R 4.6.0)
#  DBI                    1.3.0      2026-02-25 [1] CRAN (R 4.6.0)
#  dbplyr                 2.6.0      2026-06-17 [1] CRAN (R 4.6.0)
#  dbscan                 1.2.5      2026-06-09 [1] CRAN (R 4.6.0)
#  DelayedArray           0.38.2     2026-05-26 [1] https://bioc-release.r-universe.dev (R 4.6.0)
#  devtools             * 2.5.2      2026-04-30 [1] CRAN (R 4.6.0)
#  dichromat              2.0-1      2026-07-22 [1] CRAN (R 4.6.1)
#  digest                 0.6.39     2025-11-19 [1] CRAN (R 4.6.0)
#  doParallel             1.0.17     2022-02-07 [1] CRAN (R 4.6.0)
#  dplyr                  1.2.1      2026-04-03 [1] CRAN (R 4.6.0)
#  DT                     0.34.0     2025-09-02 [1] CRAN (R 4.6.0)
#  edgeR                  4.10.1     2026-05-23 [1] https://bioc-release.r-universe.dev (R 4.6.0)
#  ellipsis               0.3.3      2026-04-04 [1] CRAN (R 4.6.0)
#  ellmer                 0.4.2      2026-07-13 [1] CRAN (R 4.6.1)
#  ExperimentHub          3.2.0      2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  farver                 2.1.2      2024-05-13 [1] CRAN (R 4.6.0)
#  fastmap                1.2.0      2024-05-15 [1] CRAN (R 4.6.0)
#  filelock               1.0.3      2023-12-11 [1] CRAN (R 4.6.0)
#  foreach                1.5.2      2022-02-02 [1] CRAN (R 4.6.0)
#  Formula                1.2-6      2026-08-03 [1] CRAN (R 4.6.1)
#  fs                     2.1.0      2026-04-18 [1] CRAN (R 4.6.0)
#  generics             * 0.1.4      2025-05-09 [1] CRAN (R 4.6.0)
#  GenomicAlignments      1.48.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  GenomicRanges        * 1.64.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  GetoptLong             1.1.1      2026-04-08 [1] CRAN (R 4.6.0)
#  ggbeeswarm             0.7.3      2025-11-29 [1] CRAN (R 4.6.0)
#  ggplot2              * 4.0.3      2026-04-22 [1] CRAN (R 4.6.0)
#  ggpubr               * 1.0.0      2026-07-06 [1] CRAN (R 4.6.1)
#  ggrepel              * 0.9.8      2026-03-17 [1] CRAN (R 4.6.0)
#  ggsignif               0.6.4      2022-10-13 [1] CRAN (R 4.6.0)
#  gitcreds               0.1.2      2022-09-08 [1] CRAN (R 4.6.0)
#  GlobalOptions          0.1.4      2026-04-08 [1] CRAN (R 4.6.0)
#  glue                   1.8.1      2026-04-17 [1] CRAN (R 4.6.0)
#  golem                  1.0.1      2026-07-07 [1] CRAN (R 4.6.1)
#  gridExtra              2.3.1      2026-06-25 [1] CRAN (R 4.6.1)
#  gtable                 0.3.6      2024-10-25 [1] CRAN (R 4.6.0)
#  here                 * 1.0.2      2025-09-15 [1] CRAN (R 4.6.0)
#  htmltools              0.5.9      2025-12-04 [1] CRAN (R 4.6.0)
#  htmlwidgets            1.6.4      2023-12-06 [1] CRAN (R 4.6.0)
#  httpuv                 1.6.17     2026-03-18 [1] CRAN (R 4.6.0)
#  httr                   1.4.8      2026-02-13 [1] CRAN (R 4.6.0)
#  httr2                  1.3.0      2026-07-13 [1] CRAN (R 4.6.1)
#  igraph                 2.3.3      2026-06-26 [1] CRAN (R 4.6.1)
#  IRanges              * 2.46.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  irlba                  2.3.7      2026-01-30 [1] CRAN (R 4.6.0)
#  iterators              1.0.14     2022-02-05 [1] CRAN (R 4.6.0)
#  jquerylib              0.1.4      2021-04-26 [1] CRAN (R 4.6.0)
#  jsonlite               2.0.0      2025-03-27 [1] CRAN (R 4.6.0)
#  KEGGREST               1.52.2     2026-06-16 [1] Bioconductor 3.23 (R 4.6.1)
#  labeling               0.4.3      2023-08-29 [1] CRAN (R 4.6.0)
#  later                  1.4.8      2026-03-05 [1] CRAN (R 4.6.0)
#  lattice                0.22-9     2026-02-09 [1] CRAN (R 4.6.1)
#  leidenAlg              1.1.8      2026-05-31 [1] CRAN (R 4.6.0)
#  lifecycle              1.0.5      2026-01-08 [1] CRAN (R 4.6.0)
#  limma                  3.68.4     2026-05-31 [1] https://bioc-release.r-universe.dev (R 4.6.0)
#  locfit                 1.5-9.12   2025-03-05 [1] CRAN (R 4.6.0)
#  magick                 2.9.1      2026-02-28 [1] CRAN (R 4.6.0)
#  magrittr               2.0.5      2026-04-04 [1] CRAN (R 4.6.0)
#  Matrix                 1.7-6      2026-07-25 [1] CRAN (R 4.6.1)
#  MatrixGenerics       * 1.24.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  matrixStats          * 1.5.0      2025-01-07 [1] CRAN (R 4.6.0)
#  mclust                 6.1.3      2026-07-05 [1] CRAN (R 4.6.1)
#  memoise                2.0.1      2021-11-26 [1] CRAN (R 4.6.0)
#  mime                   0.13       2025-03-17 [1] CRAN (R 4.6.0)
#  otel                   0.2.0      2025-08-29 [1] CRAN (R 4.6.0)
#  pak                    0.11.1     2026-07-22 [1] CRAN (R 4.6.1)
#  paletteer              1.7.0      2026-01-08 [1] CRAN (R 4.6.0)
#  pillar                 1.11.1     2025-09-17 [1] CRAN (R 4.6.0)
#  pkgbuild               1.4.8      2025-05-26 [1] CRAN (R 4.6.0)
#  pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.6.0)
#  pkgload                1.5.3      2026-06-15 [1] CRAN (R 4.6.0)
#  plotly                 4.12.1     2026-07-22 [1] CRAN (R 4.6.1)
#  png                    0.1-9      2026-03-15 [1] CRAN (R 4.6.0)
#  processx               3.9.0      2026-04-22 [1] CRAN (R 4.6.0)
#  promises               1.5.0      2025-11-01 [1] CRAN (R 4.6.0)
#  prompt                 1.0.2.9000 2026-07-30 [1] Github (gaborcsardi/prompt@17bd0e1)
#  ps                     1.9.3      2026-04-20 [1] CRAN (R 4.6.0)
#  purrr                  1.2.2      2026-04-10 [1] CRAN (R 4.6.0)
#  R6                     2.6.1      2025-02-15 [1] CRAN (R 4.6.0)
#  ragg                   1.5.2      2026-03-23 [1] CRAN (R 4.6.0)
#  rappdirs               0.3.4      2026-01-17 [1] CRAN (R 4.6.0)
#  RColorBrewer         * 1.1-3      2022-04-03 [1] CRAN (R 4.6.0)
#  Rcpp                   1.1.2      2026-07-05 [1] CRAN (R 4.6.1)
#  RcppHungarian          0.3        2023-09-05 [1] CRAN (R 4.6.0)
#  RCurl                  1.98-1.19  2026-06-03 [1] CRAN (R 4.6.0)
#  rematch2               2.1.2      2020-05-01 [1] CRAN (R 4.6.0)
#  restfulr               0.0.17     2026-06-11 [1] CRAN (R 4.6.0)
#  rjson                  0.2.23     2024-09-16 [1] CRAN (R 4.6.0)
#  rlang                  1.3.0      2026-07-05 [1] CRAN (R 4.6.1)
#  rprojroot              2.1.1      2025-08-26 [1] CRAN (R 4.6.0)
#  Rsamtools              2.28.0     2026-04-29 [1] Bioconductor 3.23 (R 4.6.0)
#  RSQLite                3.53.3     2026-06-30 [1] CRAN (R 4.6.1)
#  rstatix                1.1.0      2026-07-23 [1] CRAN (R 4.6.1)
#  rsthemes               0.5.1      2025-10-29 [1] https://gadenbuie.r-universe.dev (R 4.6.1)
#  rstudioapi             0.19.0     2026-06-11 [1] CRAN (R 4.6.0)
#  rsvd                   1.0.5      2021-04-16 [1] CRAN (R 4.6.0)
#  rtracklayer            1.72.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  S4Arrays               1.12.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  S4Vectors            * 0.50.1     2026-05-03 [1] https://bioc-release.r-universe.dev (R 4.6.0)
#  S7                     0.2.2      2026-04-22 [1] CRAN (R 4.6.0)
#  sass                   0.4.10     2025-04-11 [1] CRAN (R 4.6.0)
#  ScaledMatrix           1.20.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  scales                 1.4.0      2025-04-24 [1] CRAN (R 4.6.0)
#  scater                 1.40.2     2026-06-30 [1] Bioconductor 3.23 (R 4.6.1)
#  sccore                 1.0.7      2026-04-06 [1] CRAN (R 4.6.0)
#  scuttle                1.22.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  Seqinfo              * 1.2.0      2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  sessioninfo          * 1.2.4      2026-06-04 [1] CRAN (R 4.6.0)
#  shape                  1.4.6.1    2024-02-23 [1] CRAN (R 4.6.0)
#  shiny                  1.14.0     2026-06-21 [1] CRAN (R 4.6.0)
#  shinyWidgets           0.9.1      2026-03-09 [1] CRAN (R 4.6.0)
#  SingleCellExperiment * 1.34.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  SparseArray            1.12.2     2026-05-01 [1] Bioconductor 3.23 (R 4.6.0)
#  SpatialExperiment    * 1.22.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  spatialLIBD          * 1.25.2     2026-07-31 [1] Github (LieberInstitute/spatialLIBD@94b68c6)
#  statmod                1.5.2      2026-05-17 [1] CRAN (R 4.6.0)
#  stringi                1.8.9      2026-08-04 [1] CRAN (R 4.6.1)
#  stringr                1.6.0      2025-11-04 [1] CRAN (R 4.6.0)
#  SummarizedExperiment * 1.42.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  systemfonts            1.3.2      2026-03-05 [1] CRAN (R 4.6.0)
#  testthat             * 3.3.2      2026-01-11 [1] CRAN (R 4.6.0)
#  textshaping            1.0.5      2026-03-06 [1] CRAN (R 4.6.0)
#  tibble                 3.3.1      2026-01-11 [1] CRAN (R 4.6.0)
#  tidyr                * 1.3.2      2025-12-19 [1] CRAN (R 4.6.0)
#  tidyselect             1.2.1      2024-03-11 [1] CRAN (R 4.6.0)
#  usethis              * 3.2.1      2025-09-06 [1] CRAN (R 4.6.0)
#  utf8                   1.2.6      2025-06-08 [1] CRAN (R 4.6.0)
#  uwot                   0.2.4      2025-11-10 [1] CRAN (R 4.6.0)
#  vctrs                  0.7.3      2026-04-11 [1] CRAN (R 4.6.0)
#  vipor                  0.4.7      2023-12-18 [1] CRAN (R 4.6.0)
#  viridis                0.6.5      2024-01-29 [1] CRAN (R 4.6.0)
#  viridisLite            0.4.3      2026-02-04 [1] CRAN (R 4.6.0)
#  withr                  3.0.3      2026-06-19 [1] CRAN (R 4.6.0)
#  XML                    3.99-0.23  2026-03-20 [1] CRAN (R 4.6.0)
#  xtable                 1.8-8      2026-02-22 [1] CRAN (R 4.6.0)
#  XVector                0.52.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
#  yaml                   2.3.12     2025-12-10 [1] CRAN (R 4.6.0)

#  [1] /Library/Frameworks/R.framework/Versions/4.6/Resources/library
#  * ── Packages attached to the search path.

# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────