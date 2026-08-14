suppressPackageStartupMessages(library("spatialLIBD"))
library("here")
library("Banksy")
library("tidyr")
library("ggplot2")
library("ggpubr")
library("RColorBrewer")
library("ggrepel")
library("sessioninfo")

## plots directory
dir_plots <- here::here("plots", "03_BayesSpace", "01_check_BayesSpace_versions")
dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)

## Download spatialDLPFC project data
spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")

## Locate and read in the files that Ben had created
files_ben <- dir(here::here("processed-data", "rdata", "spe", "clustering_results", "BayesSpace_versions"), "_output_cluster_assignments.csv", full.names = TRUE)
names(files_ben) <- gsub("_output_cluster_assignments.csv", "", basename(files_ben))
res_ben <- lapply(files_ben, function(x) {
  y <- read.csv(x, header = TRUE)
  y$key <- paste0(y$barcode, "_", y$sample_id)
  m <- match(spe$key, y$key)
  y[m, ]
  return(y)
})
spe$v1.5.1_ben <- res_ben[["1.5.1"]]$cluster
spe$v1.16.0_ben <- res_ben[["1.16.0"]]$cluster

bayes_res <- dir(here::here("processed-data", "rdata", "spe", "clustering_results", "BayesSpace_versions"), "BayesSpace_version_", full.names = TRUE)
names(bayes_res) <- gsub("BayesSpace_version_", "", basename(bayes_res))
bayes_res
for(i in seq_along(bayes_res)) {
  spe <- cluster_import(spe, bayes_res[i], prefix = names(bayes_res)[i])
}
colnames(colData(spe)) <- gsub("[X]*1\\.\\d+\\.\\d+BayesSpace_version_", "v", colnames(colData(spe)))

clust_interest <- c("BayesSpace_harmony_09", paste0("v", rev(names(files_ben)), "_ben"), paste0("v", names(bayes_res)))
spe_sub <- spe
colData(spe_sub) <- colData(spe)[, clust_interest]
for(i in clust_interest) {
  colData(spe_sub)[, i] <- as.factor(colData(spe_sub)[, i])
}
colnames(colData(spe_sub)) <- c(paste0("clust_", clust_interest), "sample_id")
head(colData(spe_sub))

ari <- compareClusters(spe_sub, func = "ARI")
rownames(ari) <- gsub("clust_", "", rownames(ari))
colnames(ari) <- gsub("clust_", "", colnames(ari))
ari
#                       BayesSpace_harmony_09 v1.5.1_ben v1.16.0_ben v1.14.0 v1.16.0 v1.17.0 v1.20.2
# BayesSpace_harmony_09                 1.000      1.000       0.305   1.000   0.899   0.899   0.899
# v1.5.1_ben                            1.000      1.000       0.305   1.000   0.899   0.899   0.899
# v1.16.0_ben                           0.305      0.305       1.000   0.305   0.300   0.300   0.300
# v1.14.0                               1.000      1.000       0.305   1.000   0.899   0.899   0.899
# v1.16.0                               0.899      0.899       0.300   0.899   1.000   1.000   1.000
# v1.17.0                               0.899      0.899       0.300   0.899   1.000   1.000   1.000
# v1.20.2                               0.899      0.899       0.300   0.899   1.000   1.000   1.000

nmi <- compareClusters(spe_sub, func = "NMI")
rownames(nmi) <- gsub("clust_", "", rownames(nmi))
colnames(nmi) <- gsub("clust_", "", colnames(nmi))
nmi
#                       BayesSpace_harmony_09 v1.5.1_ben v1.16.0_ben v1.14.0 v1.16.0 v1.17.0 v1.20.2
# BayesSpace_harmony_09                 1.000      1.000       0.423   1.000   0.885   0.885   0.885
# v1.5.1_ben                            1.000      1.000       0.423   1.000   0.885   0.885   0.885
# v1.16.0_ben                           0.423      0.423       1.000   0.423   0.414   0.414   0.414
# v1.14.0                               1.000      1.000       0.423   1.000   0.885   0.885   0.885
# v1.16.0                               0.885      0.885       0.414   0.885   1.000   1.000   1.000
# v1.17.0                               0.885      0.885       0.414   0.885   1.000   1.000   1.000
# v1.20.2                               0.885      0.885       0.414   0.885   1.000   1.000   1.000

## Compute them by sample
ari_sample <- lapply(unique(spe_sub$sample_id), function(x) {
  spe_sub2 <- spe_sub[, spe_sub$sample_id == x]
  x <- compareClusters(spe_sub2, func = "ARI")
  rownames(x) <- gsub("clust_", "", rownames(x))
  colnames(x) <- gsub("clust_", "", colnames(x))
  x
})

nmi_sample <- lapply(unique(spe_sub$sample_id), function(x) {
  spe_sub2 <- spe_sub[, spe_sub$sample_id == x]
  x <- compareClusters(spe_sub2, func = "NMI")
  rownames(x) <- gsub("clust_", "", rownames(x))
  colnames(x) <- gsub("clust_", "", colnames(x))
  x
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
comp_values$Target <- factor(comp_values$Target, levels = c("BayesSpace_harmony_09", "v1.5.1_ben", "v1.14.0", "v1.16.0_ben", "v1.16.0", "v1.17.0", "v1.20.2"))

## Plot ARI and NMI values by sample
## Code adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/05_ARI/02_plot_ARI.R
p_ari <- ggpaired(
    subset(comp_values, Reference == "BayesSpace_harmony_09" & ! Target %in% c("BayesSpace_harmony_09", "v1.17.0")),
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
    subset(comp_values, Reference == "BayesSpace_harmony_09" & ! Target %in% c("BayesSpace_harmony_09", "v1.17.0")),
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

pdf(file = file.path(dir_plots, "over_versions_ARI.pdf"), width = 14)
p_ari
dev.off()
pdf(file = file.path(dir_plots, "over_versions_NMI.pdf"), width = 14)
p_nmi
dev.off()

## Repeat using version 1.14.0 as the reference
p_ari_1.14.0 <- ggpaired(
    subset(comp_values, Reference == "v1.14.0" & ! Target %in% c("BayesSpace_harmony_09", "v1.17.0", "v1.14.0")),
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
    subset(comp_values, Reference == "v1.14.0" & ! Target %in% c("BayesSpace_harmony_09", "v1.17.0", "v1.14.0")),
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

pdf(file = file.path(dir_plots, "over_versions_ARI_1.14.0.pdf"), width = 14)
p_ari_1.14.0
dev.off()
pdf(file = file.path(dir_plots, "over_versions_NMI_1.14.0.pdf"), width = 14)
p_nmi_1.14.0
dev.off()

## Compare changes across versions in BayesSpace
tab <- table(
  "BayesSpace 1.14.0" = spe_sub$clust_v1.14.0,
  "BayesSpace 1.16.0" = spe_sub$clust_v1.16.0
)
addmargins(tab)
#                  BayesSpace 1.16.0
# BayesSpace 1.14.0      1      2      3      4      5      6      7      8      9    Sum
#               1     1860      5      5      7      4      1      3      1     11   1897
#               2      139   7793     65     30     30      2     32     11     50   8152
#               3        4    398  13368     34    193      0     17      7     19  14040
#               4        2      4      4  20520     43      1    206     96     13  20889
#               5        5     52    669     38  16338      0     62    278     12  17454
#               6        0      2      1     10      1   6342      9      4    259   6628
#               7        0     23     20    289    178      0  18850     26     88  19474
#               8        2      4      5    239    903      1     60  14804      9  16027
#               9        1     20     10     82     17    140    370     20   8706   9366
#               Sum   2013   8301  14147  21249  17707   6487  19609  15247   9167 113927

## Percent agreement between BayesSpace 1.14.0 and 1.16.0
round(sum(diag(tab)) / sum(tab) * 100, 2)
# [1] 95.31

## Compare changes across v1.16.0 analyses (Ben vs mine)
tab2 <- table(
  "BayesSpace 1.16.0 (by Ben)" = spe_sub$clust_v1.16.0_ben,
  "BayesSpace 1.16.0" = spe_sub$clust_v1.16.0
)
addmargins(tab2)
#                           BayesSpace 1.16.0
# BayesSpace 1.16.0 (by Ben)      1      2      3      4      5      6      7      8      9    Sum
#                        1     1764   5617    249    127    197     29    172     52    396   8603
#                        2       14   2068   3226    563   3105      1    897    509    135  10518
#                        3       36    184   9489    283   3978      0    466    327    179  14942
#                        4       55     36    167  11459    436      3   1573   1046    193  14968
#                        5       70    235    532   4851   5115      4   2930   5836    301  19874
#                        6        4     25      9     38      6   5611     55     12   1298   7058
#                        7       23     32    257   1950   1000      7  10479    403    464  14615
#                        8       11     19    147   1039   3744      0    477   6832     37  12306
#                        9       36     85     71    939    126    832   2560    230   6164  11043
#                        Sum   2013   8301  14147  21249  17707   6487  19609  15247   9167 113927

## Percent agreement between BayesSpace 1.16.0 (by Ben) and 1.16.0
round(sum(diag(tab2)) / sum(tab2) * 100, 2)
# [1] 51.77


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120) ## Makes it easier to read later
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.5.0 Patched (2025-05-21 r88220)
#  os       Rocky Linux 9.4 (Blue Onyx)
#  system   x86_64, linux-gnu
#  ui       Positron
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2026-08-14
#  pandoc   3.7.0.1 @ /jhpce/shared/community/core/conda_R/4.5/bin/pandoc
#  quarto   NA

# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-8     2024-09-12 [2] CRAN (R 4.5.0)
#  AnnotationDbi          1.70.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  AnnotationHub          3.16.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  aricode                1.1.0     2026-05-13 [1] CRAN (R 4.5.0)
#  arrow                  20.0.0    2025-05-11 [2] CRAN (R 4.5.0)
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 4.5.0)
#  attempt                0.3.1     2020-05-03 [2] CRAN (R 4.5.0)
#  backports              1.5.0     2024-05-23 [2] CRAN (R 4.5.0)
#  Banksy               * 1.4.0     2025-04-15 [1] Bioconductor 3.21 (R 4.5.0)
#  BayesSpace           * 1.17.0    2024-10-29 [2] Bioconductor 3.21 (R 4.5.0)
#  beachmat               2.24.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.5.0)
#  benchmarkme            1.0.8     2022-06-12 [2] CRAN (R 4.5.0)
#  benchmarkmeData        1.0.4     2020-04-23 [2] CRAN (R 4.5.0)
#  Biobase              * 2.68.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  BiocFileCache          2.16.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  BiocGenerics         * 0.54.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  BiocIO                 1.18.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  BiocManager            1.30.25   2024-08-28 [2] CRAN (R 4.5.0)
#  BiocNeighbors          2.2.0     2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  BiocParallel           1.42.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  BiocSingular           1.24.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  BiocVersion            3.21.1    2024-10-29 [2] Bioconductor 3.21 (R 4.5.0)
#  Biostrings             2.76.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  bit                    4.6.0     2025-03-06 [2] CRAN (R 4.5.0)
#  bit64                  4.6.0-1   2025-01-16 [2] CRAN (R 4.5.0)
#  bitops                 1.0-9     2024-10-03 [2] CRAN (R 4.5.0)
#  blob                   1.2.4     2023-03-17 [2] CRAN (R 4.5.0)
#  bluster                1.18.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  broom                  1.0.8     2025-03-28 [2] CRAN (R 4.5.0)
#  bslib                  0.9.0     2025-01-30 [2] CRAN (R 4.5.0)
#  cachem                 1.1.0     2024-05-16 [2] CRAN (R 4.5.0)
#  callr                  3.7.6     2024-03-25 [2] CRAN (R 4.5.0)
#  car                    3.1-3     2024-09-27 [2] CRAN (R 4.5.0)
#  carData                3.0-5     2022-01-06 [2] CRAN (R 4.5.0)
#  circlize               0.4.16    2024-02-20 [2] CRAN (R 4.5.0)
#  class                  7.3-23    2025-01-01 [3] CRAN (R 4.5.0)
#  cli                    3.6.5     2025-04-23 [2] CRAN (R 4.5.0)
#  clue                   0.3-66    2024-11-13 [2] CRAN (R 4.5.0)
#  cluster                2.1.8.1   2025-03-12 [3] CRAN (R 4.5.0)
#  coda                   0.19-4.1  2024-01-31 [2] CRAN (R 4.5.0)
#  codetools              0.2-20    2024-03-31 [3] CRAN (R 4.5.0)
#  colorout             * 1.3-3     2026-01-09 [1] Github (jalvesaq/colorout@64863bb)
#  colorspace             2.1-1     2024-07-26 [2] CRAN (R 4.5.0)
#  ComplexHeatmap         2.24.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  config                 0.3.2     2023-08-30 [2] CRAN (R 4.5.0)
#  cowplot                1.1.3     2024-01-22 [2] CRAN (R 4.5.0)
#  crayon                 1.5.3     2024-06-20 [2] CRAN (R 4.5.0)
#  curl                   6.2.2     2025-03-24 [2] CRAN (R 4.5.0)
#  data.table             1.17.2    2025-05-12 [2] CRAN (R 4.5.0)
#  DBI                    1.2.3     2024-06-02 [2] CRAN (R 4.5.0)
#  dbplyr                 2.5.0     2024-03-19 [2] CRAN (R 4.5.0)
#  dbscan                 1.2.5     2026-06-09 [1] CRAN (R 4.5.0)
#  DelayedArray           0.34.1    2025-04-17 [2] Bioconductor 3.21 (R 4.5.0)
#  dichromat              2.0-0.1   2022-05-02 [2] CRAN (R 4.5.0)
#  digest                 0.6.39    2025-11-19 [1] CRAN (R 4.5.0)
#  DirichletReg           0.7-2     2025-05-31 [1] CRAN (R 4.5.0)
#  doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.5.0)
#  dplyr                  1.1.4     2023-11-17 [2] CRAN (R 4.5.0)
#  dqrng                  0.4.1     2024-05-28 [2] CRAN (R 4.5.0)
#  DT                     0.33      2024-04-04 [2] CRAN (R 4.5.0)
#  edgeR                  4.6.3     2025-07-09 [1] Bioconductor 3.21 (R 4.5.0)
#  ExperimentHub          2.16.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  farver                 2.1.2     2024-05-13 [2] CRAN (R 4.5.0)
#  fastmap                1.2.0     2024-05-15 [2] CRAN (R 4.5.0)
#  filelock               1.0.3     2023-12-11 [2] CRAN (R 4.5.0)
#  foreach                1.5.2     2022-02-02 [2] CRAN (R 4.5.0)
#  Formula                1.2-5     2023-02-24 [2] CRAN (R 4.5.0)
#  generics             * 0.1.4     2025-05-09 [2] CRAN (R 4.5.0)
#  GenomeInfoDb         * 1.44.3    2025-09-21 [1] Bioconductor 3.21 (R 4.5.0)
#  GenomeInfoDbData       1.2.14    2025-05-21 [2] Bioconductor
#  GenomicAlignments      1.44.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  GenomicRanges        * 1.60.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  GetoptLong             1.0.5     2020-12-15 [2] CRAN (R 4.5.0)
#  ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.5.0)
#  ggplot2              * 4.0.3     2026-04-22 [1] CRAN (R 4.5.0)
#  ggpubr               * 0.6.0     2023-02-10 [2] CRAN (R 4.5.0)
#  ggrepel              * 0.9.6     2024-09-07 [2] CRAN (R 4.5.0)
#  ggsignif               0.6.4     2022-10-13 [2] CRAN (R 4.5.0)
#  GlobalOptions          0.1.2     2020-06-10 [2] CRAN (R 4.5.0)
#  glue                   1.8.0     2024-09-30 [2] CRAN (R 4.5.0)
#  golem                  0.5.1     2024-08-27 [2] CRAN (R 4.5.0)
#  gridExtra              2.3       2017-09-09 [2] CRAN (R 4.5.0)
#  gtable                 0.3.6     2024-10-25 [2] CRAN (R 4.5.0)
#  here                 * 1.0.2     2025-09-15 [1] CRAN (R 4.5.0)
#  htmltools              0.5.9     2025-12-04 [1] CRAN (R 4.5.0)
#  htmlwidgets            1.6.4     2023-12-06 [2] CRAN (R 4.5.0)
#  httpuv                 1.6.16    2025-04-16 [2] CRAN (R 4.5.0)
#  httr                   1.4.7     2023-08-15 [2] CRAN (R 4.5.0)
#  igraph                 2.1.4     2025-01-23 [2] CRAN (R 4.5.0)
#  IRanges              * 2.42.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.5.0)
#  iterators              1.0.14    2022-02-05 [2] CRAN (R 4.5.0)
#  jquerylib              0.1.4     2021-04-26 [2] CRAN (R 4.5.0)
#  jsonlite               2.0.0     2025-03-27 [2] CRAN (R 4.5.0)
#  KEGGREST               1.48.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  labeling               0.4.3     2023-08-29 [2] CRAN (R 4.5.0)
#  later                  1.4.8     2026-03-05 [1] CRAN (R 4.5.0)
#  lattice                0.22-7    2025-04-02 [3] CRAN (R 4.5.0)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 4.5.0)
#  leidenAlg              1.1.8     2026-05-31 [1] CRAN (R 4.5.0)
#  lifecycle              1.0.5     2026-01-08 [1] CRAN (R 4.5.0)
#  limma                  3.64.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  locfit                 1.5-9.12  2025-03-05 [2] CRAN (R 4.5.0)
#  magick                 2.8.6     2025-03-23 [2] CRAN (R 4.5.0)
#  magrittr               2.0.5     2026-04-04 [1] CRAN (R 4.5.0)
#  Matrix                 1.7-3     2025-03-11 [3] CRAN (R 4.5.0)
#  MatrixGenerics       * 1.20.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  matrixStats          * 1.5.0     2025-01-07 [2] CRAN (R 4.5.0)
#  maxLik                 1.5-2.2   2025-12-29 [1] CRAN (R 4.5.0)
#  mclust                 6.1.1     2024-04-29 [2] CRAN (R 4.5.0)
#  memoise                2.0.1     2021-11-26 [2] CRAN (R 4.5.0)
#  metapod                1.16.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  microbenchmark         1.5.0     2024-09-04 [2] CRAN (R 4.5.0)
#  mime                   0.13      2025-03-17 [2] CRAN (R 4.5.0)
#  miscTools              0.6-28    2023-05-03 [2] CRAN (R 4.5.0)
#  otel                   0.2.0     2025-08-29 [1] CRAN (R 4.5.0)
#  pak                    0.11.1    2026-07-22 [1] CRAN (R 4.5.0)
#  paletteer              1.7.0     2026-01-08 [1] CRAN (R 4.5.0)
#  pillar                 1.10.2    2025-04-05 [2] CRAN (R 4.5.0)
#  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.5.0)
#  plotly                 4.10.4    2024-01-13 [2] CRAN (R 4.5.0)
#  png                    0.1-8     2022-11-29 [2] CRAN (R 4.5.0)
#  processx               3.8.6     2025-02-21 [2] CRAN (R 4.5.0)
#  promises               1.5.0     2025-11-01 [1] CRAN (R 4.5.0)
#  ps                     1.9.1     2025-04-12 [2] CRAN (R 4.5.0)
#  purrr                  1.0.4     2025-02-05 [2] CRAN (R 4.5.0)
#  R6                     2.6.1     2025-02-15 [2] CRAN (R 4.5.0)
#  ragg                   1.4.0     2025-04-10 [2] CRAN (R 4.5.0)
#  rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.5.0)
#  RColorBrewer         * 1.1-3     2022-04-03 [2] CRAN (R 4.5.0)
#  Rcpp                   1.1.2     2026-07-05 [1] CRAN (R 4.5.0)
#  RcppHungarian          0.3       2023-09-05 [1] CRAN (R 4.5.0)
#  RCurl                  1.98-1.17 2025-03-22 [2] CRAN (R 4.5.0)
#  rematch2               2.1.2     2020-05-01 [2] CRAN (R 4.5.0)
#  restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.5.0)
#  rhdf5                  2.52.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  rhdf5filters           1.20.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  Rhdf5lib               1.30.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  rjson                  0.2.23    2024-09-16 [2] CRAN (R 4.5.0)
#  rlang                  1.3.0     2026-07-05 [1] CRAN (R 4.5.0)
#  rprojroot              2.1.1     2025-08-26 [1] CRAN (R 4.5.0)
#  Rsamtools              2.24.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  RSQLite                2.3.11    2025-05-04 [2] CRAN (R 4.5.0)
#  rstatix                0.7.2     2023-02-01 [2] CRAN (R 4.5.0)
#  rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.5.0)
#  rtracklayer            1.68.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  S4Arrays               1.8.0     2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  S4Vectors            * 0.46.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  S7                     0.2.2     2026-04-22 [1] CRAN (R 4.5.0)
#  sandwich               3.1-1     2024-09-15 [2] CRAN (R 4.5.0)
#  sass                   0.4.10    2025-04-11 [2] CRAN (R 4.5.0)
#  ScaledMatrix           1.16.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  scales                 1.4.0     2025-04-24 [2] CRAN (R 4.5.0)
#  scater                 1.36.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  sccore                 1.0.7     2026-04-06 [1] CRAN (R 4.5.0)
#  scran                  1.36.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  scuttle                1.18.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  sessioninfo          * 1.2.3     2025-02-05 [2] CRAN (R 4.5.0)
#  shape                  1.4.6.1   2024-02-23 [2] CRAN (R 4.5.0)
#  shiny                  1.10.0    2024-12-14 [2] CRAN (R 4.5.0)
#  shinyWidgets           0.9.0     2025-02-21 [2] CRAN (R 4.5.0)
#  SingleCellExperiment * 1.30.1    2025-05-07 [2] Bioconductor 3.21 (R 4.5.0)
#  SparseArray            1.8.0     2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  SpatialExperiment    * 1.18.1    2025-05-11 [2] Bioconductor 3.21 (R 4.5.0)
#  spatialLIBD          * 1.20.1    2025-05-01 [2] Bioconductor 3.21 (R 4.5.0)
#  statmod                1.5.0     2023-01-06 [2] CRAN (R 4.5.0)
#  stringi                1.8.7     2025-03-27 [2] CRAN (R 4.5.0)
#  stringr                1.5.1     2023-11-14 [2] CRAN (R 4.5.0)
#  SummarizedExperiment * 1.38.1    2025-04-30 [2] Bioconductor 3.21 (R 4.5.0)
#  systemfonts            1.3.2     2026-03-05 [1] CRAN (R 4.5.0)
#  textshaping            1.0.1     2025-05-01 [2] CRAN (R 4.5.0)
#  tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.5.0)
#  tidyr                * 1.3.1     2024-01-24 [2] CRAN (R 4.5.0)
#  tidyselect             1.2.1     2024-03-11 [2] CRAN (R 4.5.0)
#  UCSC.utils             1.4.0     2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  uwot                   0.2.4     2025-11-10 [1] CRAN (R 4.5.0)
#  vctrs                  0.6.5     2023-12-01 [2] CRAN (R 4.5.0)
#  vipor                  0.4.7     2023-12-18 [2] CRAN (R 4.5.0)
#  viridis                0.6.5     2024-01-29 [2] CRAN (R 4.5.0)
#  viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.5.0)
#  withr                  3.0.2     2024-10-28 [2] CRAN (R 4.5.0)
#  xgboost                3.2.1.1   2026-03-18 [1] CRAN (R 4.5.0)
#  XML                    3.99-0.18 2025-01-01 [2] CRAN (R 4.5.0)
#  xtable                 1.8-4     2019-04-21 [2] CRAN (R 4.5.0)
#  XVector                0.48.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
#  yaml                   2.3.12    2025-12-10 [1] CRAN (R 4.5.0)
#  zoo                    1.8-14    2025-04-10 [2] CRAN (R 4.5.0)

#  [1] /users/lcollado/R/4.5
#  [2] /jhpce/shared/community/core/conda_R/4.5/R/lib64/R/site-library
#  [3] /jhpce/shared/community/core/conda_R/4.5/R/lib64/R/library
#  * ── Packages attached to the search path.

# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
