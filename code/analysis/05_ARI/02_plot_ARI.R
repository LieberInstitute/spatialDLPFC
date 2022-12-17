library("ggpubr")
library("ggrepel")
library("here")
library("RColorBrewer")
library("sessioninfo")

## Load the ARI results
load(
    here(
        "processed-data",
        "rdata",
        "pilot_dlpfc_data",
        "05_ARI",
        "pilot_ari_clustering_across.Rdata"
    ),
    verbose = TRUE
)

## Fix the sample_id
# load(here("processed-data", "rdata", "pilot_dlpfc_data", "spe_pilot_bayesSpace_batch_corr_sampleID.Rdata"), verbose = TRUE)
# sample_ids <- unique(spe$sample_id)
# paste0(sample_ids, collapse = '", "')
sample_ids <-
    c(
        "151507",
        "151508",
        "151509",
        "151510",
        "151669",
        "151670",
        "151671",
        "151672",
        "151673",
        "151674",
        "151675",
        "151676"
    )
ari.df.long$sample_id <- rep(sample_ids, 5)

## Fix capitalization of Graph-based to Graph-Based
levels(ari.df.long$method) <-
    gsub(
        "Graph-based\\(BC\\)",
        "Graph-Based(BC)",
        levels(ari.df.long$method)
    )

pdf(here(
    "plots",
    "05_ARI",
    "ggboxplot_pilot_data_ARI_clustering_across.pdf"
))
set.seed(20221216)
ggboxplot(
    ari.df.long,
    x = "method",
    y = "ari",
    fill = "method",
    palette = brewer.pal("Paired", n = 10)[c(1:4, 10)],
    add = "jitter",
    repel = TRUE,
    font.label = list(size = 10),
    legend = "none",
    ggtheme = theme_pubr(base_size = 30),
    ylab = "Adjusted Rand Index",
    xlab = "Clustering Method",
    size = 1,
    label = "sample_id"
) +
    font("xy.text", size = 11) +
    font("xlab", size = 16) +
    font("ylab", size = 16)
dev.off()


## Make a new version with paired lines
pdf(
    here(
        "plots",
        "05_ARI",
        "ggpaired_pilot_data_ARI_clustering_across.pdf"
    ),
    width = 12,
    height = 9
)
ggpaired(
    ari.df.long,
    x = "method",
    y = "ari",
    id = "sample_id",
    xlab = "Clustering Method",
    ylab = "Adjusted Rand Index",
    palette = brewer.pal("Paired", n = 10)[c(1:4, 10)],
    fill = "method",
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
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.2 (2022-10-31)
#  os       macOS Ventura 13.0.1
#  system   aarch64, darwin20
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2022-12-17
#  pandoc   2.17.1.1 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────
#  package      * version date (UTC) lib source
#  abind          1.4-5   2016-07-21 [1] CRAN (R 4.2.0)
#  assertthat     0.2.1   2019-03-21 [1] CRAN (R 4.2.0)
#  backports      1.4.1   2021-12-13 [1] CRAN (R 4.2.0)
#  broom          1.0.2   2022-12-15 [1] CRAN (R 4.2.2)
#  car            3.1-1   2022-10-19 [1] CRAN (R 4.2.0)
#  carData        3.0-5   2022-01-06 [1] CRAN (R 4.2.0)
#  cli            3.4.1   2022-09-23 [1] CRAN (R 4.2.0)
#  colorspace     2.0-3   2022-02-21 [1] CRAN (R 4.2.0)
#  DBI            1.1.3   2022-06-18 [1] CRAN (R 4.2.0)
#  dplyr          1.0.10  2022-09-01 [1] CRAN (R 4.2.1)
#  fansi          1.0.3   2022-03-24 [1] CRAN (R 4.2.0)
#  generics       0.1.3   2022-07-05 [1] CRAN (R 4.2.0)
#  ggplot2      * 3.4.0   2022-11-04 [1] CRAN (R 4.2.0)
#  ggpubr       * 0.5.0   2022-11-16 [1] CRAN (R 4.2.2)
#  ggrepel      * 0.9.2   2022-11-06 [1] CRAN (R 4.2.0)
#  ggsignif       0.6.4   2022-10-13 [1] CRAN (R 4.2.0)
#  glue           1.6.2   2022-02-24 [1] CRAN (R 4.2.0)
#  gtable         0.3.1   2022-09-01 [1] CRAN (R 4.2.1)
#  here         * 1.0.1   2020-12-13 [1] CRAN (R 4.2.0)
#  lifecycle      1.0.3   2022-10-07 [1] CRAN (R 4.2.1)
#  magrittr       2.0.3   2022-03-30 [1] CRAN (R 4.2.0)
#  munsell        0.5.0   2018-06-12 [1] CRAN (R 4.2.0)
#  pillar         1.8.1   2022-08-19 [1] CRAN (R 4.2.0)
#  pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.2.0)
#  purrr          0.3.5   2022-10-06 [1] CRAN (R 4.2.1)
#  R6             2.5.1   2021-08-19 [1] CRAN (R 4.2.0)
#  RColorBrewer * 1.1-3   2022-04-03 [1] CRAN (R 4.2.0)
#  Rcpp           1.0.9   2022-07-08 [1] CRAN (R 4.2.1)
#  rlang          1.0.6   2022-09-24 [1] CRAN (R 4.2.0)
#  rprojroot      2.0.3   2022-04-02 [1] CRAN (R 4.2.0)
#  rstatix        0.7.1   2022-11-09 [1] CRAN (R 4.2.0)
#  scales         1.2.1   2022-08-20 [1] CRAN (R 4.2.0)
#  sessioninfo  * 1.2.2   2021-12-06 [1] CRAN (R 4.2.0)
#  tibble         3.1.8   2022-07-22 [1] CRAN (R 4.2.1)
#  tidyr          1.2.1   2022-09-08 [1] CRAN (R 4.2.0)
#  tidyselect     1.2.0   2022-10-10 [1] CRAN (R 4.2.0)
#  utf8           1.2.2   2021-07-24 [1] CRAN (R 4.2.0)
#  vctrs          0.5.1   2022-11-16 [1] CRAN (R 4.2.2)
#  withr          2.5.0   2022-03-03 [1] CRAN (R 4.2.0)
#
#  [1] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────
