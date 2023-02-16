library("here")
library("sessioninfo")

dir_rdata <- here::here(
    "processed-data", "rdata", "spe",
    "09_position_differential_expression"
)
stopifnot(file.exists(dir_rdata))

## Load Sp09 DE results
load(here("code", "deploy_app_k09", "sig_genes_subset_k09.Rdata"),
    verbose = TRUE
)
sig_domain_09 <- sig_genes
rm(sig_genes)

## Load Sp16 DE results
load(here("code", "deploy_app_k16", "sig_genes_subset_k16.Rdata"),
    verbose = TRUE
)
sig_domain_16 <- sig_genes
rm(sig_genes)

## Combine Sp09 and Sp16 results
sig_domain <- rbind(
    cbind(sig_domain_09, "spatial_domain_resolution" = "Sp09"),
    cbind(sig_domain_16, "spatial_domain_resolution" = "Sp16")
)

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

## Clean up for exporting
clean_sig <- function(sig_genes) {
    sig_genes <- subset(sig_genes, fdr < 0.05)
    sig_genes$gene_index <- NULL
    sig_genes$results <- NULL
    return(as.data.frame(sig_genes))
}

sig_domain_fdr <- clean_sig(sig_domain)
head(sig_domain_fdr)
dim(sig_domain_fdr)
# [1] 1029505       9

## Split domain results into 3 tables
write.csv(
    subset(sig_domain_fdr, model_type == "enrichment"),
    file.path(dir_rdata, "sig_genes_FDR5perc_enrichment.csv"),
    row.names = FALSE
)
write.csv(
    subset(sig_domain_fdr, model_type == "anova"),
    file.path(dir_rdata, "sig_genes_FDR5perc_anova.csv"),
    row.names = FALSE
)
write.csv(
    subset(sig_domain_fdr, model_type == "pairwise"),
    file.path(dir_rdata, "sig_genes_FDR5perc_pairwise.csv"),
    row.names = FALSE
)

## Save position results into 1 table since it's not too long
dim(clean_sig(sig_position))
# [1] 7272    8
write.csv(
    clean_sig(sig_position),
    file.path(dir_rdata, "sig_genes_FDR5perc_position.csv"),
    row.names = FALSE
)

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
#  date     2023-02-14
#  rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
#  pandoc   2.17.1.1 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date (UTC) lib source
#  brio          1.1.3   2021-11-30 [1] CRAN (R 4.2.0)
#  cachem        1.0.6   2021-08-19 [1] CRAN (R 4.2.0)
#  callr         3.7.3   2022-11-02 [1] CRAN (R 4.2.2)
#  cli           3.6.0   2023-01-09 [1] CRAN (R 4.2.0)
#  colorout      1.2-2   2022-03-01 [1] Github (jalvesaq/colorout@79931fd)
#  crayon        1.5.2   2022-09-29 [1] CRAN (R 4.2.0)
#  data.table    1.14.6  2022-11-16 [1] CRAN (R 4.2.0)
#  devtools    * 2.4.5   2022-10-11 [1] CRAN (R 4.2.0)
#  digest        0.6.31  2022-12-11 [1] CRAN (R 4.2.0)
#  ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.2.0)
#  fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.2.0)
#  fs            1.6.0   2023-01-23 [1] CRAN (R 4.2.0)
#  generics      0.1.3   2022-07-05 [1] CRAN (R 4.2.0)
#  glue          1.6.2   2022-02-24 [1] CRAN (R 4.2.0)
#  here        * 1.0.1   2020-12-13 [1] CRAN (R 4.2.0)
#  hms           1.1.2   2022-08-19 [1] CRAN (R 4.2.0)
#  htmltools     0.5.4   2022-12-07 [1] CRAN (R 4.2.0)
#  htmlwidgets   1.6.1   2023-01-07 [1] CRAN (R 4.2.0)
#  httpuv        1.6.8   2023-01-12 [1] CRAN (R 4.2.0)
#  later         1.3.0   2021-08-18 [1] CRAN (R 4.2.0)
#  lifecycle     1.0.3   2022-10-07 [1] CRAN (R 4.2.1)
#  lubridate     1.9.1   2023-01-24 [1] CRAN (R 4.2.0)
#  magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.2.0)
#  memoise       2.0.1   2021-11-26 [1] CRAN (R 4.2.0)
#  mime          0.12    2021-09-28 [1] CRAN (R 4.2.0)
#  miniUI        0.1.1.1 2018-05-18 [1] CRAN (R 4.2.0)
#  pkgbuild      1.4.0   2022-11-27 [1] CRAN (R 4.2.2)
#  pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.2.0)
#  pkgload       1.3.2   2022-11-16 [1] CRAN (R 4.2.2)
#  prettyunits   1.1.1   2020-01-24 [1] CRAN (R 4.2.0)
#  processx      3.8.0   2022-10-26 [1] CRAN (R 4.2.0)
#  profvis       0.3.7   2020-11-02 [1] CRAN (R 4.2.0)
#  promises      1.2.0.1 2021-02-11 [1] CRAN (R 4.2.0)
#  prompt        1.0.1   2022-03-01 [1] Github (gaborcsardi/prompt@7ef0f2e)
#  ps            1.7.2   2022-10-26 [1] CRAN (R 4.2.0)
#  purrr         1.0.1   2023-01-10 [1] CRAN (R 4.2.0)
#  R6            2.5.1   2021-08-19 [1] CRAN (R 4.2.0)
#  Rcpp          1.0.10  2023-01-22 [1] CRAN (R 4.2.0)
#  remotes       2.4.2   2021-11-30 [1] CRAN (R 4.2.0)
#  rlang         1.0.6   2022-09-24 [1] CRAN (R 4.2.0)
#  rprojroot     2.0.3   2022-04-02 [1] CRAN (R 4.2.0)
#  rsthemes      0.3.1   2022-03-01 [1] Github (gadenbuie/rsthemes@bbe73ca)
#  rstudioapi    0.14    2022-08-22 [1] CRAN (R 4.2.0)
#  sessioninfo * 1.2.2   2021-12-06 [1] CRAN (R 4.2.0)
#  shiny         1.7.4   2022-12-15 [1] CRAN (R 4.2.2)
#  stringi       1.7.12  2023-01-11 [1] CRAN (R 4.2.0)
#  stringr       1.5.0   2022-12-02 [1] CRAN (R 4.2.0)
#  suncalc       0.5.1   2022-09-29 [1] CRAN (R 4.2.0)
#  testthat    * 3.1.6   2022-12-09 [1] CRAN (R 4.2.0)
#  timechange    0.2.0   2023-01-11 [1] CRAN (R 4.2.0)
#  urlchecker    1.0.1   2021-11-30 [1] CRAN (R 4.2.0)
#  usethis     * 2.1.6   2022-05-25 [1] CRAN (R 4.2.0)
#  vctrs         0.5.2   2023-01-23 [1] CRAN (R 4.2.0)
#  xtable        1.8-4   2019-04-21 [1] CRAN (R 4.2.0)
#
#  [1] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
