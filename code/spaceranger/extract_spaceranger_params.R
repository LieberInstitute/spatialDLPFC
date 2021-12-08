library("here")
library("stringr")
library("sessioninfo")

## Locate the 30 _invocation files
round1 <-
    list.files(here("processed-data", "NextSeq"),
        pattern = "DLPFC",
        full.names = TRUE)
round2_4 <-
    list.files(list.files(
        here("processed-data", "NextSeq"),
        pattern = "Round",
        full.names = TRUE
    ),
        full.names = TRUE)
invocation_files <- file.path(c(round1, round2_4), "_invocation")
invocation_files <- invocation_files[file.exists(invocation_files)]
stopifnot(length(invocation_files) == length(unique(invocation_files)))
stopifnot(length(invocation_files) == 30)

names(invocation_files) <- basename(dirname(invocation_files))

update_path <- function(x) {
    gsub(
        "Images/Liebert",
        "Images/round1/Liebert",
        gsub(
            "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC",
            "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data",
            x
        )
    )
}

## Extract relevant SpaceRanger parameters
df <- do.call(rbind, lapply(invocation_files, function(i) {
    x <- readLines(i)

    data.frame(
        slide_serial_capture_area = str_match(
            str_subset(x, "slide_serial_capture_area"),
            '\\"([:graph:]+)\\"'
        )[, 2],
        loupe_alignment_file = update_path(str_match(
            str_subset(x, "loupe_alignment_file"),
            '\\"([:graph:]+)\\"'
        )[, 2]),
        tissue_image_paths = update_path(str_match(
            str_subset(x, "tissue_image_paths"),
            '\\"([:graph:]+)\\"'
        )[, 2]),
        read_paths = paste(update_path(
            str_match(
                str_subset(x, "read_path"),
                ':[:space:]*\\"([:graph:]+)\\"'
            )[, 2]
        ), collapse = ",")
    )

}))

## Use the directory names as the sample_ids to match the output
df$sample_id <- rownames(df)
rownames(df) <- NULL

## Check that all files are correct
stopifnot(all(file.exists(df$loupe_alignment_file)))
stopifnot(all(file.exists(df$tissue_image_paths)))
stopifnot(all(file.exists(unlist(
    str_split(df$read_paths, ",")
))))

## Extract capture area and slide serial
df$capture_area <- gsub(".*-", "", df$slide_serial_capture_area)
df$slide_serial <-
    str_match(df$slide_serial_capture_area, "^[:alnum:]+-[:alnum:]+")[, 1]

write.table(
    df[, c(
        "sample_id",
        "slide_serial",
        "capture_area",
        "tissue_image_paths",
        "loupe_alignment_file",
        "read_paths"
    )],
    file = here("code", "spaceranger", "spaceranger_parameters.txt"),
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t",
    quote = FALSE
)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.2 Patched (2021-11-04 r81138)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
# 2021-12-01 [2] CRAN (R 4.1.2)
#  dplyr         1.0.7   2021-06-18 [2] CRAN (R 4.1.0)
#  ellipsis      0.3.2   2021-04-29 [2] CRAN (R 4.1.0)
#  fansi         0.5.0   2021-05-25 [2] CRAN (R 4.1.0)
#  fastmap       1.1.0   2021-01-25 [2] CRAN (R 4.1.0)
#  generics      0.1.1   2021-10-25 [2] CRAN (R 4.1.2)
#  ggplot2       3.3.5   2021-06-25 [2] CRAN (R 4.1.0)
#  glue          1.5.1   2021-11-30 [2] CRAN (R 4.1.2)
#  gtable        0.3.0   2019-03-25 [2] CRAN (R 4.1.0)
#  here        * 1.0.1   2020-12-13 [1] CRAN (R 4.1.2)
#  htmltools     0.5.2   2021-08-25 [2] CRAN (R 4.1.2)
#  htmlwidgets   1.5.4   2021-09-08 [2] CRAN (R 4.1.2)
#  httpuv        1.6.3   2021-09-09 [2] CRAN (R 4.1.2)
#  jsonlite      1.7.2   2020-12-09 [2] CRAN (R 4.1.0)
#  later         1.3.0   2021-08-18 [2] CRAN (R 4.1.2)
#  lattice       0.20-45 2021-09-22 [3] CRAN (R 4.1.2)
#  lifecycle     1.0.1   2021-09-24 [2] CRAN (R 4.1.2)
#  magrittr      2.0.1   2020-11-17 [2] CRAN (R 4.1.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 4.1.0)
#  pillar        1.6.4   2021-10-18 [2] CRAN (R 4.1.2)
#  pkgconfig     2.0.3   2019-09-22 [2] CRAN (R 4.1.0)
#  png           0.1-7   2013-12-03 [2] CRAN (R 4.1.0)
#  promises      1.2.0.1 2021-02-11 [2] CRAN (R 4.1.0)
#  purrr         0.3.4   2020-04-17 [2] CRAN (R 4.1.0)
#  R6            2.5.1   2021-08-19 [2] CRAN (R 4.1.2)
#  Rcpp          1.0.7   2021-07-07 [2] CRAN (R 4.1.0)
#  rlang         0.4.12  2021-10-18 [2] CRAN (R 4.1.2)
#  rmote         0.3.4   2021-11-02 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot     2.0.2   2020-11-15 [2] CRAN (R 4.1.0)
#  scales        1.1.1   2020-05-11 [2] CRAN (R 4.1.0)
#  servr         0.24    2021-11-16 [1] CRAN (R 4.1.2)
#  sessioninfo * 1.2.2   2021-12-06 [2] CRAN (R 4.1.2)
#  stringi       1.7.6   2021-11-29 [2] CRAN (R 4.1.2)
#  stringr     * 1.4.0   2019-02-10 [2] CRAN (R 4.1.0)
#  tibble        3.1.6   2021-11-07 [2] CRAN (R 4.1.2)
#  tidyselect    1.1.1   2021-04-30 [2] CRAN (R 4.1.0)
#  utf8          1.2.2   2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs         0.3.8   2021-04-29 [2] CRAN (R 4.1.0)
#  xfun          0.28    2021-11-04 [2] CRAN (R 4.1.2)
#
#  [1] /users/lcollado/R/4.1.x
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
