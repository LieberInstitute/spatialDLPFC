## Automatically style the code in this script:
# styler::style_file("01_build_SPE.R", transformers = biocthis::bioc_style())

## This script requires R 4.1
# module load conda_R/devel

## utils
library("here")
library("sessioninfo")
library("pryr")

## reading the data
library("SpatialExperiment")
library("rtracklayer")

## vis
library("spatialLIBD")

## analysis
library("scran")
library("scater")
library("BiocParallel")


## Define some info for the samples
sample_info <- data.frame(
    sample_id = c(
        "DLPFC_Br2743_ant_manual_alignment",
        "DLPFC_Br2743_mid_manual_alignment",
        "DLPFC_Br2743_post_manual_alignment",
        "DLPFC_Br3942_ant_manual_alignment",
        "DLPFC_Br3942_mid_manual_alignment",
        "DLPFC_Br3942_post_manual_alignment",
        "DLPFC_Br6423_ant_manual_alignment",
        "DLPFC_Br6423_mid_manual_alignment",
        "DLPFC_Br6423_post_manual_alignment",
        "DLPFC_Br8492_ant_manual_alignment",
        "DLPFC_Br8492_mid_manual_alignment",
        "DLPFC_Br8492_post_manual_alignment"
    ),
    subjects = rep(
        c("Br2743", "Br3942", "Br6423", "Br8492"),
        each = 3
    ),
    regions = rep(
        c("anterior", "middle", "posterior"),
        4
    )
)
sample_info$sample_path <- file.path(
    here::here("outputs", "NextSeq"),
    sample_info$sample_id,
    "outs"
)
stopifnot(all(file.exists(sample_info$sample_path)))

## Re-shuffle to match default spatialLIBD grids
sample_info <- sample_info[rep(c(1, 4, 7, 10), 3) + rep(0:2, each = 4), ]

## For testing the code with a subset of the data
# if(FALSE) {
#     sample_info <- sample_info[1:2, ]
# }

## Read the data
Sys.time()
spe <- read10xVisium(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = "lowres",
    load = TRUE
)
Sys.time()
## About 4.5 minutes
# [1] "2021-02-12 18:57:04 EST"
# [1] "2021-02-12 19:01:27 EST"


## Add some information used by spatialLIBD
spe$key <- paste0(spe$Barcode, '_', spe$sample_id)
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)

## Add the experimental information
spe$subject <- sample_info$subjects[match(spe$sample_id, sample_info$sample_id)]
spe$region <- sample_info$regions[match(spe$sample_id, sample_info$sample_id)]


## Read in the gene information from the annotation GTF file
gtf <-
    rtracklayer::import(
        "/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(spe), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SPE object
rowRanges(spe) <- gtf[match_genes]

## Add information used by spatialLIBD
rowData(spe)$gene_search <- paste0(rowData(spe)$gene_id, "; ", rowData(spe)$gene_name)
is_mito <- which(seqnames(spe) == "chrM")
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi


## Add number of cells per spot
spe$cell_count <- NA
## Currently we don't have that information

# issue 11
## Read in the number of cells per spot
# cells <-
#     do.call(rbind, lapply(dir("Histology"), function(sampleid) {
#         x <-
#             read.csv(file.path("Histology", sampleid, "tissue_spot_counts.csv"))
#         x$key <- paste0(sampleid, "_", x$barcode)
#         return(x[, c("key", "count")])
#     }))


## Simplify sample_ids for plotting
spe$sample_id <- gsub("DLPFC_|_manual_alignment", "", spe$sample_id)

## Create two versions: one with and one without filtering by tissue spot
spe_raw <- spe
spe <- spe_raw[, which(inTissue(spe_raw))]
pryr::object_size(spe_raw)
# 971 MB
pryr::object_size(spe)
# 940 MB

dim(spe_raw)
# [1] 36601 59904
dim(spe)
# [1] 36601 50006

## Save the raw version now
dir.create(here::here("rdata"), showWarnings = FALSE)
dir.create(here::here("rdata", "spe"), showWarnings = FALSE)
save(spe_raw, file = here::here("rdata", "spe", "spe_raw.Rdata"))


## Inspect in vs outside of tissue
spe_raw$overlaps_tissue <- factor(ifelse(inTissue(spe_raw), "in", "out"))
vis_grid_clus(
    spe = spe_raw,
    clustervar = "overlaps_tissue",
    pdf = here::here("plots", "in_tissue_grid.pdf"),
    sort_clust = FALSE,
    colors = c("in" = "grey90", "out" = "orange")
)

summary(spe_raw$sum_umi[!inTissue(spe_raw)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     0.0   137.0   215.0   314.5   340.0 11177.0
head(table(spe_raw$sum_umi[!inTissue(spe_raw)]))
# 0 1 2 3 4 5
# 2 4 8 7 8 3
vis_grid_gene(
    spe = spe_raw[, which(!inTissue(spe_raw))],
    geneid = "sum_umi",
    pdf = here::here("plots", "out_tissue_sum_umi.pdf"),
    assayname = "counts"
)

summary(spe_raw$sum_gene[!inTissue(spe_raw)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0   112.0   173.0   237.1   269.0  4533.0
vis_grid_gene(
    spe = spe_raw[, which(!inTissue(spe_raw))],
    geneid = "sum_gene",
    pdf = here::here("plots", "out_tissue_sum_gene.pdf"),
    assayname = "counts"
)

summary(spe_raw$expr_chrM_ratio[!inTissue(spe_raw)])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.0000  0.1319  0.1662  0.1717  0.2061  0.5455       2
vis_grid_gene(
    spe = spe_raw[, which(!inTissue(spe_raw))],
    geneid = "expr_chrM_ratio",
    pdf = here::here("plots", "out_tissue_expr_chrM_ratio.pdf"),
    assayname = "counts"
)


vis_grid_gene(
    spe = spe,
    geneid = "sum_umi",
    pdf = here::here("plots", "in_tissue_sum_umi.pdf"),
    assayname = "counts"
)

vis_grid_gene(
    spe = spe,
    geneid = "sum_umi",
    pdf = here::here("plots", "in_tissue_sum_umi.pdf"),
    assayname = "counts"
)



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
