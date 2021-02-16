## Automatically style the code in this script:
styler::style_file(
    here::here("analysis", "01_build_SPE.R"),
    transformers = biocthis::bioc_style()
)

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
    ),
    sex = rep(
        c("M", "M", "M", "F"),
        each = 3
    ),
    age = rep(
        c(61.54, 47.53, 51.73, 53.40),
        each = 3
    ),
    diagnosis = "Control"
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
## About 3-9 minutes (depending on JHPCE load)
# [1] "2021-02-15 14:00:47 EST"
# [1] "2021-02-15 14:08:45 EST"


## Add some information used by spatialLIBD
spe$key <- paste0(spe$Barcode, '_', spe$sample_id)
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)

## Add the experimental information
### Equivalent to:
## colData(spe)$subject <- ...
spe$subject <- sample_info$subjects[match(spe$sample_id, sample_info$sample_id)]
spe$region <- sample_info$regions[match(spe$sample_id, sample_info$sample_id)]
spe$sex <- sample_info$sex[match(spe$sample_id, sample_info$sample_id)]
spe$age <- sample_info$age[match(spe$sample_id, sample_info$sample_id)]
spe$diagnosis <- sample_info$diagnosis[match(spe$sample_id, sample_info$sample_id)]


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

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 9588
length(no_expr) / nrow(spe) * 100
# [1] 26.19601
spe <- spe[-no_expr, ]

## Create two versions: one with and one without filtering by tissue spot
spe_raw <- spe
spe <- spe_raw[, which(inTissue(spe_raw))]
pryr::object_size(spe_raw)
# 970 MB
pryr::object_size(spe)
# 939 MB

dim(spe_raw)
# [1] 27013 59904
dim(spe)
# [1] 27013 50006

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


summary(spe$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    0    1319    2169    2462    3264   15870
vis_grid_gene(
    spe = spe,
    geneid = "sum_umi",
    pdf = here::here("plots", "in_tissue_sum_umi.pdf"),
    assayname = "counts"
)

summary(spe$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    0     927    1399    1470    1938    5529
vis_grid_gene(
    spe = spe,
    geneid = "sum_gene",
    pdf = here::here("plots", "in_tissue_sum_gene.pdf"),
    assayname = "counts"
)

summary(spe$expr_chrM_ratio)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.00000 0.06035 0.08099 0.09029 0.10935 1.00000       7
vis_grid_gene(
    spe = spe,
    geneid = "expr_chrM_ratio",
    pdf = here::here("plots", "in_tissue_expr_chrM_ratio.pdf"),
    assayname = "counts"
)

vis_grid_gene(
    spe = spe_raw,
    geneid = "sum_umi",
    pdf = here::here("plots", "all_sum_umi.pdf"),
    assayname = "counts"
)
vis_grid_gene(
    spe = spe_raw,
    geneid = "sum_gene",
    pdf = here::here("plots", "all_sum_gene.pdf"),
    assayname = "counts"
)
vis_grid_gene(
    spe = spe_raw,
    geneid = "expr_chrM_ratio",
    pdf = here::here("plots", "all_expr_chrM_ratio.pdf"),
    assayname = "counts"
)



## Quality control (scran)

## Remove spots without counts first
spe <- spe[, -which(colSums(counts(spe)) == 0)]
ncol(spe)
# [1] 49999

qcstats <- perCellQCMetrics(spe, subsets = list(
    Mito = which(seqnames(spe) == "chrM")
))
qcfilter <- quickPerCellQC(qcstats, sub.fields="subsets_Mito_percent")
colSums(as.matrix(qcfilter))
# low_lib_size            low_n_features high_subsets_Mito_percent                   discard
#         2774                      3055                      1763                      4389

## Prior to dropping spots with 0 counts and checking for high chrM,
## this was the output:
# low_lib_size low_n_features        discard
#         2781           3062           3062

spe$scran_discard <-
    factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
spe$scran_low_lib_size <-
    factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
spe$scran_low_n_features <-
    factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
spe$scran_high_subsets_Mito_percent <-
    factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

for(i in colnames(qcfilter)) {
    vis_grid_clus(
        spe = spe,
        clustervar = paste0("scran_", i),
        pdf = here::here("plots", paste0("scran_", i, ".pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange")
    )
}

## Find quick clusters
set.seed(20191112)
Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(4),
    block = spe$sample_id,
    block.BPPARAM = MulticoreParam(4)
)
Sys.time()


Sys.time()
## Might be needed:
# options(error = recover)
spe <-
    computeSumFactors(spe,
        clusters = spe$scran_quick_cluster,
        BPPARAM = MulticoreParam(4)
    )
Sys.time()

summary(sizeFactors(spe))

spe <- logNormCounts(spe)

## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
Sys.time()
dec <- modelGeneVar(spe,
    block = spe$sample_id,
    BPPARAM = MulticoreParam(4)
)
Sys.time()

pdf(
    here::here("plots", "scran_modelGeneVar.pdf"),
    useDingbats = FALSE
)
mapply(function(block, blockname) {
    plot(
        block$mean,
        block$total,
        xlab = "Mean log-expression",
        ylab = "Variance",
        main = blockname
    )
    curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
    )
}, dec$per.block, names(dec$per.block))
dev.off()



top.hvgs <- getTopHVGs(dec, prop = 0.1)
length(top.hvgs)


top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)


top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)


set.seed(20191112)
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs)
Sys.time()

reducedDimNames(spe)

summary(apply(reducedDim(spe, "PCA"), 2, sd))

summary(colMeans(reducedDim(spe, "PCA")))

Sys.time()
set.seed(20191206)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity50",
        perplexity = 50
    )
Sys.time()

Sys.time()
set.seed(20191206)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity5",
        perplexity = 5
    )
Sys.time()

Sys.time()
set.seed(20191206)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity20",
        perplexity = 20
    )
Sys.time()

Sys.time()
set.seed(20191206)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity80",
        perplexity = 80
    )
Sys.time()

Sys.time()
set.seed(20191206)
spe <- runUMAP(spe, dimred = "PCA", name = "UMAP_neighbors15")
Sys.time()

## Not yet implemented in spatialLIBD
# spatialLIBD::check_spe(spe)


# save
save(pce, file = here::here("rdata", "spe", "spe.rda"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
