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

## Related to https://github.com/LTLA/scuttle/issues/7#issuecomment-778710244
stopifnot(packageVersion("scuttle") >= "1.1.15")

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

## For visualizing this later with spatialLIBD
spe_raw$overlaps_tissue <- factor(ifelse(inTissue(spe_raw), "in", "out"))

## Create two versions: one with and one without filtering by tissue spot
spe_raw <- spe
pryr::object_size(spe_raw)
# 970 MB
dim(spe_raw)
# [1] 27013 59904

## Save the raw version now
dir.create(here::here("rdata"), showWarnings = FALSE)
dir.create(here::here("rdata", "spe"), showWarnings = FALSE)
save(spe_raw, file = here::here("rdata", "spe", "spe_raw.Rdata"))

## Work with SPE
spe <- spe_raw[, which(inTissue(spe_raw))]

## Remove spots without counts
spe <- spe[, -which(colSums(counts(spe)) == 0)]

dim(spe)
# [1] 27013 49999

pryr::object_size(spe)
# 939 MB

## Inspect in vs outside of tissue
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
# [1] "2021-02-17 10:23:01 EST"
# [1] "2021-02-17 10:25:47 EST"

Sys.time()
## Might be needed:
# options(error = recover)
spe <-
    computeSumFactors(spe,
        clusters = spe$scran_quick_cluster,
        BPPARAM = MulticoreParam(4)
    )
Sys.time()
# [1] "2021-02-17 10:28:34 EST"
# [1] "2021-02-17 10:35:13 EST"
## Related to https://github.com/LTLA/scuttle/issues/7
# In .rescale_clusters(clust.profile, ref.col = ref.clust, min.mean = min.mean) :
#   inter-cluster rescaling factor for cluster 64 is not strictly positive,
# reverting to the ratio of average library sizes

table(spe$scran_quick_cluster)
 #   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
 # 763  466  634  848  841  513  856  462  468  776  998  148  254  284  308  198
 #  17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32
 # 245  962  674  435  359  268  271  142  760  606  544  150  123  462  514  576
 #  33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
 # 176  290  161  183  697  662  579  199  190  235  197  424  111  186  432  259
 #  49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64
 # 684  818  271  884  586  178  102  161  920  893  380  877  675  537  104  542
 #  65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80
 # 935  461  327  207  229  366  166  110 1193  860 1219  579  679 1309 1476 1329
 #  81   82   83   84   85   86   87   88   89   90   91   92   93   94   95
 # 650  644  393  992  159  542  514  269  272  896  787  821 1141  213  760

summary(sizeFactors(spe))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000008 0.510786 0.862512 1.000000 1.330087 6.854382

spe <- logNormCounts(spe)
pryr::object_size(spe)

## Temporarely save this object (eventually it'll be done at the end only)
save(spe, file = here::here("rdata", "spe", "spe.rda"))
# 1.53 GB

## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
    block = spe$sample_id,
    BPPARAM = MulticoreParam(4)
)

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
# [1] 1769

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# [1] 15343

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
# [1] 14484

save(top.hvgs,
    top.hvgs.fdr5,
    top.hvgs.fdr1,
    file = here::here("rdata", "spe", "top.hvgs.rda"))

set.seed(20191112)
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs)
Sys.time()
# [1] "2021-02-17 10:41:57 EST"
# [1] "2021-02-17 10:43:06 EST"

reducedDimNames(spe)

summary(apply(reducedDim(spe, "PCA"), 2, sd))
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.8947  0.9127  0.9308  1.1226  0.9826  3.7028
summary(colMeans(reducedDim(spe, "PCA")))
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -2.466e-14 -3.442e-15  2.650e-15  2.456e-15  9.457e-15  2.697e-14

## From https://github.com/davismcc/scater/blob/master/R/runTSNE.R#L85
## I see that the default perplexity will be 50
# > mat <- scater:::.get_mat_from_sce(spe, exprs_values = 'logcounts', dimred = 'PCA', n_dimred = NULL)
# > dim(mat)
# [1] 49999    50
# > min(50, floor(nrow(mat) / 5))
# [1] 50

Sys.time()
set.seed(20191206)
spe <-
    runTSNE(spe,
        dimred = "PCA",
        name = "TSNE_perplexity50",
        perplexity = 50
    )
Sys.time()
# [1] "2021-02-17 10:45:30 EST"

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
save(spe, file = here::here("rdata", "spe", "spe.rda"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
