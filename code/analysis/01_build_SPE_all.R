# start screen before qrsh ex. screen -S nameIwant
# check cluster resources  qpic -q bluejay
#spatialDLPFC $ qrsh -l bluejay,mem_free=150G,h_vmem=150G
#~ $ cd /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
# R


## Automatically style the code in this script:
styler::style_file(
  here::here("analysis", "01_build_SPE_all.R"),
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
library("scran") ## requires uwot for UMAP
library("scater")
library("BiocParallel")
library("PCAtools")

## Related to https://github.com/LTLA/scuttle/issues/7#issuecomment-778710244
stopifnot(packageVersion("scuttle") >= "1.1.15")

## Define some info for the samples
sample_info <- data.frame(
  sample_id = c(
    "DLPFC_Br2743_ant_manual_alignment",
    "DLPFC_Br2743_mid_manual_alignment_extra_reads",
    "DLPFC_Br2743_post_manual_alignment",
    "DLPFC_Br3942_ant_manual_alignment",
    "DLPFC_Br3942_mid_manual_alignment",
    "DLPFC_Br3942_post_manual_alignment",
    "DLPFC_Br6423_ant_manual_alignment_extra_reads",
    "DLPFC_Br6423_mid_manual_alignment",
    "DLPFC_Br6423_post_manual_alignment",
    "DLPFC_Br8492_ant_manual_alignment",
    "DLPFC_Br8492_mid_manual_alignment_extra_reads",
    "DLPFC_Br8492_post_manual_alignment",
    "Round2/DLPFC_Br2720_ant_manual_alignment",
    "Round2/DLPFC_Br2720_mid_manual_alignment",
    "Round2/DLPFC_Br2720_post_manual_alignment",
    "Round2/DLPFC_Br6432_ant_manual_alignment",
    "Round2/DLPFC_Br6432_mid_manual_alignment",
    "Round2/DLPFC_Br6432_post_manual_alignment",
    "Round3/DLPFC_Br6471_ant_manual_alignment_all",
    "Round3/DLPFC_Br6471_mid_manual_alignment_all",
    "Round3/DLPFC_Br6471_post_manual_alignment_all",
    "Round3/DLPFC_Br6522_ant_manual_alignment_all",
    "Round3/DLPFC_Br6522_mid_manual_alignment_all",
    "Round3/DLPFC_Br6522_post_manual_alignment_all",
    "Round3/DLPFC_Br8325_ant_manual_alignment_all",
    "Round3/DLPFC_Br8325_mid_manual_alignment_all",
    "Round3/DLPFC_Br8325_post_manual_alignment_all",
    "Round3/DLPFC_Br8667_ant_manual_alignment_all",
    "Round3/DLPFC_Br8667_mid_manual_alignment_all",
    "Round3/DLPFC_Br8667_post_manual_alignment_all"
  ),
  subjects = rep(
    c("Br2743", "Br3942", "Br6423", "Br8492", "Br2720","Br6432","Br6471","Br6522","Br8325","Br8667"),
    each = 3
  ),
  regions = rep(
    c("anterior", "middle", "posterior"),
    10
  ),
  sex = rep(
    c("M", "M", "M", "F","F","M","M","M","F","F"),
    each = 3
  ),
  age = rep(
    c(61.54, 47.53, 51.73, 53.40,48.22,48.88,55.46,33.39,57.62,37.33),
    each = 3
  ),
  diagnosis = "Control"
)
sample_info$sample_path <- file.path(
  here::here("processed-data", "NextSeq"),
  sample_info$sample_id,
  "outs"
)
stopifnot(all(file.exists(sample_info$sample_path)))

## Re-shuffle to match default spatialLIBD grids
sample_info <- sample_info[rep(c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28), 3) + rep(0:2, each = 10), ]

## For testing the code with a subset of the data
# if(FALSE) {
#     sample_info <- sample_info[1:2, ]
# }

## Build SPE object
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
# [1] "2021-09-08 09:04:34 EDT"
# [1] "2021-09-08 09:08:37 EDT"

##at this point the sample info is messed up

#clean up sample_id
spe$sample_id <- gsub("Round2/", "", spe$sample_id)
spe$sample_id <- gsub("Round3/", "", spe$sample_id)
spe$sample_id <- gsub("_all", "", spe$sample_id)
spe$sample_id <- gsub("_extra_reads", "", spe$sample_id)

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
# [1] 7627
length(no_expr) / nrow(spe) * 100
# [1] 20.83823

#remove genes that are not expressed in any spots
spe <- spe[-no_expr, ]

## Create two versions: one with and one without filtering by tissue spot
spe_raw <- spe
pryr::object_size(spe_raw)
# 2,359,020,128 B
dim(spe_raw)
# [1] 28974 149757

## Save the raw version now
dir.create(here::here("rdata"), showWarnings = FALSE)
dir.create(here::here("rdata", "spe"), showWarnings = FALSE)
save(spe_raw, file = here::here("processed-data", "rdata", "spe", "spe_raw_090821.Rdata"))

## Work with SPE
spe <- spe_raw[, which(spatialData(spe_raw)$in_tissue=="TRUE")]
dim(spe)
#[1]  28974 118010

## Remove spots without counts
spe <- spe[, -which(colSums(counts(spe)) == 0)]

dim(spe)
# [1] 28974 118003

save(spe, file = here::here("processed-data", "rdata", "spe", "spe_090821.Rdata"))

#load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/spe_090821.Rdata")
#load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/spe_raw_090221.Rdata")

pryr::object_size(spe)
# 2,276,730,184 B

## Inspect in vs outside of tissue
vis_grid_clus(
  spe = spe_raw,
  clustervar = "in_tissue",
  pdf = here::here("plots", "in_tissue_grid.pdf"),
  sort_clust = FALSE,
  colors = c("TRUE" = "grey90", "FALSE" = "orange")
)

summary(spe_raw$sum_umi[which(spatialData(spe_raw)$in_tissue=="FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     0.0    84.0   181.0   292.8   320.0 39053.0

head(table(spe_raw$sum_umi[which(spatialData(spe_raw)$in_tissue=="FALSE")]))
# 0 1 2 3 4 5
# 9 15 35 40 37 54

vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "sum_umi",
  pdf = here::here("plots", "out_tissue_sum_umi_all.pdf"),
  assayname = "counts"
)

summary(spe_raw$sum_gene[which(spatialData(spe_raw)$in_tissue=="FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0    69.0   141.0   207.6   245.0  7956.0
vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "sum_gene",
  pdf = here::here("plots", "out_tissue_sum_gene_all.pdf"),
  assayname = "counts"
)

summary(spe_raw$expr_chrM_ratio[which(spatialData(spe_raw)$in_tissue=="FALSE")])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.0000  0.1542  0.1973  0.2101  0.2519  1.0000       9
vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "out_tissue_expr_chrM_ratio_all.pdf"),
  assayname = "counts"
)


summary(spe$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1    1409    2328    2725    3577   46881
vis_grid_gene(
  spe = spe,
  geneid = "sum_umi",
  pdf = here::here("plots", "in_tissue_sum_umi_all.pdf"),
  assayname = "counts"
)

summary(spe$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1     927    1433    1529    2014    8344
vis_grid_gene(
  spe = spe,
  geneid = "sum_gene",
  pdf = here::here("plots", "in_tissue_sum_gene_all.pdf"),
  assayname = "counts"
)

summary(spe$expr_chrM_ratio)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.00000 0.06035 0.08099 0.09029 0.10935 1.00000       7
vis_grid_gene(
  spe = spe,
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "in_tissue_expr_chrM_ratio_all.pdf"),
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
# low_lib_size            low_n_features high_subsets_Mito_percent
# 3869                      4578                      3745
# discard
# 7776
## Prior to dropping spots with 0 counts and checking for high chrM,
## this was the output:

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

save(spe, file = here::here("processed-data", "rdata", "spe", "spe_090821.Rdata"))

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
# [1] "2021-09-08 09:52:20 EDT"
# [1] "2021-09-08 09:56:31 EDT"

Sys.time()
## Might be needed:
# options(error = recover)
spe <-
  computeSumFactors(spe,
                    clusters = spe$scran_quick_cluster,
                    BPPARAM = MulticoreParam(4)
  )
Sys.time()
# [1] "2021-09-08 09:57:55 EDT"
# [1] "2021-09-02 13:15:54 EDT"
## Related to https://github.com/LTLA/scuttle/issues/7
# In .rescale_clusters(clust.profile, ref.col = ref.clust, min.mean = min.mean) :
#   inter-cluster rescaling factor for cluster 64 is not strictly positive,
# reverting to the ratio of average library sizes

table(spe$scran_quick_cluster)


summary(sizeFactors(spe))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.000034  0.478379  0.838123  1.000000  1.324870 16.630121

spe <- logNormCounts(spe)
pryr::object_size(spe)
# 3,723,808,272 B

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
# [1] 2039

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# [1] 18525

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
# [1] 17900

save(top.hvgs,
     top.hvgs.fdr5,
     top.hvgs.fdr1,
     file = here::here("processed-data", "rdata", "spe", "top.hvgs_all.Rdata"))


set.seed(20191112)
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = 50)
Sys.time()
# [1] "2021-09-21"
# [1] "2021-09-21"
load(file = here::here("processed-data", "rdata", "spe", "spe_090821.Rdata"))
load(file = here::here("processed-data", "rdata", "spe", "top.hvgs_all.Rdata"))

reducedDimNames(spe)
# [1] "PCA"
dim(reducedDim(spe, "PCA"))

#deafult is 50 PCs, check on OSTA for the number of PCs I should use. Make github issue for this. 

#make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(spe,"PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow

pdf(
  here::here("plots", "pca_elbow.pdf"),
  useDingbats = FALSE
)
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")
dev.off()

#I chose 8 as the elbow

# use getClusteredPCs as another way to determine the correct number of PCs to use
get.PCs <- getClusteredPCs(reducedDim(spe))
save(get.PCs, file = here::here("processed-data", "rdata", "spe", "get_PCs_092121.Rdata"))

#####redo later
summary(apply(reducedDim(spe, "PCA"), 2, sd))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9129  0.9292  0.9476  1.1828  1.0557  4.1074
summary(colMeans(reducedDim(spe, "PCA")))
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -3.300e-14 -9.055e-15 -2.285e-16  2.575e-17  7.583e-15  2.568e-14

## From https://github.com/davismcc/scater/blob/master/R/runTSNE.R#L85
## I see that the default perplexity will be 50
# > mat <- scater:::.get_mat_from_sce(spe, exprs_values = 'logcounts', dimred = 'PCA', n_dimred = NULL)
# > dim(mat)
# [1] 49999    50
# > min(50, floor(nrow(mat) / 5))
# [1] 50

save(spe, file = here::here("processed-data", "rdata", "spe", "spe_092121.Rdata"))
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
# [1] "2021-02-17 11:02:59 EST"

Sys.time()
set.seed(20191206)
spe <-
  runTSNE(spe,
          dimred = "PCA",
          name = "TSNE_perplexity5",
          perplexity = 5
  )
Sys.time()
# [1] "2021-09-08 11:14:41 EDT"
# [1] "2021-09-08 12:03:04 EDT"

Sys.time()
set.seed(20191206)
spe <-
  runTSNE(spe,
          dimred = "PCA",
          name = "TSNE_perplexity20",
          perplexity = 20
  )
Sys.time()
# [1] "2021-09-08 12:17:38 EDT"
# [1] "2021-09-08 13:09:17 EDT"

Sys.time()
set.seed(20191206)
spe <-
  runTSNE(spe,
          dimred = "PCA",
          name = "TSNE_perplexity80",
          perplexity = 80
  )
Sys.time()
# [1] "2021-09-08 13:17:46 EDT"
# [1] "2021-02-17 12:58:33 EST"

Sys.time()
set.seed(20191206)
spe <- runUMAP(spe, dimred = "PCA", name = "UMAP_neighbors15")
Sys.time()
# [1] "2021-09-08 17:36:28 EDT"
# [1] "2021-09-08 17:42:38 EDT"

## Not yet implemented in spatialLIBD
# spatialLIBD::check_spe(spe)


# save
save(spe, file = here::here("rdata", "spe", "spe_090821.Rdata"))


#######end redo later 

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R Under development (unstable) (2021-02-17 r80017)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-02-17
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version  date       lib source
#  AnnotationDbi            1.53.1   2021-02-04 [2] Bioconductor
#  AnnotationHub            2.23.2   2021-02-05 [2] Bioconductor
#  assertthat               0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  attempt                  0.3.1    2020-05-03 [1] CRAN (R 4.1.0)
#  beachmat                 2.7.6    2021-01-15 [2] Bioconductor
#  beeswarm                 0.2.3    2016-04-25 [1] CRAN (R 4.1.0)
#  benchmarkme              1.0.5    2021-02-09 [1] CRAN (R 4.1.0)
#  benchmarkmeData          1.0.4    2020-04-23 [1] CRAN (R 4.1.0)
#  Biobase                * 2.51.0   2020-10-27 [2] Bioconductor
#  BiocFileCache            1.15.1   2020-11-09 [2] Bioconductor
#  BiocGenerics           * 0.37.1   2021-02-04 [2] Bioconductor
#  BiocIO                   1.1.2    2020-12-05 [2] Bioconductor
#  BiocManager              1.30.10  2019-11-16 [2] CRAN (R 4.1.0)
#  BiocNeighbors            1.9.4    2020-12-17 [1] Bioconductor
#  BiocParallel           * 1.25.4   2021-02-04 [2] Bioconductor
#  BiocSingular             1.7.2    2021-01-23 [1] Bioconductor
#  BiocVersion              3.13.1   2020-10-27 [2] Bioconductor
#  Biostrings               2.59.2   2020-12-18 [2] Bioconductor
#  bit                      4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
#  bit64                    4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
#  bitops                   1.0-6    2013-08-17 [2] CRAN (R 4.1.0)
#  blob                     1.2.1    2020-01-20 [2] CRAN (R 4.1.0)
#  bluster                  1.1.5    2021-01-14 [1] Bioconductor
#  cachem                   1.0.4    2021-02-13 [2] CRAN (R 4.1.0)
#  cli                      2.3.0    2021-01-31 [2] CRAN (R 4.1.0)
#  codetools                0.2-18   2020-11-04 [3] CRAN (R 4.1.0)
#  colorout                 1.2-2    2021-02-11 [1] Github (jalvesaq/colorout@726d681)
#  colorspace               2.0-0    2020-11-11 [2] CRAN (R 4.1.0)
#  config                   0.3.1    2020-12-17 [1] CRAN (R 4.1.0)
#  cowplot                  1.1.1    2020-12-30 [1] CRAN (R 4.1.0)
#  crayon                   1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
#  curl                     4.3      2019-12-02 [2] CRAN (R 4.1.0)
#  data.table               1.13.6   2020-12-30 [2] CRAN (R 4.1.0)
#  DBI                      1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
#  dbplyr                   2.1.0    2021-02-03 [2] CRAN (R 4.1.0)
#  DelayedArray             0.17.7   2020-12-26 [2] Bioconductor
#  DelayedMatrixStats       1.13.5   2021-02-04 [2] Bioconductor
#  desc                     1.2.0    2018-05-01 [2] CRAN (R 4.1.0)
#  digest                   0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
#  dockerfiler              0.1.3    2019-03-19 [1] CRAN (R 4.1.0)
#  doParallel               1.0.16   2020-10-16 [2] CRAN (R 4.1.0)
#  dotCall64                1.0-1    2021-02-11 [1] CRAN (R 4.1.0)
#  dplyr                    1.0.4    2021-02-02 [2] CRAN (R 4.1.0)
#  dqrng                    0.2.1    2019-05-17 [1] CRAN (R 4.1.0)
#  DropletUtils             1.11.10  2021-02-04 [1] Bioconductor
#  DT                       0.17     2021-01-06 [2] CRAN (R 4.1.0)
#  edgeR                    3.33.1   2021-01-11 [2] Bioconductor
#  ellipsis                 0.3.1    2020-05-15 [2] CRAN (R 4.1.0)
#  ExperimentHub            1.17.1   2021-02-08 [2] Bioconductor
#  fastmap                  1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
#  fields                   11.6     2020-10-09 [2] CRAN (R 4.1.0)
#  filelock                 1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
#  foreach                  1.5.1    2020-10-15 [2] CRAN (R 4.1.0)
#  fs                       1.5.0    2020-07-31 [2] CRAN (R 4.1.0)
#  generics                 0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
#  GenomeInfoDb           * 1.27.6   2021-02-04 [2] Bioconductor
#  GenomeInfoDbData         1.2.4    2020-11-03 [2] Bioconductor
#  GenomicAlignments        1.27.2   2020-12-12 [2] Bioconductor
#  GenomicRanges          * 1.43.3   2021-01-14 [2] Bioconductor
#  ggbeeswarm               0.6.0    2017-08-07 [1] CRAN (R 4.1.0)
#  ggplot2                * 3.3.3    2020-12-30 [2] CRAN (R 4.1.0)
#  glue                     1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
#  golem                    0.2.1    2020-03-05 [1] CRAN (R 4.1.0)
#  gridExtra                2.3      2017-09-09 [2] CRAN (R 4.1.0)
#  gtable                   0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  HDF5Array                1.19.4   2021-02-14 [2] Bioconductor
#  here                   * 1.0.1    2020-12-13 [1] CRAN (R 4.1.0)
#  htmltools                0.5.1.1  2021-01-22 [2] CRAN (R 4.1.0)
#  htmlwidgets              1.5.3    2020-12-10 [2] CRAN (R 4.1.0)
#  httpuv                   1.5.5    2021-01-13 [2] CRAN (R 4.1.0)
#  httr                     1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
#  igraph                   1.2.6    2020-10-06 [2] CRAN (R 4.1.0)
#  interactiveDisplayBase   1.29.0   2020-10-27 [2] Bioconductor
#  IRanges                * 2.25.6   2020-12-18 [2] Bioconductor
#  irlba                    2.3.3    2019-02-05 [2] CRAN (R 4.1.0)
#  iterators                1.0.13   2020-10-15 [2] CRAN (R 4.1.0)
#  jsonlite                 1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
#  KEGGREST                 1.31.1   2020-11-23 [2] Bioconductor
#  knitr                    1.31     2021-01-27 [2] CRAN (R 4.1.0)
#  later                    1.1.0.1  2020-06-05 [2] CRAN (R 4.1.0)
#  lattice                  0.20-41  2020-04-02 [3] CRAN (R 4.1.0)
#  lazyeval                 0.2.2    2019-03-15 [2] CRAN (R 4.1.0)
#  lifecycle                1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
#  limma                    3.47.7   2021-02-15 [2] Bioconductor
#  locfit                   1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
#  magick                   2.6.0    2021-01-13 [2] CRAN (R 4.1.0)
#  magrittr                 2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
#  maps                     3.3.0    2018-04-03 [2] CRAN (R 4.1.0)
#  Matrix                 * 1.3-2    2021-01-06 [3] CRAN (R 4.1.0)
#  MatrixGenerics         * 1.3.1    2021-02-01 [2] Bioconductor
#  matrixStats            * 0.58.0   2021-01-29 [2] CRAN (R 4.1.0)
#  memoise                  2.0.0    2021-01-26 [2] CRAN (R 4.1.0)
#  metapod                  0.99.5   2020-12-14 [1] Bioconductor
#  mime                     0.10     2021-02-13 [2] CRAN (R 4.1.0)
#  munsell                  0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                   1.4.7    2020-11-20 [2] CRAN (R 4.1.0)
#  pkgconfig                2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  pkgload                  1.1.0    2020-05-29 [2] CRAN (R 4.1.0)
#  plotly                   4.9.3    2021-01-10 [2] CRAN (R 4.1.0)
#  png                      0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
#  Polychrome               1.2.6    2020-11-11 [1] CRAN (R 4.1.0)
#  promises                 1.2.0.1  2021-02-11 [1] CRAN (R 4.1.0)
#  pryr                   * 0.1.4    2018-02-18 [2] CRAN (R 4.1.0)
#  purrr                    0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R.methodsS3              1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
#  R.oo                     1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
#  R.utils                  2.10.1   2020-08-26 [2] CRAN (R 4.1.0)
#  R6                       2.5.0    2020-10-28 [2] CRAN (R 4.1.0)
#  rappdirs                 0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
#  RColorBrewer             1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
#  Rcpp                     1.0.6    2021-01-15 [2] CRAN (R 4.1.0)
#  RCurl                    1.98-1.2 2020-04-18 [2] CRAN (R 4.1.0)
#  remotes                  2.2.0    2020-07-21 [2] CRAN (R 4.1.0)
#  restfulr                 0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
#  rhdf5                    2.35.0   2020-10-27 [2] Bioconductor
#  rhdf5filters             1.3.3    2020-12-07 [2] Bioconductor
#  Rhdf5lib                 1.13.0   2020-10-27 [2] Bioconductor
#  rjson                    0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
#  rlang                    0.4.10   2020-12-30 [2] CRAN (R 4.1.0)
#  roxygen2                 7.1.1    2020-06-27 [2] CRAN (R 4.1.0)
#  rprojroot                2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  Rsamtools                2.7.1    2021-01-18 [2] Bioconductor
#  RSQLite                  2.2.3    2021-01-24 [2] CRAN (R 4.1.0)
#  rstudioapi               0.13     2020-11-12 [2] CRAN (R 4.1.0)
#  rsvd                     1.0.3    2020-02-17 [1] CRAN (R 4.1.0)
#  rtracklayer            * 1.51.4   2021-01-14 [2] Bioconductor
#  S4Vectors              * 0.29.7   2021-02-04 [2] Bioconductor
#  ScaledMatrix             0.99.2   2021-01-14 [1] Bioconductor
#  scales                   1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  scater                 * 1.19.9   2021-02-01 [1] Bioconductor
#  scatterplot3d            0.3-41   2018-03-14 [1] CRAN (R 4.1.0)
#  scran                  * 1.19.13  2021-02-13 [1] Bioconductor
#  scuttle                * 1.1.15   2021-02-14 [1] Bioconductor
#  sessioninfo            * 1.1.1    2018-11-05 [2] CRAN (R 4.1.0)
#  shiny                    1.6.0    2021-01-25 [2] CRAN (R 4.1.0)
#  shinyWidgets             0.5.7    2021-02-03 [1] CRAN (R 4.1.0)
#  SingleCellExperiment   * 1.13.10  2021-02-10 [1] Bioconductor
#  spam                     2.6-0    2020-12-14 [2] CRAN (R 4.1.0)
#  sparseMatrixStats        1.3.6    2021-02-04 [2] Bioconductor
#  SpatialExperiment      * 1.1.432  2021-02-17 [1] Github (drighelli/SpatialExperiment@8407ee8)
#  spatialLIBD            * 1.3.5    2021-02-17 [1] Github (LieberInstitute/spatialLIBD@3fbc875)
#  statmod                  1.4.35   2020-10-19 [2] CRAN (R 4.1.0)
#  stringi                  1.5.3    2020-09-09 [2] CRAN (R 4.1.0)
#  stringr                  1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
#  SummarizedExperiment   * 1.21.1   2020-12-12 [2] Bioconductor
#  testthat                 3.0.2    2021-02-14 [2] CRAN (R 4.1.0)
#  tibble                   3.0.6    2021-01-29 [2] CRAN (R 4.1.0)
#  tidyr                    1.1.2    2020-08-27 [2] CRAN (R 4.1.0)
#  tidyselect               1.1.0    2020-05-11 [2] CRAN (R 4.1.0)
#  usethis                  2.0.1    2021-02-10 [2] CRAN (R 4.1.0)
#  uwot                   * 0.1.10   2020-12-15 [1] CRAN (R 4.1.0)
#  vctrs                    0.3.6    2020-12-17 [2] CRAN (R 4.1.0)
#  vipor                    0.4.5    2017-03-22 [1] CRAN (R 4.1.0)
#  viridis                  0.5.1    2018-03-29 [2] CRAN (R 4.1.0)
#  viridisLite              0.3.0    2018-02-01 [2] CRAN (R 4.1.0)
#  withr                    2.4.1    2021-01-26 [2] CRAN (R 4.1.0)
#  xfun                     0.21     2021-02-10 [2] CRAN (R 4.1.0)
#  XML                      3.99-0.5 2020-07-23 [2] CRAN (R 4.1.0)
#  xml2                     1.3.2    2020-04-23 [2] CRAN (R 4.1.0)
#  xtable                   1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
#  XVector                  0.31.1   2020-12-12 [2] Bioconductor
#  yaml                     2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
#  zlibbioc                 1.37.0   2020-10-27 [2] Bioconductor
#
# [1] /users/lcollado/R/devel
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/library

Sys.time()
g_k50 <- buildSNNGraph(spe, k = 50, use.dimred = 'PCA')
Sys.time()
## About 12 minutes
# [1] "2019-11-13 15:20:32 EST"
# [1] "2019-11-13 15:31:50 EST"
save(g_k50, file=here::here("processed-data", "rdata","spe","g_k50_090921.Rdata"))

Sys.time()
g_walk_k50 <- igraph::cluster_walktrap(g_k50)
Sys.time()
save(g_walk_k50, file = here::here("processed-data", "rdata","spe","g_walk_k50_091421.Rdata"))
## For me, failed at 6 days

##try fewwer pcs (8) and lower K (20) to decrease runtime
Sys.time()
g_k20 <- buildSNNGraph(spe, k = 20, use.dimred = 'PCA')
Sys.time()
## 092121 About 11 minutes 
save(g_k20, file=here::here("processed-data", "rdata","spe","g_k20_092121.Rdata"))

Sys.time()
g_walk_k20 <- igraph::cluster_walktrap(g_k20)
Sys.time()
save(g_walk_k20, file = here::here("processed-data", "rdata","spe","g_walk_k20_092121.Rdata"))

clust_k20 <- sort_clusters(g_walk_k20$membership)
save(clust_k20, file = here::here("processed-data", "rdata","spe","clust_k20_092121.Rdata"))
#finished in 4.5 hours

### For the SNN graph with K = 50, find which nested subset best matches
## the clusters from 10x Genomics labeled by Kristen Maynard and Keri Martinowich
clust_k5_list <- lapply(4:28, function(n) {
  message(paste(Sys.time(), 'n =', n))
  sort_clusters(igraph::cut_at(g_walk_k20, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)
save(clust_k5_list, file = here::here("processed-data", "rdata", "clust_k5_092121_list.Rdata"))

## Add clusters to spe colData
col.names <- paste0("SNN_k10_k",4:28)
for (i in seq_along(col.names)){
  colData(spe) <- cbind(colData(spe),clust_k5_list[i])
}
colnames(colData(spe))[19:43] <- col.names

save(spe, file = here::here("processed-data", "rdata","spe", "spe_snn_clusters_092121.Rdata"))


#Now do 8 PCs and k = 50
Sys.time()
g_k50 <- buildSNNGraph(spe, k = 50, use.dimred = 'PCA')
Sys.time()
## 092121 About 11 minutes 
save(g_k50, file=here::here("processed-data", "rdata","spe","g_k50_092121.Rdata"))

Sys.time()
g_walk_k50 <- igraph::cluster_walktrap(g_k50)
Sys.time()
save(g_walk_k50, file = here::here("processed-data", "rdata","spe","g_walk_k50_092121.Rdata"))

clust_k50 <- sort_clusters(g_walk_k50$membership)
save(clust_k50, file = here::here("processed-data", "rdata","spe","clust_k50_092121.Rdata"))

### For the SNN graph with K = 50, find which nested subset best matches
## the clusters from 10x Genomics labeled by Kristen Maynard and Keri Martinowich
clust_k5_list <- lapply(4:28, function(n) {
  message(paste(Sys.time(), 'n =', n))
  sort_clusters(igraph::cut_at(g_walk_k50, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)
save(clust_k5_list, file = here::here("processed-data", "rdata", "clust_k5_092121_list_k50.Rdata"))

## Add clusters to spe colData
col.names <- paste0("SNN_k10_k",4:28)
for (i in seq_along(col.names)){
  colData(spe) <- cbind(colData(spe),clust_k5_list[i])
}
colnames(colData(spe))[19:43] <- col.names

save(spe, file = here::here("processed-data", "rdata","spe", "spe_snn_clusters_092121_k50.Rdata"))

#now do 50 pcs and k = 10
load(file=here::here("processed-data","rdata", "spe", "spe_090821.Rdata"))

Sys.time()
g_k20 <- buildSNNGraph(spe, k = 10, use.dimred = 'PCA')
Sys.time()
save(g_k20, file=here::here("processed-data", "rdata","spe","g_pc50_k10_092221.Rdata"))

Sys.time()
g_walk_k20 <- igraph::cluster_walktrap(g_k20)
Sys.time()
save(g_walk_k20, file = here::here("processed-data", "rdata","spe","g_walk_pc50_k10_092221.Rdata"))

clust_k20 <- sort_clusters(g_walk_k20$membership)
save(clust_k20, file = here::here("processed-data", "rdata","spe","clust_pc50_k10_092221.Rdata"))

### For the SNN graph with K = 50, find which nested subset best matches
## the clusters from 10x Genomics labeled by Kristen Maynard and Keri Martinowich
clust_k5_list <- lapply(4:28, function(n) {
  message(paste(Sys.time(), 'n =', n))
  sort_clusters(igraph::cut_at(g_walk_k20, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)
save(clust_k5_list, file = here::here("processed-data", "rdata", "clust_k5_pc50_k10_092221_list.Rdata"))

## Add clusters to spe colData
col.names <- paste0("SNN_k10_k",4:28)
for (i in seq_along(col.names)){
  colData(spe) <- cbind(colData(spe),clust_k5_list[i])
}
colnames(colData(spe))[19:43] <- col.names

save(spe, file = here::here("processed-data", "rdata","spe", "spe_snn_clusters_pc50_k10_092221.Rdata"))



