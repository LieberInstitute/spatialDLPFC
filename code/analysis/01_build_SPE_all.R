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
    "DLPFC_Br6423_post_extra_reads",
    "DLPFC_Br8492_ant_manual_alignment",
    "DLPFC_Br8492_mid_manual_alignment_extra_reads",
    "DLPFC_Br8492_post_manual_alignment",
    "Round2/DLPFC_Br2720_ant_manual_alignment",
    "Round2/DLPFC_Br2720_mid_manual_alignment",
    "Round2/DLPFC_Br2720_post_extra_reads",
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
    "Round3/DLPFC_Br8667_ant_extra_reads",
    "Round3/DLPFC_Br8667_mid_manual_alignment_all",
    "Round3/DLPFC_Br8667_post_manual_alignment_all",
    "Round4/DLPFC_Br2720_ant_2",
    "Round4/DLPFC_Br6432_ant_2",
    "Round4/DLPFC_Br8325_mid_2"
  ),
  subjects = c(rep(
    c("Br2743", "Br3942", "Br6423", "Br8492", "Br2720","Br6432","Br6471","Br6522","Br8325","Br8667"),
    each = 3
  ),"Br2720","Br6432","Br8325"),
  regions = c(rep(
    c("anterior", "middle", "posterior"),
    10
  ),"anterior","anterior","middle"),
  sex = c(rep(
    c("M", "M", "M", "F","F","M","M","M","F","F"),
    each = 3
  ),"F","M","F"),
  age = c(rep(
    c(61.54, 47.53, 51.73, 53.40,48.22,48.88,55.46,33.39,57.62,37.33),
    each = 3
  ),48.22,48.88,57.62),
  diagnosis = "Control"
)
sample_info$sample_path <- file.path(
  here::here("processed-data", "NextSeq"),
  sample_info$sample_id,
  "outs"
)
stopifnot(all(file.exists(sample_info$sample_path)))

#clean up sample_id
sample_info$sample_id <- gsub("_all|_extra_reads|DLPFC_|_manual_alignment", "", basename(sample_info$sample_id))

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
# [1] "2021-09-28 14:38:13 EDT"
# [1] "2021-09-28 14:42:47 EDT"

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

dim(spe)
# 36601 164733

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)

length(no_expr)
# [1] 7468
length(no_expr) / nrow(spe) * 100
# [1] 20.40381

## Create two versions: one with and one without filtering by tissue spot
spe_raw <- spe

pryr::object_size(spe_raw)
# 2,566,897,847 B
dim(spe_raw)
# [1] 36601 164733

## Save the raw version now
dir.create(here::here("rdata"), showWarnings = FALSE)
dir.create(here::here("rdata", "spe"), showWarnings = FALSE)
save(spe_raw, file = here::here("processed-data", "rdata", "spe", "spe_raw_final.Rdata"))

#remove genes that are not expressed in any spots
spe <- spe_raw[-no_expr, ]

dim(spe)
# 29133 164733

## Work with SPE
spe <- spe[, which(spatialData(spe_raw)$in_tissue=="TRUE")]
dim(spe)
#[1]  29133 129437

## Remove spots without counts
spe <- spe[, -which(colSums(counts(spe)) == 0)]
dim(spe)
# [1] 29133129430

save(spe, file = here::here("processed-data", "rdata", "spe", "spe_final.Rdata"))

pryr::object_size(spe)
# 2,473,727,880 B

## Inspect in vs outside of tissue
vis_grid_clus(
  spe = spe_raw,
  clustervar = "in_tissue",
  pdf = here::here("plots", "in_tissue_grid_final.pdf"),
  sort_clust = FALSE,
  colors = c("TRUE" = "grey90", "FALSE" = "orange")
)

summary(spe_raw$sum_umi[which(spatialData(spe_raw)$in_tissue=="FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0    66.0   170.5   287.0   318.0 39053.0

head(table(spe_raw$sum_umi[which(spatialData(spe_raw)$in_tissue=="FALSE")]))
# 0  1  2  3  4  5
# 12 18 50 59 67 87

vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "sum_umi",
  pdf = here::here("plots", "out_tissue_sum_umi_all_finall.pdf"),
  assayname = "counts"
)

summary(spe_raw$sum_gene[which(spatialData(spe_raw)$in_tissue=="FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0    55.0   133.0   203.4   242.0  7956.0

vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "sum_gene",
  pdf = here::here("plots", "out_tissue_sum_gene_all_final.pdf"),
  assayname = "counts"
)

summary(spe_raw$expr_chrM_ratio[which(spatialData(spe_raw)$in_tissue=="FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#  0.0000  0.1515  0.1963  0.2084  0.2518  1.0000      12
vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "out_tissue_expr_chrM_ratio_all_all.pdf"),
  assayname = "counts"
)

summary(spe$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1    1346    2315    2701    3579   46881
vis_grid_gene(
  spe = spe,
  geneid = "sum_umi",
  pdf = here::here("plots", "in_tissue_sum_umi_all_final.pdf"),
  assayname = "counts"
)

summary(spe$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1     908    1425    1516    2012    8344
vis_grid_gene(
  spe = spe,
  geneid = "sum_gene",
  pdf = here::here("plots", "in_tissue_sum_gene_all_final.pdf"),
  assayname = "counts"
)

summary(spe$expr_chrM_ratio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.07488 0.10776 0.12361 0.15866 1.00000
vis_grid_gene(
  spe = spe,
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "in_tissue_expr_chrM_ratio_all_final.pdf"),
  assayname = "counts"
)

vis_grid_gene(
  spe = spe_raw,
  geneid = "sum_umi",
  pdf = here::here("plots", "all_sum_umi_final.pdf"),
  assayname = "counts"
)

vis_grid_gene(
  spe = spe_raw,
  geneid = "sum_gene",
  pdf = here::here("plots", "all_sum_gene_final.pdf"),
  assayname = "counts"
)
vis_grid_gene(
  spe = spe_raw,
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "all_expr_chrM_ratio_final.pdf"),
  assayname = "counts"
)



## Quality control (scran)
qcstats <- perCellQCMetrics(spe, subsets = list(
  Mito = which(seqnames(spe) == "chrM")
))
qcfilter <- quickPerCellQC(qcstats, sub.fields="subsets_Mito_percent")
colSums(as.matrix(qcfilter))
# low_lib_size            low_n_features high_subsets_Mito_percent
# 5126                    6105                      3490
# discard
# 9035

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
save(spe, file = here::here("processed-data", "rdata", "spe", "spe_final.Rdata"))

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
# [1] "2021-10-28 13:05:32 EDT"
# [1]  "2021-10-28 13:14:52 EDT"

Sys.time()
## Might be needed:
# options(error = recover)
spe <-
  computeSumFactors(spe,
                    clusters = spe$scran_quick_cluster,
                    BPPARAM = MulticoreParam(4)
  )
Sys.time()
# [1] "2021-10-28 13:14:52 EDT"
# [1] "2021-10-28 13:16:36 EDT"

## Related to https://github.com/LTLA/scuttle/issues/7
# In .rescale_clusters(clust.profile, ref.col = ref.clust, min.mean = min.mean) :
#   inter-cluster rescaling factor for cluster 64 is not strictly positive,
# reverting to the ratio of average library sizes

table(spe$scran_quick_cluster)


summary(sizeFactors(spe))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.000034  0.471505  0.840405  1.000000  1.333201 16.712342

spe <- logNormCounts(spe)
pryr::object_size(spe)
# 4,047,023,760 B, 3.8gb

## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
                    block = spe$sample_id,
                    BPPARAM = MulticoreParam(4)
)

pdf(
  here::here("plots", "scran_modelGeneVar_final.pdf"),
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
# [1] 2059

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# [1] 18830

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
# [1] 17900

save(top.hvgs,
     top.hvgs.fdr5,
     top.hvgs.fdr1,
     file = here::here("processed-data", "rdata", "spe", "top.hvgs_all_final.Rdata"))


set.seed(20191112)
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = 50)
Sys.time()
#  "2021-09-29 15:24:44 EDT"
# "2021-09-29 15:28:19 EDT"

reducedDimNames(spe)
# [1] "PCA"
dim(reducedDim(spe, "PCA"))
## 129430   50

#make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(spe,"PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow
# 50

pdf(
  here::here("plots", "pca_elbow_final.pdf"),
  useDingbats = FALSE
)
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")
dev.off()



summary(apply(reducedDim(spe, "PCA"), 2, sd))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9191 0.9321  0.9488  1.1913  1.0570  4.2217

summary(colMeans(reducedDim(spe, "PCA")))

## From https://github.com/davismcc/scater/blob/master/R/runTSNE.R#L85
## I see that the default perplexity will be 50
# > mat <- scater:::.get_mat_from_sce(spe, exprs_values = 'logcounts', dimred = 'PCA', n_dimred = NULL)
# > dim(mat)
# [1] 49999    50
# > min(50, floor(nrow(mat) / 5))
# [1] 50

save(spe, file = here::here("processed-data", "rdata", "spe", "spe_final.Rdata"))
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
save(spe, file = here::here("processed-data","rdata", "spe", "spe_final.Rdata"))


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


#make downsampled version of spe object 
if(spatialData(spe)$array_col < 64 & spatialData(spe)$array_row < 39) {
  spe$quadrant<-"topleft"} else if(spatialData(spe)$array_col >= 64 & spatialData(spe)$array_row < 39){
  spe$quadrant<-"topright"}else if(spatialData(spe)$array_col < 64 & spatialData(spe)$array_row >= 39){
  spe$quadrant<-"bottomleft"}else if(spatialData(spe)$array_col >= 64 & spatialData(spe)$array_row >= 39){
  spe$quadrant<-"bottomright"}else{
  spe$quadrant<-"NA"}


spe$sample_quadrant <- paste0(spe$sampleid, "_", spe$quadrant)
sample_quad_list <- rafalib::splitit(spe$sample_quadrant)
selected_sample_quad <- unlist(lapply(sample_squad_list, sample, n = 250))
spe_subsampled <- spe[, selected_sample_quad]

# use getClusteredPCs as another way to determine the correct number of PCs to use
# get.PCs <- getClusteredPCs(reducedDim(spe))
# save(get.PCs, file = here::here("processed-data", "rdata", "spe", "get_PCs_092121.Rdata"))


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



