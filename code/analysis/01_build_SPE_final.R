# start screen before qrsh ex. screen -S nameIwant
# check cluster resources  qpic -q bluejay
#spatialDLPFC $ qrsh -l bluejay,mem_free=150G,h_vmem=150G
#~ $ cd /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
# R

## This script requires R 4.1
# module load conda_R/devel

## utils
library("here")
library("sessioninfo")
library("pryr")

## reading the data
library("SpatialExperiment")
library("rtracklayer")
library(readr)
library(readxl)
library(dplyr)
library(tidyr)

## vis
library("spatialLIBD")
library("RColorBrewer")
library(nlme)
library(ggplot2)

## analysis
library("scran") ## requires uwot for UMAP
library("scater")
library("BiocParallel")
library("PCAtools")
library(harmony)
library(uwot)
library(mclust)
library(aricode)
library(BayesSpace)
library(Polychrome)
library(patchwork)
library(broom)
library(magick)



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
    "DLPFC_Br2720_ant_2",
    "DLPFC_Br2720_mid_manual_alignment",
    "DLPFC_Br2720_post_extra_reads",
    "DLPFC_Br6432_ant_2",
    "DLPFC_Br6432_mid_manual_alignment",
    "DLPFC_Br6432_post_manual_alignment",
    "DLPFC_Br6471_ant_manual_alignment_all",
    "DLPFC_Br6471_mid_manual_alignment_all",
    "DLPFC_Br6471_post_manual_alignment_all",
    "DLPFC_Br6522_ant_manual_alignment_all",
    "DLPFC_Br6522_mid_manual_alignment_all",
    "DLPFC_Br6522_post_manual_alignment_all",
    "DLPFC_Br8325_ant_manual_alignment_all",
    "DLPFC_Br8325_mid_2",
    "DLPFC_Br8325_post_manual_alignment_all",
    "DLPFC_Br8667_ant_extra_reads",
    "DLPFC_Br8667_mid_manual_alignment_all",
    "DLPFC_Br8667_post_manual_alignment_all"
    
  ),
  subjects = c(rep(
    c("Br2743", "Br3942", "Br6423", "Br8492", "Br2720","Br6432","Br6471","Br6522","Br8325","Br8667"),
    each = 3
  )),
  regions = c(rep(
    c("anterior", "middle", "posterior"),
    10
  )),
  sex = c(rep(
    c("M", "M", "M", "F","F","M","M","M","F","F"),
    each = 3
  )),
  age = c(rep(
    c(61.54, 47.53, 51.73, 53.40,48.22,48.88,55.46,33.39,57.62,37.33),
    each = 3
  )),
  diagnosis = "Control"
)
sample_info$sample_path <- file.path(
  here::here("processed-data", "rerun_spaceranger"),
  sample_info$sample_id,
  "outs"
)
stopifnot(all(file.exists(sample_info$sample_path)))

#clean up sample_id
sample_info$sample_id <- gsub("_all|_extra_reads|DLPFC_|_manual_alignment", "", basename(sample_info$sample_id))

## Build SPE object
Sys.time()
spe <- read10xVisiumWrapper(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres", "hires", "detected", "aligned"),
  load = TRUE,
  reference_gtf = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
)
Sys.time()

#[1] "2021-12-15 10:59:02 EST"
#[1] "2021-12-15 11:10:57 EST"

## Add the experimental information
### Equivalent to:
## colData(spe)$subject <- ...
spe$subject <- sample_info$subjects[match(spe$sample_id, sample_info$sample_id)]
spe$region <- sample_info$regions[match(spe$sample_id, sample_info$sample_id)]
spe$sex <- sample_info$sex[match(spe$sample_id, sample_info$sample_id)]
spe$age <- sample_info$age[match(spe$sample_id, sample_info$sample_id)]
spe$diagnosis <- sample_info$diagnosis[match(spe$sample_id, sample_info$sample_id)]


## Add information used by spatialLIBD
is_mito <- which(seqnames(spe) == "chrM")
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi

## Add number of cells per spot
spe$cell_count <- NA
# ## Read in cell counts and segmentation results
# segmentations_list <- lapply(sample_info$sample_id, function(sampleid) {
#   current<-sample_info$sample_path[sample_info$sample_id==sampleid]
#   file <- file.path(current, "spatial", "tissue_spot_counts.csv")
#   if(!file.exists(file)) return(NULL)
#   x <- read.csv(file)
#   x$key <- paste0(x$barcode, "_", sampleid)
#   return(x)
# })
# ## Merge them (once the these files are done, this could be replaced by an rbind)
# segmentations <- Reduce(function(...) merge(..., all = TRUE), segmentations_list[lengths(segmentations_list) > 0])
# 
# ## Add the information
# segmentation_match <- match(spe_wrapper$key, segmentations$key)
# segmentation_info <- segmentations[segmentation_match, - which(colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")), drop=FALSE]
# colData(spe_wrapper) <- cbind(colData(spe_wrapper), segmentation_info)

dim(spe)
# [1]  36601 149757

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)

length(no_expr)
# [1] 7685
length(no_expr) / nrow(spe) * 100
# [1] 20.99669

## Create two versions: one with and one without filtering by tissue spot
spe_raw <- spe

pryr::object_size(spe_raw)
# 4,754,346,424 B
dim(spe_raw)
# [1] 36601 149757

## Save the raw version now
dir.create(here::here("rdata"), showWarnings = FALSE)
dir.create(here::here("rdata", "spe"), showWarnings = FALSE)
save(spe_raw, file = here::here("processed-data", "rdata", "spe", "spe_raw_final.Rdata"))

#remove genes that are not expressed in any spots
spe <- spe_raw[-no_expr, ]

dim(spe)
# 28916 149757

## Work with SPE
spe <- spe[, which(spatialData(spe_raw)$in_tissue=="TRUE")]
dim(spe)
#[1]  28916 118800

## Remove spots without counts
spe <- spe[, -which(colSums(counts(spe)) == 0)]
dim(spe)
# [1] 28916 118793

pryr::object_size(spe)
# 4,664,781,528 B

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
# 0.0    66.0   170.0   291.6   326.0 39053.0 

head(table(spe_raw$sum_umi[which(spatialData(spe_raw)$in_tissue=="FALSE")]))
# 0  1  2  3  4  5 
# 12 17 50 59 65 81 

vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "sum_umi",
  pdf = here::here("plots", "out_tissue_sum_umi_all_finall.pdf"),
  assayname = "counts"
)

summary(spe_raw$sum_gene[which(spatialData(spe_raw)$in_tissue=="FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    55.0   133.0   206.6   249.0  7956.0 

vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "sum_gene",
  pdf = here::here("plots", "out_tissue_sum_gene_all_final.pdf"),
  assayname = "counts"
)

summary(spe_raw$expr_chrM_ratio[which(spatialData(spe_raw)$in_tissue=="FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.1500  0.1947  0.2072  0.2500  1.0000      12 

vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "out_tissue_expr_chrM_ratio_all_all.pdf"),
  assayname = "counts"
)

summary(spe$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    1358    2310    2674    3540   46881 
vis_grid_gene(
  spe = spe,
  geneid = "sum_umi",
  pdf = here::here("plots", "in_tissue_sum_umi_all_final.pdf"),
  assayname = "counts"
)

summary(spe$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1     910    1426    1508    1999    8344 

vis_grid_gene(
  spe = spe,
  geneid = "sum_gene",
  pdf = here::here("plots", "in_tissue_sum_gene_all_final.pdf"),
  assayname = "counts"
)

summary(spe$expr_chrM_ratio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.07279 0.10440 0.12106 0.15513 1.00000 
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
# 5043                      6036                      3740 
# discard 
# 9183 


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
# [1] "2021-12-15 13:06:04 EST"
# [1] "2021-12-15 13:17:12 EST"

Sys.time()
#[1] "2021-12-15 13:19:48 EST"
## Might be needed:
# options(error = recover)
spe <-
  computeSumFactors(spe,
                    clusters = spe$scran_quick_cluster,
                    BPPARAM = MulticoreParam(4)
  )
Sys.time()
#"2021-12-07 13:24:55 EST"
# 1] "2021-12-15 13:40:18 EST"


## Related to https://github.com/LTLA/scuttle/issues/7
# In .rescale_clusters(clust.profile, ref.col = ref.clust, min.mean = min.mean) :
#   inter-cluster rescaling factor for cluster 64 is not strictly positive,
# reverting to the ratio of average library sizes

table(spe$scran_quick_cluster)


summary(sizeFactors(spe))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000034  0.477444  0.847973  1.000000  1.331274 16.813313 

spe <- logNormCounts(spe)
pryr::object_size(spe)
#6,101,374,224 B



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
# [1] 2016

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# [1] 18417

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
# [1] 17830

save(top.hvgs,
     top.hvgs.fdr5,
     top.hvgs.fdr1,
     file = here::here("processed-data", "rdata", "spe", "top.hvgs_all_final.Rdata"))

#Run PCA with BayesSpace 
set.seed(20211207)
Sys.time()
spe = spatialPreprocess(spe, n.PCs = 50)
Sys.time()
#2021-12-15 13:51:04 EST"
#"2021-12-15 14:54:10 EST"


dim(reducedDim(spe, "PCA"))
# 118793     50

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
# 0.9184  0.9343  0.9509  1.1872  1.0561  4.1944 

# RunUMAP
Sys.time()
spe = runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) = c("UMAP1", "UMAP2")
Sys.time()

#[1] "2021-12-15 15:22:25 EST"
#[1] "2021-12-15 15:28:23 EST"

#Run TSNE
Sys.time()
set.seed(20211207)
spe <-
  runTSNE(spe,
          dimred = "PCA",
          name = "TSNE_perplexity80",
          perplexity = 80
  )
Sys.time()
#[1] "2021-12-15 15:30:48 EST"
#[1] "2021-12-15 16:49:28 EST"

Sys.time()
save(spe, file = here::here("processed-data","rdata", "spe", "spe_final.Rdata"))
Sys.time()
#[1] "2021-12-15 16:59:47 EST"
# [1] "2021-12-15 17:10:13 EST"


#make plots ofUMAP
pdf(file=here::here("plots", "UMAP_subject.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe$subject))) +
  geom_point() +
  labs(color = "Subject") +
  theme_bw()
dev.off()

pdf(file=here::here("plots", "UMAP_sample_id.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()
dev.off()


###harmony batch correction
spe = RunHarmony(spe, "sample_id", verbose = F)
spe = runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")


pdf(file=here::here("plots", "UMAP_harmony_sample_id.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")),
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()
dev.off()

Sys.time()
save(spe, file = here::here("processed-data","rdata", "spe", "spe_final.Rdata"))
Sys.time()
#[1] "2021-12-16 15:02:28 EST"


#graph-based clustering, no batch correction
Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = 'PCA')
Sys.time()
# [1] "2021-12-16 01:09:06 EST"
#  "2021-12-16 01:24:11 EST"

save(g_k10, file=here::here("processed-data", "rdata","spe","g_k10.Rdata"))

Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
# "2021-12-16 01:27:49 EST"
# [1] "2021-12-16 09:26:14 EST"

save(g_walk_k10, file = here::here("processed-data", "rdata","spe","g_walk_k10.Rdata"))

clust_k10 <- sort_clusters(g_walk_k10$membership)
save(clust_k10, file = here::here("processed-data", "rdata","spe","clust_k10.Rdata"))
#use cluster_export

clust_k5_list <- lapply(4:28, function(n) {
  message(paste(Sys.time(), 'n =', n))
  sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)
save(clust_k5_list, file = here::here("processed-data", "rdata", "spe", "clust_k5_list.Rdata"))

## Add clusters to spe colData
cluster_colNames <- paste0("SNN_k10_k",4:28)
for (i in seq_along(cluster_colNames)){
  colData(spe) <- cbind(colData(spe),clust_k5_list[i])
}
colnames(colData(spe))[34:58] <- cluster_colNames

##make plot
nb.cols <- 25
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
sample_ids <- unique(colData(spe)$sample_id)
pdf(file = here::here("plots","vis_clus_graph_based_pca.pdf"))
for (i in seq_along(sample_ids)){
  for(j in seq_along(cluster_colNames)){
    my_plot <- vis_clus(
      spe = spe,
      clustervar = cluster_colNames[j],
      sampleid = sample_ids[i],
      colors =  mycolors,
      ... = paste0(" ",cluster_colNames[j])
    )
    print(my_plot)
  }
  
}
dev.off()


##graph-based on batch correct
Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = 'HARMONY')
Sys.time()
#"2021-12-16 16:17:31 EST"
#"2021-12-16 16:34:28 EST"
save(g_k10, file=here::here("processed-data", "rdata","spe","g_k10_harmony.Rdata"))

Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
#[1] "2021-12-17 13:05:59 EST"

save(g_walk_k10, file = here::here("processed-data", "rdata","spe","g_walk_k10_harmony.Rdata"))

clust_k10 <- sort_clusters(g_walk_k10$membership)
save(clust_k10, file = here::here("processed-data", "rdata","spe","clust_k10_harmony.Rdata"))

clust_k5_list <- lapply(4:28, function(n) {
  message(paste(Sys.time(), 'n =', n))
  sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)
save(clust_k5_list, file = here::here("processed-data", "rdata", "spe", "clust_k5_list_harmony.Rdata"))

## Add clusters to spe colData
cluster_colNames <- paste0("SNN_k10_k",4:28)
for (i in seq_along(cluster_colNames)){
  colData(spe) <- cbind(colData(spe),clust_k5_list[i])
}
colnames(colData(spe))[34:58] <- cluster_colNames

##make plot
sample_ids <- unique(colData(spe)$sample_id)
pdf(file = here::here("plots","vis_clus_graph_based_har.pdf"))
for (i in seq_along(sample_ids)){
  for(j in seq_along(cluster_colNames)){
    my_plot <- vis_clus(
      spe = spe,
      clustervar = cluster_colNames[j],
      sampleid = sample_ids[i],
      colors =  mycolors,
      ... = paste0(" ",cluster_colNames[j])
    )
    print(my_plot)
  }
  
}
dev.off()


##do offset so we can run BayesSpace
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <-unique(spe$sample_id)
spe$row <- spatialData(spe)$array_row + auto_offset_row[spe$sample_id]
spe$col <- spatialData(spe)$array_col

pdf(file=here::here("plots", "bayesSpace_offset_check.pdf"))
clusterPlot(spe, "subject", color = NA) + #make sure no overlap between samples
  labs(fill = "Subject", title = "Offset check")
dev.off()

#non-batch corrected bayesSpace terminal 2
spe = spatialCluster(spe, use.dimred = "PCA", q = 7, nrep = 10000)

pdf(file=here::here("plots", "bayesSpace_clusterPlot_PCA.pdf"))
    clusterPlot(spe, color = NA) + #plot clusters
      labs(title = "BayesSpace joint clustering")
  dev.off()

  cluster_export(
    spe,
    "spatial.cluster",
    cluster_dir = here::here("processed-data", "rdata", "spe", "bayesSpace_PCs.csv" )
  )


##bayesSpace on batch corrected data terminal 3
spe = spatialCluster(spe, use.dimred = "HARMONY", q = 7, nrep = 10000) 

pdf(file=here::here("plots", "bayesSpace_clusterPlot_harmony.pdf"))
clusterPlot(spe, color = NA) + #plot clusters
  labs(title = "BayesSpace joint clustering")
dev.off()

cluster_export(
  spe,
  "spatial.cluster",
  cluster_dir = here::here("processed-data", "rdata", "spe", "bayesSpace_harmony.csv" )
)

spe.enhanced <- spatialEnhance(spe, use.dimred = "HARMONY", q = 7, nrep = 10000,  burn.in=100)

pdf(file=here::here("plots", "bayesSpace_clusterPlot_harmony.pdf"))
clusterPlot(spe.enhaced, color = NA) + #plot clusters
  labs(title = "BayesSpace enhanced clustering")
dev.off()

cluster_export(
  spe.enhanced,
  "spatial.cluster",
  cluster_dir = here::here("processed-data", "rdata", "spe", "bayesSpace_enhanced_harmony" )
)

#semi_supervised 
#adapted from Luka's script https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/SpatialDE_clustering.Rmd

# sample names
sample_names <- paste0("sample_", unique(colData(spe)$sample_id))
sample_names

## Load pseudobulk genes (from Leo's analyses)

# load spreadsheet of significant genes for pseudobulk layers (from Leo's analyses)
sig_genes <- read_csv("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/sig_genes.csv")
sig_genes

genes_pseudobulk <- sig_genes[, c("ensembl", "gene")]
colnames(genes_pseudobulk) <- c("gene_id", "gene_name")
dim(genes_pseudobulk)
#[1] 490   2

#remove duplicates (i.e. genes identified for multiple layers)
genes_pseudobulk <- distinct(genes_pseudobulk)
dim(genes_pseudobulk)
# [1] 198   2

##Clustering

# -clustering on top 50 PCs on pseudobulk layer genes (from Leo's analyses; 198 genes)
# -clustering on top 10 UMAPs on pseudobulk layer genes (from Leo's analyses; 198 genes)

# parameters
n_umap <- 10
max_spatial <- 1
n_neighbors <- 10
n_clus <- 7 #changed from 8

d_plot <- data.frame()

#filter out pseudobulk genes that aren't in data, removes 1 gene
genes_pseudobulk <- genes_pseudobulk[which(genes_pseudobulk$gene_id %in% rownames(spe)),]

# ---------------------------------------------
# extract and calculate features (PCA and UMAP)
# ---------------------------------------------

### pseudobulk layer genes (from Leo's analyses; 198 genes)

# run PCA on pseudobulk layer genes (from Leo's analyses; 198 genes)
logcounts_pseudobulk <- logcounts(spe[genes_pseudobulk$gene_id, ])

# note: use 'prcomp' instead of 'calculatePCA' due to small number of genes
out_pca_pseudobulk <- prcomp(t(as.matrix(logcounts_pseudobulk)))$x[, 1:50]

dims_pseudobulk_PCA <- out_pca_pseudobulk
rownames(dims_pseudobulk_PCA) <- colnames(spe)
dim(dims_pseudobulk_PCA)
#[1] 298  50
stopifnot(nrow(dims_pseudobulk_PCA) == ncol(spe))

#add pseudobulk PCA to spe object 
reducedDims(spe)$pseudobulk_PCA <- dims_pseudobulk_PCA


# run UMAP on pseudobulk layer genes (from Leo's analyses; 197 genes)
set.seed(1234)
out_umap_pseudobulk <- umap(dims_pseudobulk_PCA, scale = TRUE, n_components = n_umap)

dims_pseudobulk_UMAP <- out_umap_pseudobulk
colnames(dims_pseudobulk_UMAP) <- paste0("UMAP", seq_len(n_umap))
rownames(dims_pseudobulk_UMAP) <- colnames(spe)
dim(dims_pseudobulk_UMAP)
stopifnot(nrow(dims_pseudobulk_UMAP) == ncol(spe))

reducedDims(spe)$pseudobulk_UMAP <- dims_pseudobulk_UMAP

#plot UMAP before batch correction
# ggplot(data.frame(reducedDim(sce.combined, "UMAP")), 
#        aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$sample_name))) +
#   geom_point() +
#   labs(color = "Sample") +
#   theme_bw()



###harmony batch correction
spe = RunHarmony(spe, "sample_id",reduction = "pseudobulk_PCA", verbose = F)

spe = runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")


pdf(file=here::here("plots", "UMAP_pseudobulk_harmony_sample_id.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")),
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()
dev.off()

Sys.time()
save(spe, file = here::here("processed-data","rdata", "spe", "spe_final_pseudobulk.Rdata"))
Sys.time()


# --------------------------------------------------------------------------------------------
# run clustering and calculate Adjusted Rand Index (ARI) / Normalized Mutual Information (NMI)
# --------------------------------------------------------------------------------------------

# using graph-based clustering (see Bioconductor OSCA book)

# convenience function; note uses some external variables from above
run_clustering <- function(input, method) {
  dims_clus <- input
  
  set.seed(1234)
  g <- buildSNNGraph(t(dims_clus), k = n_neighbors, d = ncol(dims_clus))
  g_walk <- igraph::cluster_walktrap(g)
  clus <- igraph::cut_at(g_walk, n = n_clus)
  clus <- sort_clusters(clus)
  
  table(clus)
  stopifnot(length(clus) == nrow(dims_clus))
  
  data.frame(
    spot_name = as.character(rownames(dims_clus)), 
    sample_name = as.character(colData(spe)$sample_id), 
    method = as.character(method), 
    cluster = as.numeric(clus),  
    stringsAsFactors = FALSE
  )
}

d_plot <- rbind(d_plot, run_clustering(reducedDim(spe, "HARMONY"), method = "pseudobulk_PCA"))
d_plot <- rbind(d_plot, run_clustering(reducedDim(spe, "UMAP.HARMONY"), method = "pseudobulk_UMAP"))


library(tidyr)
#divide d_plot by method
d_plot_wide <-as.data.frame(pivot_wider(d_plot, names_from = method, values_from = cluster))
#make key and add to d_plot
d_plot_wide$key <-gsub("sample_","", with(d_plot_wide,paste0(spot_name,"_",sample_name)))
#drop two columns we used to make the key
d_plot_wide$spot_name <- NULL
d_plot_wide$sample_name <-NULL
#match keys and reorder
#https://github.com/LieberInstitute/spatialLIBD/blob/master/R/cluster_import.R#L51-L64
merged_info <-
  merge(
    colData(spe),
    d_plot_wide,
    by = "key",
    sort = FALSE,
    all = TRUE
  )
m <- match(spe$key, merged_info$key)
merged_info <- merged_info[m, ]
spot_names <- rownames(colData(spe))

colData(spe) <- DataFrame(merged_info, check.names = FALSE)
colnames(spe) <- spot_names

##make plot
sample_ids <- unique(colData(spe)$sample_id)
cluster_colNames <- c("pseudobulk_PCA.y","pseudobulk_UMAP")
mycolors <- brewer.pal(7, "Dark2")

pdf(file = here::here("plots","vis_clus_semi_supervised_across_samples_2.pdf"))
for (i in seq_along(sample_ids)){
  for(j in seq_along(cluster_colNames)){
    my_plot <- vis_clus(
      spe = spe,
      clustervar = cluster_colNames[j],
      sampleid = sample_ids[i],
      colors =  mycolors,
      ... = paste0(" ",cluster_colNames[j])
    )
    print(my_plot)
  }
  
}
dev.off()

with(colData(spe),addmargins(table(spatial.cluster,pseudobulk_PCA.y,sample_id)))


