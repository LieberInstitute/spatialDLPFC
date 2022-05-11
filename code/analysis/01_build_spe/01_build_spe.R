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

## reading the data
library("SpatialExperiment")
library("rtracklayer")


## vis
library("spatialLIBD")
library("RColorBrewer")
library(ggplot2)

## analysis
library("scran") ## requires uwot for UMAP
library("scater")
library("BiocParallel")
library("PCAtools")
library(harmony)


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

#[1] "2022-03-01 11:22:49 EST"

## Add the experimental information
spe$key <- paste0(colnames(spe), '_', spe$sample_id)
spe$subject <- sample_info$subjects[match(spe$sample_id, sample_info$sample_id)]
spe$region <- sample_info$regions[match(spe$sample_id, sample_info$sample_id)]
spe$sex <- sample_info$sex[match(spe$sample_id, sample_info$sample_id)]
spe$age <- sample_info$age[match(spe$sample_id, sample_info$sample_id)]
spe$diagnosis <- sample_info$diagnosis[match(spe$sample_id, sample_info$sample_id)]
spe$sample_id_complete <- spe$sample_id
spe$sample_id <- gsub("_2", "", spe$sample_id)


## Add information used by spatialLIBD
is_mito <- which(seqnames(spe) == "chrM")
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi

## Read in cell counts and segmentation results
segmentations_list <- lapply(sample_info$sample_id, function(sampleid) {
  current<-sample_info$sample_path[sample_info$sample_id==sampleid]
  file <- file.path(current, "spatial", "tissue_spot_counts.csv")
  if(!file.exists(file)) return(NULL)
  x <- read.csv(file)
  x$key <- paste0(x$barcode, "_", sampleid)
  return(x)
})
## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <- Reduce(function(...) merge(..., all = TRUE), segmentations_list[lengths(segmentations_list) > 0])

## Add the information
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <- segmentations[segmentation_match, - which(colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")), drop=FALSE]
colData(spe) <- cbind(colData(spe), segmentation_info)

mean(spe$count)
# [1] 5.092243

pdf(here::here("plots","01_build_spe","cells_per_spot.pdf"))
boxplot(spe$count)
dev.off()

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
save(spe_raw, file = here::here("processed-data", "rdata", "spe", "01_build_spe","spe_raw_final.Rdata"))

#remove genes that are not expressed in any spots
spe <- spe_raw[-no_expr, ]

dim(spe)
# 28916 149757

## Work with SPE
spe <- spe[, which(colData(spe_raw)$in_tissue)]
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
  pdf = here::here("plots", "01_build_spe","in_tissue_grid.pdf"),
  sort_clust = FALSE,
  colors = c("TRUE" = "grey90", "FALSE" = "orange")
)

summary(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    66.0   170.0   291.6   326.0 39053.0 

head(table(spe_raw$sum_umi[which(!colData(spe_raw)$in_tissue)]))
# 0  1  2  3  4  5 
# 12 17 50 59 65 81 

vis_grid_gene(
  spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
  geneid = "sum_umi",
  pdf = here::here("plots", "01_build_spe","out_tissue_sum_umi_all.pdf"),
  assayname = "counts"
)

summary(spe_raw$sum_gene[which(!colData(spe_raw)$in_tissue)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    55.0   133.0   206.6   249.0  7956.0 

vis_grid_gene(
  spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
  geneid = "sum_gene",
  pdf = here::here("plots", "01_build_spe", "out_tissue_sum_gene_all.pdf"),
  assayname = "counts"
)

summary(spe_raw$expr_chrM_ratio[which(!colData(spe_raw)$in_tissue)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.1500  0.1947  0.2072  0.2500  1.0000      12 

vis_grid_gene(
  spe = spe_raw[, which(!colData(spe_raw)$in_tissue)],
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "01_build_spe","out_tissue_expr_chrM_ratio_all.pdf"),
  assayname = "counts"
)

summary(spe$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    1358    2310    2674    3540   46881 
vis_grid_gene(
  spe = spe,
  geneid = "sum_umi",
  pdf = here::here("plots", "01_build_spe","in_tissue_sum_umi_all.pdf"),
  assayname = "counts"
)

summary(spe$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1     910    1426    1508    1999    8344 

vis_grid_gene(
  spe = spe,
  geneid = "sum_gene",
  pdf = here::here("plots", "01_build_spe","in_tissue_sum_gene_all.pdf"),
  assayname = "counts"
)

summary(spe$expr_chrM_ratio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.07279 0.10440 0.12106 0.15513 1.00000 
vis_grid_gene(
  spe = spe,
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots","01_build_spe", "in_tissue_expr_chrM_ratio_all.pdf"),
  assayname = "counts"
)

vis_grid_gene(
  spe = spe_raw,
  geneid = "sum_umi",
  pdf = here::here("plots", "01_build_spe","all_sum_umi.pdf"),
  assayname = "counts"
)

vis_grid_gene(
  spe = spe_raw,
  geneid = "sum_gene",
  pdf = here::here("plots","01_build_spe", "all_sum_gene.pdf"),
  assayname = "counts"
)
vis_grid_gene(
  spe = spe_raw,
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots","01_build_spe", "all_expr_chrM_ratio.pdf"),
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
## Compute the low library size thresholds in a sample-aware way
spe$scran_low_lib_size <-
  factor(
    isOutlier(
      spe$sum_umi,
      type = "lower",
      log = TRUE,
      batch = spe$sample_id
    ),
    levels = c("TRUE", "FALSE")
  )
spe$scran_low_n_features <-
  factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
spe$scran_high_subsets_Mito_percent <-
  factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))
#save
save(spe, file = here::here("processed-data", "rdata", "spe", "01_build_spe","spe_final.Rdata"))


for(i in colnames(qcfilter)) {
  vis_grid_clus(
    spe = spe,
    clustervar = paste0("scran_", i),
    pdf = here::here("plots","01_build_spe", paste0("scran_", i, ".pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange")
  )
}

## Sample-aware low-library size
vis_grid_clus(
  spe = spe,
  clustervar = "scran_low_lib_size",
  pdf = here("plots", "01_build_spe",paste0("vis_clus_sample_aware_low_lib_size.pdf")),
  sort_clust = FALSE,
  colors = c("FALSE" = "grey90", "TRUE" = "orange"),
  spatial = FALSE,
  point_size = 2,
  height = 24,
  width = 90
)

#edit low library size plot 
sample_order = unique(spe$sample_id)

p_scran <-vis_grid_clus(
  spe = spe,
  clustervar = "scran_low_lib_size",
  sort_clust = FALSE,
  colors = c("FALSE" = "grey90", "TRUE" = "orange"),
  spatial = FALSE,
  point_size = 2,
  height = 24,
  width = 90, 
  return_plots = TRUE,
  image_id = "lowres"
)

p_scran <- lapply(sample_order, function(sampleid){
  p <- p_scran[[sampleid]]
  p + theme(legend.position = "none")
})
names(p_scran) <- sample_order

pdf(file = here::here("plots", "01_build_spe",paste0("vis_clus_sample_aware_low_lib_size_sfigur.pdf")), height = 24, width = 36)
print(cowplot::plot_grid(plotlist = p_scran))
dev.off()


#discard low library size and save spe
dim(spe)
# [1]  36601 149757

## Remove genes with no data
low_lib_remove <- which(colData(spe)$scran_low_lib_size == TRUE)

length(low_lib_remove)
# [1] 7685

#remove remove spots that have low lib size
spe <- spe[,-low_lib_remove]


#look at spots with large cell counts
## Let's look at cluster 12 from BayesSpace k = 15
spe$large_cell_count <- factor(spe$count > 50,
                                 levels = c("TRUE", "FALSE"))
vis_grid_clus(
  spe = spe,
  clustervar = "large_cell_count",
  pdf = here("plots",  "01_build_spe",paste0("vis_clus_large_cell_count.pdf")),
  sort_clust = FALSE,
  colors = c("FALSE" = "grey90", "TRUE" = "orange"),
  spatial = FALSE,
  point_size = 2,
  height = 24,
  width = 90
)

p_count <- vis_grid_gene(
  spe = spe,
  geneid = "count",
  assayname = "counts",
  point_size = 2,
  height = 48,
  width = 90,
  return_plots = TRUE
)
# p_count <- lapply(sample_order, function(sampleid){
#   p <- p_count[[sampleid]]
#   p + scale_color_gradientn(colors = viridisLite::viridis(21),limits=c(0,300))
# })
names(p_count) <- sample_order

pdf(file = here::here("plots", "01_build_spe", paste0("vis_clus_cell_count.pdf")), height = 24, width = 36)
print(cowplot::plot_grid(plotlist = p_count))
dev.off()

save(spe, file = here::here("processed-data", "rdata", "spe", "01_build_spe","spe_filtered_final.Rdata"))

##stopping here
## Find quick clusters
set.seed(030122)
Sys.time()
spe$scran_quick_cluster <- quickCluster(
  spe,
  BPPARAM = MulticoreParam(4),
  block = spe$sample_id,
  block.BPPARAM = MulticoreParam(4)
)
Sys.time()


Sys.time()
spe <-
  computeSumFactors(spe,
                    clusters = spe$scran_quick_cluster,
                    BPPARAM = MulticoreParam(4)
  )
Sys.time()


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
  here::here("plots", "01_build_spe","scran_modelGeneVar_final.pdf"),
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
     file = here::here("processed-data", "rdata", "spe","01_build_spe", "top.hvgs_all.Rdata"))

set.seed(020122)
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = 50)
Sys.time()

dim(reducedDim(spe, "PCA"))
# 118793     50

#make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(spe,"PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow
# 50

pdf(
  here::here("plots","01_build_spe", "pca_elbow_final.pdf"),
  useDingbats = FALSE
)
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")
dev.off()


summary(apply(reducedDim(spe, "PCA"), 2, sd))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9184  0.9343  0.9509  1.1872  1.0561  4.1944 

# RunUMAP
set.seed(030122)
Sys.time()
spe = runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) = c("UMAP1", "UMAP2")
Sys.time()

#Run TSNE
set.seed(030122)
Sys.time()
spe <-
  runTSNE(spe,
          dimred = "PCA",
          name = "TSNE_perplexity80",
          perplexity = 80
  )
Sys.time()
#[1] "2021-12-15 15:30:48 EST"
#[1] "2021-12-15 16:49:28 EST"


#make plots of UMAP
pdf(file=here::here("plots", "01_build_spe","UMAP_subject.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe$subject))) +
  geom_point() +
  labs(color = "Subject") +
  theme_bw()
dev.off()

pdf(file=here::here("plots", "01_build_spe","UMAP_sample_id.pdf"))
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


pdf(file=here::here("plots", "01_build_spe","UMAP_harmony_sample_id.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")),
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()
dev.off()


###########for bayesSpace

##do offset so we can run BayesSpace
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <-unique(spe$sample_id)
spe$row <- colData(spe)$array_row + auto_offset_row[spe$sample_id]
spe$col <- colData(spe)$array_col

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message("Running spatialCluster()")

Sys.time()
save(spe, file = here::here("processed-data","rdata", "spe", "01_build_spe","spe_filtered_final.Rdata"))
Sys.time()

# library(spatialLIBD)
# load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final.Rdata")
# spe <- cluster_import(
#   spe,
#   cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"), 
#   prefix = ""
# )
# save(spe, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

