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

## vis
library("spatialLIBD")

## analysis
library("scran") ## requires uwot for UMAP
library("scater")
library("BiocParallel")
library("PCAtools")
library(harmony)
library(BayesSpace)

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
    "Round4/DLPFC_Br2720_ant_2",
    "Round2/DLPFC_Br2720_mid_manual_alignment",
    "Round2/DLPFC_Br2720_post_extra_reads",
    "Round4/DLPFC_Br6432_ant_2",
    "Round2/DLPFC_Br6432_mid_manual_alignment",
    "Round2/DLPFC_Br6432_post_manual_alignment",
    "Round3/DLPFC_Br6471_ant_manual_alignment_all",
    "Round3/DLPFC_Br6471_mid_manual_alignment_all",
    "Round3/DLPFC_Br6471_post_manual_alignment_all",
    "Round3/DLPFC_Br6522_ant_manual_alignment_all",
    "Round3/DLPFC_Br6522_mid_manual_alignment_all",
    "Round3/DLPFC_Br6522_post_manual_alignment_all",
    "Round3/DLPFC_Br8325_ant_manual_alignment_all",
    "Round4/DLPFC_Br8325_mid_2",
    "Round3/DLPFC_Br8325_post_manual_alignment_all",
    "Round3/DLPFC_Br8667_ant_extra_reads",
    "Round3/DLPFC_Br8667_mid_manual_alignment_all",
    "Round3/DLPFC_Br8667_post_manual_alignment_all"
    
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
  here::here("processed-data", "NextSeq"),
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

##"2021-12-07 12:06:15 EST"

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
# [1] 7468
length(no_expr) / nrow(spe) * 100
# [1] 20.40381

## Create two versions: one with and one without filtering by tissue spot
spe_raw <- spe

pryr::object_size(spe_raw)
# 24,754,498,680 B
dim(spe_raw)
# [1] 36601 149757

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
# 4,663,066,584 B

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
# 0.0    67.0   171.0   297.5   329.0 39053.0 

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
# 0      55     133     210     250    7956

vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "sum_gene",
  pdf = here::here("plots", "out_tissue_sum_gene_all_final.pdf"),
  assayname = "counts"
)

summary(spe_raw$expr_chrM_ratio[which(spatialData(spe_raw)$in_tissue=="FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.1497  0.1944  0.2068  0.2500  1.0000      12 

vis_grid_gene(
  spe = spe_raw[, which(spatialData(spe_raw)$in_tissue=="FALSE")],
  geneid = "expr_chrM_ratio",
  pdf = here::here("plots", "out_tissue_expr_chrM_ratio_all_all.pdf"),
  assayname = "counts"
)

summary(spe$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    1362    2312    2677    3542   46881 
vis_grid_gene(
  spe = spe,
  geneid = "sum_umi",
  pdf = here::here("plots", "in_tissue_sum_umi_all_final.pdf"),
  assayname = "counts"
)

summary(spe$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1     912    1427    1510    2000    8344 

vis_grid_gene(
  spe = spe,
  geneid = "sum_gene",
  pdf = here::here("plots", "in_tissue_sum_gene_all_final.pdf"),
  assayname = "counts"
)

summary(spe$expr_chrM_ratio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.07273 0.10433 0.12102 0.15511 1.00000 
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
# 4991                      5984                      3739 
# discard 
# 9127 


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
# [1]"2021-12-07 13:16:29 EST"
# [1] "2021-12-07 13:16:29 EST"

Sys.time()
#"2021-12-07 13:16:29 EST"
## Might be needed:
# options(error = recover)
spe <-
  computeSumFactors(spe,
                    clusters = spe$scran_quick_cluster,
                    BPPARAM = MulticoreParam(4)
  )
Sys.time()
#"2021-12-07 13:24:55 EST"


## Related to https://github.com/LTLA/scuttle/issues/7
# In .rescale_clusters(clust.profile, ref.col = ref.clust, min.mean = min.mean) :
#   inter-cluster rescaling factor for cluster 64 is not strictly positive,
# reverting to the ratio of average library sizes

table(spe$scran_quick_cluster)


summary(sizeFactors(spe))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000034  0.478415  0.847992  1.000000  1.330493 16.795995 

spe <- logNormCounts(spe)
pryr::object_size(spe)
#6,098,469,888 B



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
# [1]  2005

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# [1] 18299

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
# [1] 17712

save(top.hvgs,
     top.hvgs.fdr5,
     top.hvgs.fdr1,
     file = here::here("processed-data", "rdata", "spe", "top.hvgs_all_final.Rdata"))

#Run PCA with BayesSpace 
set.seed(20211207)
Sys.time()
spe = spatialPreprocess(spe, n.PCs = 50)
Sys.time()
#"2021-12-07 13:39:03 EST"
#"2021-12-07 14:07:02 EST"

dim(reducedDim(spe, "PCA"))
#118584     50

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
# 0.9192  0.9339  0.9503  1.1889  1.0560  4.2002 

# RunUMAP
spe = runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) = c("UMAP1", "UMAP2")

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
#"2021-12-07 14:57:57 EST"
#"2021-12-07 15:56:37 EST"

save(spe, file = here::here("processed-data","rdata", "spe", "spe_final.Rdata"))

#make plots of TSNE and UMAP
#Run harmony batch correction and produce plots 


