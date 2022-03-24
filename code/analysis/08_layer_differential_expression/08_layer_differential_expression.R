
library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(pheatmap)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load pseudobulked sce object
load(file = here::here("processed-data","rdata","spe","07_spatial_registration",paste0("sce_pseudobulk_bayesSpace_k",k,".Rdata")))

#load data (eb0_list_cell)
load(file = here::here("processed-data","rdata","spe","07_spatial_registration",paste0("dlpfc_pseudobulked_bayesSpace_specific_Ts_k",k,".Rdata")))

###############################
##### get mean expression  ####
###imported from 07_spatial_registration.R
mat <- assays(sce_pseudobulk_bayesSpace)$logcounts #make matrix of just the log normalized counts

## filter
gIndex = rowMeans(mat) > 0.2 # find the genes for which the mean expression is greater than 0.2
mat_filter = mat[gIndex, ] #subset matrix on just those genes.  want to remove lowly expressed genes. 


##########
## Extract the p-values
###imported from 07_spatial_registration.R
pvals0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cell) = rownames(mat_filter)

t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) = rownames(mat_filter)

### pick top 5 genes per cluster:sample
cluster_specific_indices = mapply(function(t, p) {
  oo = order(t, decreasing = TRUE)[1:5]
},
as.data.frame(t0_contrasts_cell),
as.data.frame(pvals0_contrasts_cell))
cluster_ind = unique(as.numeric(cluster_specific_indices))

exprs_heatmap <- assays(sce_pseudobulk_bayesSpace)[[2]][cluster_ind,]
rownames(exprs_heatmap) <- rowData(sce_pseudobulk_bayesSpace)$gene_name[cluster_ind]

pdf(file = here::here("plots","08_layer_differential_expression",paste0("enrichment_heatmap_k",k,".pdf")), width = k*5, height = k*1.5)
pheatmap(
  exprs_heatmap,
  show_rownames = TRUE,
  cluster_rows = FALSE
)
dev.off()

