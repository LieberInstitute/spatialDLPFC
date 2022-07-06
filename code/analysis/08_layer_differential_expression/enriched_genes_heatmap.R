# code from Anthony https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/code/Pseudobulk/Plot_adult_enrichmentDE.R

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(dplyr)
  library(sessioninfo)
  library(tidyr)
})

# Load SPE
spe_pseudo <- readRDS(here::here("processed-data","rdata","spe","pseudo_bulked_spe", "spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k9.RDS"))

# Load modeling results
load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression","parsed_modeling_results_k9.Rdata"))

# Get mean expression
mat <- assays(spe_pseudo)$logcounts

# Extract the p-values
pvals <- modeling_results$enrichment[,grep("p_value_", colnames(modeling_results$enrichment))]
rownames(pvals) = rownames(mat)

# Extract the t-statistics
t_stat <- modeling_results$enrichment[,grep("[f|t]_stat_", colnames(modeling_results$enrichment))]
rownames(t_stat) = rownames(mat)

# Extract the FDRs
fdrs <- modeling_results$enrichment[,grep("fdr", colnames(modeling_results$enrichment))]
rownames(fdrs) = rownames(mat)

#get ensembl ids of mitochondrial genes
drop_mt <- rownames(rowData(spe_pseudo)[grep("^MT",rowData(spe_pseudo)$gene_name),])
length(drop_mt)
# [1] 13

dim(spe_pseudo)
# [1] 12225   268
spe_pseudo <- spe_pseudo[!(rownames(spe_pseudo)%in%drop_mt),]
dim(spe_pseudo)
# [1] 12212   268

dim(mat)
# [1] 12225   268
mat <- mat[!(rownames(mat)%in%drop_mt),]
dim(mat)
# [1] 12212   268

pvals <- pvals[!(rownames(pvals)%in%drop_mt),]
t_stat <- t_stat[!(rownames(t_stat)%in%drop_mt),]
fdrs <- fdrs[!(rownames(fdrs)%in%drop_mt),]


### pick top 10 genes per cluster:sample
cluster_specific_indices = mapply(function(t, p, f) {
  oo = order(t, decreasing = TRUE)[1:5]
},
as.data.frame(t_stat),
as.data.frame(pvals),
as.data.frame(fdrs))

#cluster_specific_indices <- gather(as.data.frame(cluster_specific_indices),t_stat,index,t_stat_1:t_stat_9,factor_key = TRUE)
#turn it into a long dataframe with a column for index and a column for cluster
cluster_ind = unique(as.numeric(cluster_specific_indices))

# Add logcounts from indexed from top genes
exprs_heatmap <- assays(spe_pseudo)[[2]][cluster_ind,]
rownames(exprs_heatmap) <- rowData(spe_pseudo)$gene_name[cluster_ind]

# Add column/sample annotations for pheatmap
cluster_labels <- as.vector(c(rep("Cluster_1", 28), rep("Cluster_2", 30), rep("Cluster_3", 30), rep("Cluster_4", 30),
                              rep("Cluster_5", 30), rep("Cluster_6", 30), rep("Cluster_7", 30), rep("Cluster_8", 30),
                              rep("Cluster_9",30)))

annotation_col <- data.frame(BayesSpace = factor(c(cluster_labels)))
rownames(annotation_col) = colnames(exprs_heatmap)
ann_colors = list(BayesSpace = brewer.pal(9, "Set1"))
names(ann_colors$BayesSpace) <- unique(annotation_col$BayesSpace)

# Add row/gene annotations
gene_labels <- as.vector(c(rep("Cluster_1", 5), rep("Cluster_2", 5), rep("Cluster_3", 5), rep("Cluster_4", 5),
                                     rep("Cluster_5", 5), rep("Cluster_6", 5), rep("Cluster_7", 5), rep("Cluster_8", 5),
                                     rep("Cluster_9",4)))
annotation_row <- data.frame(Genes = factor(c(gene_labels)))
rownames(annotation_row) = rownames(exprs_heatmap)

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots","08_layer_differential_expression","top5_enrichment_heatmap.pdf"), width = 8, height = 8)
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  color = inferno(20),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row = 9,
  main = "logcounts from enrichment model of 5 genes/cluster",
)
dev.off()

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots","08_layer_differential_expression","top5_enrichment_heatmap_clusRows.pdf"), width = 8, height = 8)
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = inferno(20),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row = 9,
  main = "logcounts from enrichment model of 5 genes/cluster",
)
dev.off()

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots","08_layer_differential_expression","top5_enrichment_heatmap_clusCols.pdf"), width = 8, height = 8)
pheatmap(
  exprs_heatmap,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  color = inferno(20),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row = 9,
  main = "logcounts from enrichment model of 5 genes/cluster",
)
dev.off()

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots","08_layer_differential_expression","top5_enrichment_heatmap_noClus.pdf"), width = 8, height = 8)
pheatmap(
  exprs_heatmap,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = inferno(20),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  annotation_row = annotation_row,
  fontsize_row = 9,
  main = "logcounts from enrichment model of 5 genes/cluster",
)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()