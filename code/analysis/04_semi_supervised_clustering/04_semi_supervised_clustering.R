library(SingleCellExperiment)
library(scran)
library(scater)
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(uwot)
library(mclust)
library(aricode)
library(nlme)
library(ggplot2)
library(RColorBrewer)
library(Polychrome)
library(patchwork)
library(broom)
library(magick)
library(harmony)


## perform within and across?
load(file=here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"))

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


###do harmony batch correction for PCA calculated using pseudobulk genes
spe = RunHarmony(spe, "sample_id",reduction = "pseudobulk_PCA", verbose = F)

spe = runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")


pdf(file=here::here("plots","04_semi_supervised_clustering","UMAP_pseudobulk_harmony_sample_id.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")),
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()
dev.off()

Sys.time()
save(spe, file = here::here("processed-data","rdata", "spe","01_build_spe","spe_filtered_final_pseudobulk.Rdata"))
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

save(d_plot, file = here::here("processed-data","rdata","spe","clustering_results","d_plot_semi_supervised_across.Rdata"))

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
cluster_colNames <- c("pseudobulk_PCA","pseudobulk_UMAP")
mycolors <- brewer.pal(7, "Dark2")

pdf(file = here::here("plots","vis_clus_semi_supervised_across_samples.pdf"))
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

cluster_export(
  spe,
  "pseduobulk_PCA_across",
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

cluster_export(
  spe,
  "pseduobulk_PCA_across",
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

with(colData(spe),addmargins(table(spatial.cluster,pseudobulk_PCA.y,sample_id)))