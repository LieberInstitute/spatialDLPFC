#adapted from https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/SpatialDE_clustering.Rmd

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

library(here)
library(SpatialExperiment)

##load spe object
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/spe_tiny.Rdata")

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

# function to sort cluster labels by descending frequency (from Leo)
sort_clusters <- function(clusters, map_subset = NULL) {
  if (is.null(map_subset)) {
    map_subset <- rep(TRUE, length(clusters))
  }
  map <- rank(length(clusters[map_subset]) - table(clusters[map_subset]), ties.method = "first")
  res <- map[clusters]
  factor(res)
}

# run once per sample
for (i in seq_along(sample_names)) {
  # select spots from this sample
  spe_sub <- spe[, colData(spe)$sample_id == gsub("^sample_", "", sample_names[i])]
  dim(spe_sub)
  
  # ---------------------------------------------
  # extract and calculate features (PCA and UMAP)
  # ---------------------------------------------
  
  ### pseudobulk layer genes (from Leo's analyses; 198 genes)
  
  #filter out pseudobulk genes that aren't in data, removes 1 gene
  genes_pseudobulk <- genes_pseudobulk[which(genes_pseudobulk$gene_id %in% rownames(spe_sub)),]
  
  # run PCA on pseudobulk layer genes (from Leo's analyses; 198 genes)
  logcounts_pseudobulk <- logcounts(spe_sub[genes_pseudobulk$gene_id, ])
  
  # note: use 'prcomp' instead of 'calculatePCA' due to small number of genes
  out_pca_pseudobulk <- prcomp(t(as.matrix(logcounts_pseudobulk)))$x[, 1:50]
  
  dims_pseudobulk_PCA <- out_pca_pseudobulk
  rownames(dims_pseudobulk_PCA) <- colnames(spe_sub)
  dim(dims_pseudobulk_PCA)
  #[1] 298  50
  stopifnot(nrow(dims_pseudobulk_PCA) == ncol(spe_sub))
  
  
  # run UMAP on pseudobulk layer genes (from Leo's analyses; 197 genes)
  set.seed(1234)
  out_umap_pseudobulk <- umap(dims_pseudobulk_PCA, scale = TRUE, n_components = n_umap)
  
  dims_pseudobulk_UMAP <- out_umap_pseudobulk
  colnames(dims_pseudobulk_UMAP) <- paste0("UMAP", seq_len(n_umap))
  rownames(dims_pseudobulk_UMAP) <- colnames(spe_sub)
  dim(dims_pseudobulk_UMAP)
  stopifnot(nrow(dims_pseudobulk_UMAP) == ncol(spe_sub))
  
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
      sample_name = as.character(sample_names[i]), 
      method = as.character(method), 
      cluster = as.numeric(clus),  
      stringsAsFactors = FALSE
    )
  }
  
  d_plot <- rbind(d_plot, run_clustering(dims_pseudobulk_PCA, method = "pseudobulk_PCA"))
  d_plot <- rbind(d_plot, run_clustering(dims_pseudobulk_UMAP, method = "pseudobulk_UMAP"))
  
}





