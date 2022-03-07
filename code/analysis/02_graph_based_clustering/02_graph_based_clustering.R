#perform within and across?

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

#load spe object
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"))

###########################################
#graph-based clustering, no batch correction
##########################################
Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = 'PCA')
Sys.time()
# [1] "2021-12-16 01:09:06 EST"
#  "2021-12-16 01:24:11 EST"

Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
# "2021-12-16 01:27:49 EST"
# [1] "2021-12-16 09:26:14 EST"

clust_k10 <- sort_clusters(g_walk_k10$membership)

clust_k5_list <- lapply(4:28, function(n) {
  message(paste(Sys.time(), 'n =', n))
  sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)
save(clust_k5_list,file = here::here("processed-data","rdata","spe","02_graph_based_clustering","clust_k5_list_PCA.Rdata"))

## Add clusters to spe colData
cluster_colNames <- paste0("PCA_SNN_k10_k",4:28)
for (i in seq_along(cluster_colNames)){
  colData(spe) <- cbind(colData(spe),clust_k5_list[i])
}
colnames(colData(spe))[33:57] <- cluster_colNames  ##check here column range

cluster_export(
  spe,
  "PCA_SNN_k10_k7",
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

##make plot
nb.cols <- 25
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
sample_ids <- unique(colData(spe)$sample_id)
pdf(file = here::here("plots","02_graph_based_clustering","vis_clus_graph_based_pca.pdf"))
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

############################
##graph-based on batch corrected
############################


Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = 'HARMONY')
Sys.time()
#"2021-12-16 16:17:31 EST"
#"2021-12-16 16:34:28 EST"

Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
#[1] "2021-12-17 13:05:59 EST"

clust_k10 <- sort_clusters(g_walk_k10$membership)

clust_k5_list <- lapply(4:28, function(n) {
  message(paste(Sys.time(), 'n =', n))
  sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)

save(clust_k5_list,file = here::here("processed-data","rdata","spe","02_graph_based_clustering","clust_k5_list_Harmony.Rdata"))

## Add clusters to spe colData
cluster_colNames <- paste0("Harmony_SNN_k10_k",4:28)
for (i in seq_along(cluster_colNames)){
  colData(spe) <- cbind(colData(spe),clust_k5_list[i])
}
colnames(colData(spe))[35:59] <- cluster_colNames ### this needs to be edited

cluster_export(
  spe,
  "Harmony_SNN_k10_k7",
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

##make plot
sample_ids <- unique(colData(spe)$sample_id)
pdf(file = here::here("plots","02_graph_based_clustering","vis_clus_graph_based_harmony.pdf"))
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


######################################3
#graph-based clustering within samples
#####################################
#adapted from Luka's script https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/SpatialDE_clustering.Rmd

sample_names <- paste0("sample_", unique(colData(spe)$sample_id))
sample_names
load(file = here::here("processed-data","rdata","spe","01_build_spe","top.hvgs_all.Rdata"))
top.hvgs

genes<- rowData(spe)[rownames(rowData(spe))%in%top.hvgs,c("gene_id", "gene_name")]

dim(genes)
#[1] 2016    2

##Clustering
# -clustering on top 50 PCs 
# -clustering on top 10 UMAPs 

# parameters
n_umap <- 10
max_spatial <- 1
n_neighbors <- 10
n_clus <- 7 #changed from 8
d_plot <- data.frame()

#filter out genes that aren't in data, doesn't remove any genes
genes <- genes[which(genes$gene_id %in% rownames(spe)),]
dim(genes)

d_plot <- data.frame()
# run once per sample
for (i in seq_along(sample_names)) {
  # select spots from this sample
  spe_sub <- spe[, colData(spe)$sample_id == gsub("^sample_", "", sample_names[i])]
  dim(spe_sub)
  
  # ---------------------------------------------
  # extract and calculate features (PCA and UMAP)
  # ---------------------------------------------

  
  # run PCA on each sample
  logcounts_gb <- logcounts(spe_sub[genes$gene_id, ])
  
  # note: use 'prcomp' instead of 'calculatePCA' due to small number of genes
  out_pca <- prcomp(t(as.matrix(logcounts_gb)))$x[, 1:50]
  
  dims_PCA <- out_pca
  rownames(dims_PCA) <- colnames(spe_sub)
  dim(dims_PCA)
  stopifnot(nrow(dims_PCA) == ncol(spe_sub))
  
  
  # run UMAP 
  set.seed(1234)
  out_umap <- umap(dims_PCA, scale = TRUE, n_components = n_umap)
  
  dims_UMAP <- out_umap
  colnames(dims_UMAP) <- paste0("UMAP", seq_len(n_umap))
  rownames(dims_UMAP) <- colnames(spe_sub)
  dim(dims_UMAP)
  stopifnot(nrow(dims_UMAP) == ncol(spe_sub))
  
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
  
  d_plot <- rbind(d_plot, run_clustering(dims_PCA, method = "graph_based_PCA_within"))
  
}

save(d_plot, file = here::here("processed-data","rdata","spe","02_graph_based_clustering","d_plot_harmony_within.Rdata"))

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

cluster_export(
  spe,
  "graph_based_PCA_within",
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

##make plot
sample_ids <- unique(colData(spe)$sample_id)
cluster_colNames <- c("graph_based_PCA_within")
mycolors <- brewer.pal(7, "Dark2")

pdf(file = here::here("plots","02_graph_based_clustering","vis_clus_graph_based_within_samples.pdf"))
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

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

