library(SpatialExperiment)
library(here)
library(pheatmap)
library(spatialLIBD)

#load spe object
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)

#load clusters
spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"), 
  prefix = ""
)

#make table of bayesSpace k = 9 vs k = 28
my.table <- with(colData(spe), table(bayesSpace_harmony_9, bayesSpace_harmony_28))

pdf(file = here::here("plots","07a_bad_clusters","clusters_k9_k28.pdf"))
pheatmap(
  my.table,
  show_rownames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)
dev.off()

my.table.log <- with(colData(spe), log10(table(bayesSpace_harmony_9, bayesSpace_harmony_28)+1))

pdf(file = here::here("plots","07a_bad_clusters","clusters_k9_k28_log.pdf"))
pheatmap(
  my.table.log,
  show_rownames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)
dev.off()


