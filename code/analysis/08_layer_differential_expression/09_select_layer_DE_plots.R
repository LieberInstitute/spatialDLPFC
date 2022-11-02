# library(edgeR)
library("SpatialExperiment")
library("scater")
# library(scran)
library("here")
library("sessioninfo")

plot_dir <- here("plots", "08_layer_differential_expression", "09_select_layer_DE_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

# data_dir <- here("processed-data", "rdata", "spe", "08_layer_differential_expression")
# if(!dir.exists(data_dir)) dir.create(data_dir)

## Load data
message(Sys.time(), " - Loading spe")
load(here("processed-data", "rdata","spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

## prep objet for plotting
rownames(spe) <- rowData(spe)$gene_name
spe$bayesSpace_harmony_9 <- as.factor(spe$bayesSpace_harmony_9)
spe$bayesSpace_harmony_16 <- as.factor(spe$bayesSpace_harmony_16)

source("my_plotExpression.R")

#### K9 violin plots ####
# load(file = here("processed-data", "rdata", "spe", "pseudo_bulked_spe","sce_pseudobulk_bayesSpace_k9.Rdata"), verbose = TRUE) ## file doesn't exist?
spe_k9 <- readRDS(file = here("processed-data", "rdata", "spe", "pseudo_bulked_spe","spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k9.RDS"))
rownames(spe_k9) <- rowData(spe_k9)$gene_name

## establish color pallet
k9_colors <- Polychrome::palette36.colors(9)
names(k9_colors) <- c(1:9)

k9_genes <- c("CLDN5", "TAGLN", "MYL9", "ACTA2", "SLC2A1", "HBA1", "EPAS1")

k9_plot <- my_plotExpression(spe_k9, genes = k9_genes, assay = "logcounts", cat = "BayesSpace", fill_colors = k9_colors)
ggsave(k9_plot, filename = here(plot_dir, "k9_expression_meninges.png"))

#### K9 violin plots ####
# load(file = here("processed-data", "rdata", "spe", "pseudo_bulked_spe","sce_pseudobulk_bayesSpace_k9.Rdata"), verbose = TRUE) ## file doesn't exist?
spe_k16 <- readRDS(file = here("processed-data", "rdata", "spe", "pseudo_bulked_spe","spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k16.RDS"))
rownames(spe_k16) <- rowData(spe_k16)$gene_name

## establish color pallet
k16_colors <- Polychrome::palette36.colors(16)
names(k16_colors) <- c(1:16)

## 1a vs. 1b
k16_genes_1ab <- c("SPARC", "MSX1", "RELN", "APOE")
k16_genes %in% rownames(spe_k16)

k16_1ab_plot <- my_plotExpression(spe_k16, genes = k16_genes, assay = "logcounts", cat = "BayesSpace", fill_colors = k16_colors)
ggsave(k16_1ab_plot, filename = here(plot_dir, "k16_expression_1a-1b.png"))

## 6a 
k16_6a_plot <- my_plotExpression(spe_k16, genes = c("SMIM32", "DACH1", "KIF1A", "GALNT14"), 
                                  assay = "logcounts", cat = "BayesSpace", fill_colors = k16_colors)

ggsave(k16_6a_plot, filename = here(plot_dir, "k16_expression_6a.png"))

## 6b 
k16_6b_plot <- my_plotExpression(spe_k16, genes = c("KRT17", "DIRAS2", "SEMA3E"), 
                                 assay = "logcounts", cat = "BayesSpace", fill_colors = k16_colors)

ggsave(k16_6b_plot, filename = here(plot_dir, "k16_expression_6b.png"))


