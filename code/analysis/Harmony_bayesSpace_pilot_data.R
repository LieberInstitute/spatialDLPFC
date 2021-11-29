#tutorial from: https://edward130603.github.io/BayesSpace/articles/joint_clustering.html

library("spatialLIBD")
library("SpatialExperiment")
library("here")
library("scran")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scater)
  library(harmony)
  library(BayesSpace)
})

sce <- fetch_data(type = "sce")
spe <- sce_to_spe(sce)
#spe_sub <- spe[,colData(spe)$position == "0" & colData(spe)$replicate == "1"]
#table(colData(spe_sub)$subject_position, colData(spe_sub)$replicate)

spe = spatialPreprocess(spe, n.PCs = 50)
spe = runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) = c("UMAP1", "UMAP2")

pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/UMAP_pilot_data_subject.pdf")
ggplot(data.frame(reducedDim(spe, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe$subject))) +
  geom_point() +
  labs(color = "Subject") +
  theme_bw()
dev.off()

pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/UMAP_pilot_data_subject_position.pdf")
ggplot(data.frame(reducedDim(spe, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe$subject_position))) +
  geom_point() +
  labs(color = "Subject Position") +
  theme_bw()
dev.off()

pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/UMAP_pilot_data_sample_id.pdf")
ggplot(data.frame(reducedDim(spe, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()
dev.off()

# install.packages("devtools")
# devtools::install_github("immunogenomics/harmony")

spe = RunHarmony(spe, "sample_id", verbose = F)
spe = runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

save(spe, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/pilot_dlpfc_data/spe_pilot_102121.Rdata")

# pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/UMAP_harmony_pilot_data_subject.pdf")
# ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")), 
#        aes(x = UMAP1, y = UMAP2, color = factor(spe$subject))) +
#   geom_point() +
#   labs(color = "Subject") +
#   theme_bw()
# dev.off()

# pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/UMAP_harmony_pilot_data_subject_position.pdf")
# ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")), 
#        aes(x = UMAP1, y = UMAP2, color = factor(spe$subject_position))) +
#   geom_point() +
#   labs(color = "Subject Position") +
#   theme_bw()
# dev.off()
# 
pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/UMAP_harmony_pilot_data_sample_id.pdf")
ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")),
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()
dev.off()

summary(spatialData(spe)$array_row)
summary(spatialData(spe)$array_col)


auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <-unique(spe$sample_id)
spe$row <- spatialData(spe)$array_row + auto_offset_row[spe$sample_id]
spe$col <- spatialData(spe)$array_col

summary(colData(spe)$row)
summary(colData(spe)$col)

pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/bayesSpace_offset_check_pilot.pdf")
clusterPlot(spe, "subject", color = NA) + #make sure no overlap between samples
  labs(fill = "Subject", title = "Offset check")
dev.off()

spe = spatialCluster(spe, use.dimred = "HARMONY", q = 7, nrep = 10000) 
spe = spatialCluster(spe, use.dimred = "PCA", q = 7, nrep = 10000)#use HARMONY
#spe = spatialCluster(spe, use.dimred = "HARMONY", q = 14, nrep = 10000)
spe = spatialCluster(spe, use.dimred = "HARMONY", q = 8, nrep = 10000) 

#save(spe, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/pilot_dlpfc_data/spe_pilot_102121.Rdata")
save(spe, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/pilot_dlpfc_data/spe_pilot_bayesSpace8.Rdata")
save(spe, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/pilot_dlpfc_data/spe_pilot_bayesSpace_pcs.Rdata")


# pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/bayesSpace_clusterPlot_pilot.pdf")
# clusterPlot(spe, color = NA) + #plot clusters
#   labs(title = "BayesSpace joint clustering")
# dev.off()

pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/bayesSpace_clusterPlot_pilot_sample_id.pdf")
clusterPlot(spe, color = NA) + #plot clusters
  labs(title = "BayesSpace joint clustering")
dev.off()

###graph-based clustering on harmony dimensions
load(file=here::here("processed-data","rdata", "pilot_dlpfc_data", "spe_snn_clusters_pc50_k10_pilot.Rdata"))

Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = 'HARMONY')
Sys.time()
save(g_k10, file=here::here("processed-data","rdata", "pilot_dlpfc_data","g_pc50_k10__harmony_pilot.Rdata"))

Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
save(g_walk_k10, file = here::here("processed-data", "rdata","pilot_dlpfc_data","g_walk_pc50_k10_harmony_pilot.Rdata"))

clust_k10 <- sort_clusters(g_walk_k10$membership)
save(clust_k10, file = here::here("processed-data", "rdata","pilot_dlpfc_data","clust_pc50_k10_harmony_pilot.Rdata"))

### For the SNN graph with K = 50, find which nested subset best matches
## the clusters from 10x Genomics labeled by Kristen Maynard and Keri Martinowich
clust_k5_list <- lapply(4:28, function(n) {
  message(paste(Sys.time(), 'n =', n))
  sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)
save(clust_k5_list, file = here::here("processed-data", "rdata", "pilot_dlpfc_data","clust_k5_pc50_k10_harmony_pilot_list.Rdata"))

## Add clusters to spe colData
col.names <- paste0("batch_corr_SNN_k10_k",4:28)
for (i in seq_along(col.names)){
  colData(spe) <- cbind(colData(spe),clust_k5_list[i])
}
colnames(colData(spe))[72:96] <- col.names

save(spe, file = here::here("processed-data", "rdata","pilot_dlpfc_data", "spe_snn_clusters_pilot.Rdata"))




