#tutorial from: https://edward130603.github.io/BayesSpace/articles/joint_clustering.html

#  qrsh -l bluejay,mem_free=20G,h_vmem=20G

suppressPackageStartupMessages({
  library("spatialLIBD")
  library("here")
  library(ggplot2)
  library(patchwork)
  library(scater)
  library(harmony)
  library(BayesSpace)
  library(scran)
})

load(file=here::here("processed-data","rdata", "spe", "for_talk","spe_snn_clusters_pc50_k10_101421.Rdata"))
dim(spe)
# [1]  28974 118003

spe_sub <- spe[,colData(spe)$sample_id == "Br2743_ant" | colData(spe)$sample_id == "Br6423_ant" | colData(spe)$sample_id == "Br3942_post" |colData(spe)$sample_id == "Br8492_post"]
dim(spe_sub)
# [1] 28974 16989
table(spe_sub$sample_id)
# Br2743_ant Br3942_post  Br6423_ant Br8492_post
# 4065        4400        3906        4618


spe_sub = spatialPreprocess(spe_sub, n.PCs = 50)
spe_sub = runUMAP(spe_sub, dimred = "PCA")

colnames(reducedDim(spe_sub, "UMAP")) = c("UMAP1", "UMAP2")

pdf(file=here::here("plots", "UMAP_dlpfc_spe.pdf"))
ggplot(data.frame(reducedDim(spe_sub, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe_sub$sample_id))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()
dev.off()

pdf(file=here::here("plots", "UMAP_dlpfc_spe.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()
dev.off()

# install.packages("devtools")
# devtools::install_github("immunogenomics/harmony")

spe_sub = RunHarmony(spe_sub, "sample_id", verbose = F)
spe_sub = runUMAP(spe_sub, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe_sub, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

pdf(file=here::here("plots", "UMAP_harmony_spe_sub.pdf"))
ggplot(data.frame(reducedDim(spe_sub, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe_sub$sample_id))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()
dev.off()

pdf(file=here::here("plots", "UMAP_harmony_spe.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()
dev.off()

save(spe_sub, file = here::here("processed-data","rdata","spe","spe_sub_harmony_dim.Rdata"))

###graph-based clustering with PC dimensions
Sys.time()
g_k10 <- buildSNNGraph(spe_sub, k = 10, use.dimred = 'PCA')
Sys.time()
#took a few minutes
save(g_k10, file=here::here("processed-data", "rdata","spe","g_k10_pca.Rdata"))

Sys.time() #12:55
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time() #1:09
save(g_walk_k10, file = here::here("processed-data", "rdata","spe","g_walk_k10_pca.Rdata"))

clust_k10 <- sort_clusters(g_walk_k10$membership)
save(clust_k10, file = here::here("processed-data", "rdata","spe","clust_k10_pca.Rdata"))

clust_k5_list <- lapply(4:28, function(n) {
  message(paste(Sys.time(), 'n =', n))
  sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)
save(clust_k5_list, file = here::here("processed-data", "rdata", "clust_k5_pca.Rdata"))

col.names <- paste0("SNN_k10_k",4:28)
for (i in seq_along(col.names)){
  colData(spe_sub) <- cbind(colData(spe_sub),clust_k5_list[i])
}
colnames(colData(spe_sub))[19:43] <- col.names

save(spe_sub, file = here::here("processed-data", "rdata","spe", "spe_sub_snn_clusters_pca.Rdata"))

pdf(file = here::here("plots","UMAP_pc50_k10_pca.pdf"))
for(i in seq_along(col.names)){
  myplot <- plotReducedDim(spe_sub, "UMAP", colour_by=col.names[i])
  print(myplot)
}
dev.off()

###graph-based clustering with harmony dimensions
Sys.time()
g_k10 <- buildSNNGraph(spe, k = 10, use.dimred = 'HARMONY')
Sys.time()
#took a few minutes
save(g_k10, file=here::here("processed-data", "rdata","spe","g_k10_harmony.Rdata"))

Sys.time() #12:55
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time() #1:09
save(g_walk_k10, file = here::here("processed-data", "rdata","spe","g_walk_k10_harmony.Rdata"))

clust_k10 <- sort_clusters(g_walk_k10$membership)
save(clust_k10, file = here::here("processed-data", "rdata","spe","clust_k10_harmony.Rdata"))

clust_k5_list <- lapply(4:28, function(n) {
  message(paste(Sys.time(), 'n =', n))
  sort_clusters(igraph::cut_at(g_walk_k10, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)
save(clust_k5_list, file = here::here("processed-data", "rdata", "clust_k5_harmony.Rdata"))

col.names <- paste0("SNN_k10_k",4:28)
for (i in seq_along(col.names)){
  colData(spe) <- cbind(colData(spe),clust_k5_list[i])
}
colnames(colData(spe))[19:43] <- col.names

save(spe, file = here::here("processed-data", "rdata","spe", "spe_snn_clusters_harmony.Rdata"))

##
pdf(file = here::here("plots","UMAP_pc50_k10_harmony.pdf"))
for(i in seq_along(col.names)){
  myplot <- plotReducedDim(spe_sub, "UMAP.HARMONY", colour_by=col.names[i])
  print(myplot)
}
dev.off()

##start here on 10/5/21
summary(spatialData(spe_sub)$array_row)
summary(spatialData(spe_sub)$array_col)

auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <-unique(spe$sample_id)
spe$row <- spatialData(spe)$array_row + auto_offset_row[spe$sample_id]
spe$col <- spatialData(spe)$array_col

summary(spe_sub$row)
summary(spe_sub$col)


pdf(file=here::here("plots", "bayesSpace_offset_check_dlpfc.pdf"),height = 7*length(unique(spe$sample_id)))
clusterPlot(spe, "sample_id", color = NA) + #make sure no overlap between samples
  labs(fill = "Sample", title = "Offset check")
dev.off()

save(spe_sub, file=here::here("processed-data","rdata", "spe", "spe_sub.Rdata"))

Sys.time() #4:30pm
spe = spatialCluster(spe, use.dimred = "HARMONY", q = 7, nrep = 10000) #use HARMONY
Sys.time()

pdf(file=here::here("plots","bayesSpace_clusterPlot_dlpfc.pdf"))
clusterPlot(spe, color = NA) + #plot clusters
  labs(title = "BayesSpace joint clustering")
dev.off()

save(spe, file=here::here("processed-data","rdata", "spe", "spe_bayesSpace.Rdata"))

