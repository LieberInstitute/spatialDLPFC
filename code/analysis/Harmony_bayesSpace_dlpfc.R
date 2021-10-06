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
})

load(file=here::here("processed-data","rdata", "spe", "spe_snn_clusters_pc50_k10_092221.Rdata"))
dim(spe)
# [1]  28974 118003

spe_sub <- spe[,colData(spe)$sample_id == "Br2743_ant" | colData(spe)$sample_id == "Br6423_ant" | colData(spe)$sample_id == "Br3942_post" |colData(spe)$sample_id == "Br8492_post"]
dim(spe_sub)
# [1] 28974 16989
table(colData(spe_sub)$sample_id)
# Br2743_ant Br3942_post  Br6423_ant Br8492_post
# 4065        4400        3906        4618


spe_sub = spatialPreprocess(spe_sub, n.PCs = 50)
spe_sub = runUMAP(spe_sub, dimred = "PCA")

colnames(reducedDim(spe_sub, "UMAP")) = c("UMAP1", "UMAP2")

pdf(file=here::here("plots", "UMAP_dlpfc_spe_sub.pdf"))
ggplot(data.frame(reducedDim(spe_sub, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe_sub$sample_id))) +
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

save(spe_sub, file=here::here("processed-data","rdata", "spe", "spe_sub.Rdata"))
##start here on 10/5/21
summary(spatialData(spe_sub)$array_row)
summary(spatialData(spe_sub)$array_col)

auto_offset_row <- as.numeric(factor(unique(colData(spe_sub)$sample_id))) * 100
names(auto_offset_row) <-factor(unique(colData(spe_sub)$sample_id))
colData(spe_sub)$row <- colData(spe_sub)$array_row + auto_offset_row[as.numeric(factor(colData(spe_sub)$sample_id))]

spe_sub$row[spe_sub$sample_id == "Br2743_ant"] = 
  100 + spatialData(spe)$array_row[spe_sub$sample_id == "Br2743_ant"]
spe_sub$col[spe_sub$sample_id == "Br2743_ant"] = 
  spatialData(spe)$array_col[spe_sub$sample_id == "Br2743_ant"]

spe_sub$row[spe_sub$sample_id == "Br6423_ant"] = 
  200 + spatialData(spe_sub)$array_row[spe_sub$sample_id == "Br6423_ant"]
spe_sub$col[spe_sub$sample_id == "Br6423_ant"] =
  spatialData(spe_sub)$array_col[spe_sub$sample_id == "Br6423_ant"]

spe_sub$row[spe_sub$sample_id == "Br3942_post"] = 
  300 + spatialData(spe_sub)$array_row[spe_sub$sample_id == "Br3942_post"]
spe_sub$col[spe_sub$sample_id == "Br3942_post"] = 
  spatialData(spe_sub)$array_col[spe_sub$sample_id == "Br3942_post"]

spe_sub$row[spe_sub$sample_id == "Br8492_post"] = 
  400 + spatialData(spe_sub)$array_row[spe_sub$sample_id == "Br8492_post"]
spe_sub$col[spe_sub$sample_id == "Br8492_post"] = 
  spatialData(spe_sub)$array_col[spe_sub$sample_id == "Br8492_post"]

summary(colData(spe_sub)$row)
summary(colData(spe_sub)$col)


pdf(file=here::here("plots", "bayesSpace_offset_check_dlpfc.pdf"))
clusterPlot(spe_sub, "sample_id", color = NA) + #make sure no overlap between samples
  labs(fill = "Sample", title = "Offset check")
dev.off()

save(spe_sub, file=here::here("processed-data","rdata", "spe", "spe_sub.Rdata"))

Sys.time() #4:30pm
spe_sub = spatialCluster(spe_sub, use.dimred = "HARMONY", q = 7, nrep = 10000) #use HARMONY
Sys.time()

pdf(file=here::here("plots","bayesSpace_clusterPlot_dlpfc.pdf"))
clusterPlot(spe_sub, color = NA) + #plot clusters
  labs(title = "BayesSpace joint clustering")
dev.off()
