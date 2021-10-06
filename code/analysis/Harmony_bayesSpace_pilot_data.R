#tutorial from: https://edward130603.github.io/BayesSpace/articles/joint_clustering.html

library("spatialLIBD")


suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scater)
  library(harmony)
  library(BayesSpace)
})

sce <- fetch_data(type = "sce")
spe <- sce_to_spe(sce)
spe_sub <- spe[,colData(spe)$position == "0" & colData(spe)$replicate == "1"]
table(colData(spe_sub)$subject_position, colData(spe_sub)$replicate)

spe_sub = spatialPreprocess(spe_sub, n.PCs = 50)

spe_sub = runUMAP(spe_sub, dimred = "PCA")
colnames(reducedDim(spe_sub, "UMAP")) = c("UMAP1", "UMAP2")

pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/UMAP_pilot_data.pdf")
ggplot(data.frame(reducedDim(spe_sub, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe_sub$subject_position))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()
dev.off()

# install.packages("devtools")
# devtools::install_github("immunogenomics/harmony")

spe_sub = RunHarmony(spe_sub, "subject_position", verbose = F)
spe_sub = runUMAP(spe_sub, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe_sub, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/UMAP_harmony_pilot_data.pdf")
ggplot(data.frame(reducedDim(spe_sub, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe_sub$subject_position))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()
dev.off()

summary(spatialData(spe_sub)$array_row)
summary(spatialData(spe_sub)$array_col)
summary(colData(spe_sub)$row)
summary(colData(spe_sub)$col)

spe_sub$row[spe_sub$subject_position == "Br5292_pos0"] = 
  100 + spatialData(spe)$array_row[spe_sub$subject_position == "Br5292_pos0"]
spe_sub$col[spe_sub$subject_position == "Br5292_pos0"] = 
  spatialData(spe)$array_col[spe_sub$subject_position == "Br5292_pos0"]

spe_sub$row[spe_sub$subject_position == "Br5595_pos0"] = 
  spatialData(spe_sub)$array_row[spe_sub$subject_position == "Br5595_pos0"]
spe_sub$col[spe_sub$subject_position == "Br5595_pos0"] = 
  150 + spatialData(spe_sub)$array_col[spe_sub$subject_position == "Br5595_pos0"]

spe_sub$row[spe_sub$subject_position == "Br8100_pos0"] = 
  spatialData(spe_sub)$array_row[spe_sub$subject_position == "Br8100_pos0"]
spe_sub$col[spe_sub$subject_position == "Br8100_pos0"] = 
  spatialData(spe_sub)$array_col[spe_sub$subject_position == "Br8100_pos0"]

pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/bayesSpace_offset_check_pilot.pdf")
clusterPlot(spe_sub, "subject_position", color = NA) + #make sure no overlap between samples
  labs(fill = "Sample", title = "Offset check")
dev.off()

save(spe_sub, file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/pilot_dlpfc_data/spe_sub_100121.pdf")

spe_sub = spatialCluster(spe_sub, use.dimred = "HARMONY", q = 7, nrep = 10000) #use HARMONY

pdf(file="/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/bayesSpace_clusterPlot_pilot.pdf")
clusterPlot(spe_sub, color = NA) + #plot clusters
  labs(title = "BayesSpace joint clustering")
dev.off()
