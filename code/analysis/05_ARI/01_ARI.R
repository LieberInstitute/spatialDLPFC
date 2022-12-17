library("SpatialExperiment")
library("mclust")
library("spatialLIBD")
library("tidyr")
library("ggplot2")
library("ggpubr")
library("viridis")
library("here")
library("sessioninfo")

# plot ARI for pilot data comparing different clustering algorithms to Kristen's manual annnotations

# load pilot data with non batch-corrected BayesSpace k - 7
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/pilot_dlpfc_data/05_ARI/spe_pilot_bayesSpace_pcs.Rdata")
spe.bspc <- spe

# load pilot data with batch-corrected BayesSpace k - 7
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/pilot_dlpfc_data/spe_pilot_bayesSpace_batch_corr_sampleID.Rdata")

spe$spatial.cluster.pcs <- spe.bspc$spatial.cluster


sample_ids <- unique(spe$sample_id)
ari.df <- data.frame(matrix(ncol = 5, nrow = 12))
row.names(ari.df) <- sample_ids
colnames(ari.df) <- c("sample_id", "SNN_k10_k7", "batch_corr_SNN_k10_k7", "bayesSpace_pc", "bayesSpace")

for (i in seq_along(sample_ids)) {
    spe_sub <- spe[, colData(spe)$sample_id == sample_ids[i]]
    ari.df$sample_id <- sample_ids[i]
    ari.df[sample_ids[i], "SNN_k10_k7"] <- adjustedRandIndex(spe_sub$layer_guess_reordered, spe_sub$SNN_k50_k7)
    ari.df[sample_ids[i], "batch_corr_SNN_k10_k7"] <- adjustedRandIndex(spe_sub$layer_guess_reordered, spe_sub$batch_corr_SNN_k10_k7)
    ari.df[sample_ids[i], "bayesSpace_pc"] <- adjustedRandIndex(spe_sub$layer_guess_reordered, spe_sub$spatial.cluster.pcs)
    ari.df[sample_ids[i], "bayesSpace"] <- adjustedRandIndex(spe_sub$layer_guess_reordered, spe_sub$spatial.cluster)
}

ari.df.long <- gather(ari.df, method, ari, SNN_k10_k7:bayesSpace, factor_key = TRUE)
save(ari.df.long, file = here::here("processed-data", "rdata", "pilot_dlpfc_data", "05_ARI", "pilot_ari_clustering_across.Rdata"))


#### make plots for fiure 2a, clustering across pilot data
# pdf(here::here("plots","05_ARI","pilot_data_ARI_clustering_across_2.pdf"))
# ggplot(ari.df.long, aes(x = method, y=ari)) +
#   geom_boxplot()+
#   theme_bw()+
#   geom_jitter(color="black", size=0.4, alpha=0.9)+
#   ylim(0,0.6)+
#   theme(text = element_text(size = 40))
# dev.off()



#### spaGCN ARI for pilot data

# load old ARI results
# load(file = here::here("processed-data", "rdata", "pilot_dlpfc_data","05_ARI", "pilot_ari_clustering_across.Rdata"))
# ari.df.long <-ari.df.long[1:48,]

# load pilot data
# load(file = here::here("processed-data","rdata","pilot_dlpfc_data","spe_pilot_bayesSpace_batch_corr_sampleID.Rdata"))

# sample_ids <- unique(spe$sample_id)
# df <- data.frame()
# for(i in seq_along(sample_ids)){
#   x <-read.csv(file = here::here("..","spython","spagcn","processed-data","03-our_data_analysis",sample_ids[i],"7_clusters","clusters.csv"))
#   df <- rbind(df,x)
# }
#
# write.csv(df,file = here::here("processed-data","rdata","pilot_dlpfc_data","clustering_results","spaGCN_k7","clusters.csv"))

# import spaGCN clusters
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "rdata", "pilot_dlpfc_data", "clustering_results"),
    prefix = "spaGCN_"
)

ari.df <- data.frame(matrix(ncol = 2, nrow = 12))
row.names(ari.df) <- sample_ids
colnames(ari.df) <- c("sample_id", "spaGCN_refined_cluster")

for (i in seq_along(sample_ids)) {
    spe_sub <- spe[, colData(spe)$sample_id == sample_ids[i]]
    ari.df$sample_id <- sample_ids[i]
    ari.df[sample_ids[i], "spaGCN_refined_cluster"] <- adjustedRandIndex(spe_sub$layer_guess_reordered, spe_sub$spaGCN_refined_cluster)
}
# append to ari.df.long
ari.df.long.2 <- gather(ari.df, method, ari, spaGCN_refined_cluster, factor_key = TRUE)

dim(ari.df.long)
# 48  3

dim(ari.df.long.2)
# 12  3

ari.df.long <- rbind(ari.df.long, ari.df.long.2)

dim(ari.df.long)
# 60  3

levels(ari.df.long$method) <- c(levels(ari.df.long$method), "Graph-Based", "Graph-based(BC)", "BayesSpace", "BayesSpace(BC)", "SpaGCN")
ari.df.long$method[ari.df.long$method == "SNN_k10_k7"] <- "Graph-Based"
ari.df.long$method[ari.df.long$method == "batch_corr_SNN_k10_k7"] <- "Graph-based(BC)"
ari.df.long$method[ari.df.long$method == "bayesSpace_pc"] <- "BayesSpace"
ari.df.long$method[ari.df.long$method == "bayesSpace"] <- "BayesSpace(BC)"
ari.df.long$method[ari.df.long$method == "spaGCN_refined_cluster"] <- "SpaGCN"

ari.df.long$general_method[c(1:12)] <- "GB"
ari.df.long$general_method[c(25:48)] <- "BS"
ari.df.long$general_method[c(49:60)] <- "SG"

save(ari.df.long, file = here::here("processed-data", "rdata", "pilot_dlpfc_data", "05_ARI", "pilot_ari_clustering_across.Rdata"))

# level_order <- c("Graph-Based","Graph-based(BC)","BayesSpace","BayesSpace(BC)","SpaGCN")
# pdf(here::here("plots","05_ARI","pilot_data_ARI_clustering_across.pdf"))
# ggplot(ari.df.long, aes(x = factor(method,level =  level_order), y=ari)) +
#   geom_boxplot(outlier.shape = NA)+
#   theme_bw()+
#   geom_jitter(color="black", size=1.0, alpha=0.9)+
#   ylim(0,0.6)+
#   theme(axis.text.x = element_text(size = 20,angle = 90, vjust = 0.5, hjust = 1, colour = c("blue","blue","red","red","green")),text = element_text(size = 30),axis.title = element_text(size = 30))+
#   ylab("Adjusted Rand Index")+
#   xlab("Clustering Method")
# dev.off()
