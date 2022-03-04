library("SpatialExperiment")
library("mclust")
library("spatialLIBD")
library(tidyr)
library(ggplot2)

load(file = here::here("processed-data","rdata","pilot_dlpfc_data","spe_pilot_bayesSpace_batch_corr_sampleID.Rdata"))

##import clusters from bayesSpace, graph-based, spaGCN. Might have to remake pilot spe objects and export clusters 

#run ARI and create plot

#plot ARI for pilot data comparing different clustering algorithms to Kristen's manual annnotations
Method <- c("BayesSpace", "BayesSpace Batch Corrected", "Graph-Based", "Graph-Based Batch Corrected")
ARI<- c(0.21, 0.34, 0.13, 0.16)

df <- data.frame(Method, ARI)

print (df)

pdf(here::here("plots","pilot_ARI.pdf"))
ggplot(df, aes(x=Method, y=ARI)) +
  geom_point(size=2) +
  theme_bw()
dev.off()

#plot ARI for my data comparing different clustering methods against BayesSpace batch corrected

spe<- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","graph_based_harmony"),
  prefix = "harmony_graph_based_"
)

spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","semi_supervised_harmony_across_samples"),
  prefix = "harmony_semi_supervised_"
)

spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","graph_based_pcs"),
  prefix = "pcs_graph_based_"
)


adjustedRandIndex(spe$spatial.cluster,spe$harmony_semi_supervised_SNN_k10_k7)
# 0.1958784

adjustedRandIndex(spe$spatial.cluster,spe$harmony_graph_based_SNN_k10_k7)
# 0.1958784

adjustedRandIndex(spe$spatial.cluster,spe$pcs_graph_based_SNN_k10_k7)
# 0.213768

Method <- c("Graph-Based PCs", "Graph-Based Batch Corrected", "Semi-supervised Batch Corrected")
ARI<- c(0.213768, 0.1958784, 0.1958784)

df <- data.frame(Method, ARI)

print (df)

pdf(here::here("plots","my_data_ARI.pdf"))
ggplot(df, aes(x=Method, y=ARI)) +
  geom_point(size=2) +
  ylim(0,0.6)+
  theme_bw()
dev.off()


#### spaGCN ARI for pilot data
spe.temp <- cluster_import(
  spe,
  cluster_dir = here::here("..","spython","spagcn","processed-data","03-our_data_analysis"),
  prefix = "spaGCN_"
)


with(colData(spe),addmargins(table(spatial.cluster,pseudobulk_PCA.y,sample_id)))

### ari for semi-supervised within sample clustering 
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/clustering_results/semi_supervised_pcs_within_samples/d_plot_semi_supervised.Rdata")
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/spe_final.Rdata")


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

spe$semi_supervised_pcs_within_samples <- spe$pseudobulk_PCA 
spe$semi_supervised_UMAP_within_samples <- spe$pseudobulk_UMAP 

cluster_export(
  spe,
  "semi_supervised_pcs_within_samples",
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

cluster_export(
  spe,
  "semi_supervised_UMAP_within_samples",
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","graph_based_within_samples"),
  prefix = ""
)

sample_ids <- unique(spe$sample_id)
ari.df <- data.frame(matrix(ncol = 4, nrow = 30))
row.names(ari.df) <- sample_ids
colnames(ari.df)<-c("sample_id","ari.semi.supervised.pca","ari.semi.supervised.umap","ari.graph.based.pca")

for (i in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[i]]
  ari.df$sample_id <-sample_ids[i]
  ari.df[sample_ids[i],"ari.semi.supervised.pca"]<-adjustedRandIndex(spe_sub$spatial.cluster,spe_sub$semi_supervised_within_PCA)
  ari.df[sample_ids[i],"ari.semi.supervised.umap"]<-adjustedRandIndex(spe_sub$spatial.cluster,spe_sub$semi_supervised_within_UMAP)
  ari.df[sample_ids[i],"ari.graph.based.pca"]<-adjustedRandIndex(spe_sub$spatial.cluster,spe_sub$graph_based_PCA)
  
}

# The arguments to gather():
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)
ari.df.long <- gather(ari.df, method, ari, ari.semi.supervised.pca:ari.graph.based.pca, factor_key=TRUE)

pdf(here::here("plots","my_data_ARI_semi_supervised.pdf"))
ggplot(ari.df.long, aes(x = method, y=ari)) + 
  geom_boxplot()+
  theme_bw()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  ylim(0,0.6)
dev.off()

save(ari.df.long, file = here::here("processed-data", "rdata", "spe", "ari_semi_supervised_within.Rdata" ))

#### redo ari calculations for across sample clustering
spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","graph_based_harmony"),
  prefix = "graph_based_corrected_across"
)

spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","graph_based_pcs"),
  prefix = "graph_based_pca_across"
)

spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","semi_supervised_harmony_across_samples"),
  prefix = "semi_supervised_corrected_acrosss"
)

sample_ids <- unique(spe$sample_id)
ari.df <- data.frame(matrix(ncol = 4, nrow = 30))
row.names(ari.df) <- sample_ids
colnames(ari.df)<-c("sample_id","graph_based_pca_acrossSNN_k10_k7","semi_supervised_corrected_acrosssSNN_k10_k7","graph_based_corrected_acrossSNN_k10_k7.x")

for (i in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[i]]
  ari.df$sample_id <-sample_ids[i]
  ari.df[sample_ids[i],"graph_based_pca_acrossSNN_k10_k7"]<-adjustedRandIndex(spe_sub$spatial.cluster,spe_sub$graph_based_pca_acrossSNN_k10_k7)
  ari.df[sample_ids[i],"semi_supervised_corrected_acrosssSNN_k10_k7"]<-adjustedRandIndex(spe_sub$spatial.cluster,spe_sub$semi_supervised_corrected_acrosssSNN_k10_k7)
  ari.df[sample_ids[i],"graph_based_corrected_acrossSNN_k10_k7.x"]<-adjustedRandIndex(spe_sub$spatial.cluster,spe_sub$graph_based_corrected_acrossSNN_k10_k7.x)
  
}

ari.df.long <- gather(ari.df, method, ari, graph_based_pca_acrossSNN_k10_k7:graph_based_corrected_acrossSNN_k10_k7.x, factor_key=TRUE)

pdf(here::here("plots","my_data_ARI_clustering_across.pdf"))
ggplot(ari.df.long[31:90,], aes(x = method, y=ari)) + 
  geom_boxplot()+
  theme_bw()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  ylim(0,0.6)
dev.off()

save(ari.df.long, file = here::here("processed-data", "rdata", "spe", "ari_clustering_across.Rdata"))

#plot ARI for pilot data comparing different clustering algorithms to Kristen's manual annnotations
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/pilot_dlpfc_data/spe_pilot_bayesSpace_pcs.Rdata")
spe.bspc <- spe
spe$spatial.cluster.pcs <- spe.bspc$spatial.cluster
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/pilot_dlpfc_data/spe_pilot_bayesSpace_batch_corr_sampleID.Rdata")


sample_ids <- unique(spe$sample_id)
ari.df <- data.frame(matrix(ncol = 5, nrow = 12))
row.names(ari.df) <- sample_ids
colnames(ari.df)<-c("sample_id","SNN_k10_k7","batch_corr_SNN_k10_k7","bayesSpace_pc","bayesSpace")

for (i in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[i]]
  ari.df$sample_id <-sample_ids[i]
  ari.df[sample_ids[i],"SNN_k10_k7"]<-adjustedRandIndex(spe_sub$layer_guess_reordered,spe_sub$SNN_k50_k7)
  ari.df[sample_ids[i],"batch_corr_SNN_k10_k7"]<-adjustedRandIndex(spe_sub$layer_guess_reordered,spe_sub$batch_corr_SNN_k10_k7)
  ari.df[sample_ids[i],"bayesSpace_pc"]<-adjustedRandIndex(spe_sub$layer_guess_reordered,spe_sub$spatial.cluster.pcs)
  ari.df[sample_ids[i],"bayesSpace"]<-adjustedRandIndex(spe_sub$layer_guess_reordered,spe_sub$spatial.cluster)
  
}

ari.df.long <- gather(ari.df, method, ari, SNN_k10_k7:bayesSpace, factor_key=TRUE)

pdf(here::here("plots","pilot_ARI_clustering_across.pdf"))
ggplot(ari.df.long, aes(x = method, y=ari)) + 
  geom_boxplot()+
  theme_bw()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  ylim(0,0.6)
dev.off()

save(ari.df.long, file = here::here("processed-data", "rdata", "spe", "pilot_ari_clustering_across.Rdata"))

####remaking plots for fiure 2b, clustering across
pdf(here::here("plots","my_data_ARI_clustering_across_2.pdf"))
ggplot(ari.df.long[31:90,], aes(x = method, y=ari)) + 
  geom_boxplot()+
  theme_bw()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  ylim(0,0.6)+
  theme(text = element_text(size = 10)) 
dev.off()

####remaking plots for fiure 2b, clustering within
pdf(here::here("plots","my_data_ARI_clustering_within_2.pdf"))
ggplot(ari.df.long, aes(x = method, y=ari)) + 
  geom_boxplot()+
  theme_bw()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  ylim(0,0.6)+
  theme(text = element_text(size = 40)) 
dev.off()


####remaking plots for fiure 2a, clustering across pilot data
pdf(here::here("plots","pilot_data_ARI_clustering_across_2.pdf"))
ggplot(ari.df.long, aes(x = method, y=ari)) + 
  geom_boxplot()+
  theme_bw()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  ylim(0,0.6)+
  theme(text = element_text(size = 40)) 
dev.off()