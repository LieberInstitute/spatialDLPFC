library('spatialLIBD')
library('SpatialExperiment')
library('here')
library('ggpubr')
library('mclust')
library('sessioninfo')

processed_dir = here('spagcn', 'processed-data', '03-our_data_analysis')

spe = fetch_data("spe")
ari_list = list()

for(id in unique(spe$sample_id)) {
    #   Subset object to this sample only
    sub_spe = spe[,colData(spe)$sample_id == id]
    
    #   Import 'clusters.csv'
    cluster_dir = file.path(processed_dir, id, '7_clusters')
    sub_spe = cluster_import(sub_spe, cluster_dir = cluster_dir)
    
    #   Compute ARI
    ari_list[[id]] = adjustedRandIndex(
        sub_spe$Cluster,sub_spe$layer_guess_reordered_short
    )
}

ari_df = data.frame(
    'sample_id' = names(ari_list),
    'ARI' = as.numeric(ari_list),
    'method' = 'SpaGCN_raw_cluster'
)

p = ggboxplot(
    ari_df, x = 'method', y = 'ARI', add = "jitter", label = "sample_id",
    repel = TRUE, font.label = list(size = 10), legend = "none",
    ggtheme = theme_pubr(base_size = 20)
)

print(p)
