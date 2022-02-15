library('spatialLIBD')
library('SpatialExperiment')
library('here')
library('ggpubr')
library('mclust')
library('sessioninfo')

processed_dir = here('spagcn', 'processed-data', '03-our_data_analysis')

spe = fetch_data("spe")
ari_list = list()

#   Compute ARI using 'clusters.csv' files for all samples
for(id in unique(spe$sample_id)) {
    #   Subset object to this sample only
    sub_spe = spe[,colData(spe)$sample_id == id]
    
    #   Import 'clusters.csv'
    cluster_dir = file.path(processed_dir, id, '7_clusters')
    sub_spe = cluster_import(sub_spe, cluster_dir = cluster_dir)
    
    #   Compute ARI (raw and refined clusters vs. manual label)
    ari_list[[id]] = c(
        adjustedRandIndex(
            sub_spe$imported_raw_cluster, sub_spe$layer_guess_reordered_short
        ),
        adjustedRandIndex(
            sub_spe$imported_refined_cluster,
            sub_spe$layer_guess_reordered_short
        )
    )
}

#   Form a data frame of ARI info, which will be useful for plotting
ari_raw = as.numeric(sapply(ari_list, function(x) x[1]))
ari_refined = as.numeric(sapply(ari_list, function(x) x[2]))

ari_df = data.frame(
    'sample_id' = names(ari_list),
    'ARI' = c(ari_raw, ari_refined),
    'method' = rep(c('SpaGCN_raw', 'SpaGCN_refined'), each = length(ari_list))
)

#   Boxplots showing ARI for each sample and method (raw or refined)
p = ggboxplot(
    ari_df, x = 'method', y = 'ARI', add = "jitter", label = "sample_id",
    color = 'method', palette = 'Dark2', repel = TRUE,
    font.label = list(size = 10), legend = "none",
    ggtheme = theme_pubr(base_size = 20)
)

pdf(file.path(processed_dir, 'ARI_boxplots.pdf'))
print(p)
dev.off()
