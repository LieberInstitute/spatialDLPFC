library("spatialLIBD")
library("SpatialExperiment")
library("here")
library("ggpubr")
library("mclust")
library("sessioninfo")

processed_dir <- here("spagcn", "processed-data", "03-our_data_analysis")
plot_dir <- here("spagcn", "plots", "03-our_data_analysis")
n_clusters <- 7

spe <- fetch_data("spe")
ari_list <- list()
spe_list <- list()

#   Compute ARI using 'clusters.csv' files for all samples
for (id in unique(spe$sample_id)) {
    #   Subset object to this sample only
    sub_spe <- spe[, colData(spe)$sample_id == id]

    #   Import 'clusters.csv'
    cluster_dir <- file.path(processed_dir, id, "7_clusters")
    sub_spe <- cluster_import(sub_spe, cluster_dir = cluster_dir)

    #   Compute ARI (raw and refined clusters vs. manual label)
    ari_list[[id]] <- c(
        adjustedRandIndex(
            sub_spe$imported_raw_cluster, sub_spe$layer_guess_reordered_short
        ),
        adjustedRandIndex(
            sub_spe$imported_refined_cluster,
            sub_spe$layer_guess_reordered_short
        )
    )

    spe_list[[id]] <- sub_spe
}

#   Form a data frame of ARI info, which will be useful for plotting
ari_raw <- as.numeric(sapply(ari_list, function(x) x[1]))
ari_refined <- as.numeric(sapply(ari_list, function(x) x[2]))

ari_df <- data.frame(
    "sample_id" = names(ari_list),
    "ARI" = c(ari_raw, ari_refined),
    "method" = rep(c("SpaGCN_raw", "SpaGCN_refined"), each = length(ari_list))
)

#   Boxplots showing ARI for each sample and method (raw or refined)
p <- ggboxplot(
    ari_df,
    x = "method", y = "ARI", add = "jitter", label = "sample_id",
    color = "method", palette = "Dark2", repel = TRUE,
    font.label = list(size = 10), legend = "none",
    ggtheme = theme_pubr(base_size = 20)
)

pdf(file.path(plot_dir, "ARI_boxplots.pdf"))
print(p)
dev.off()

#   We formed a list of SpatialExperiment objects, one object per sample, each
#   of which had cluster info imported from 'cluster_import'. We want this as
#   a single object to use with 'vis_grid_clus'
spe <- do.call(cbind, spe_list)

cols <- Polychrome::palette36.colors(n_clusters)
names(cols) <- sort(unique(spe$imported_raw_cluster))

#   Plot a grid of spatial cluster info for all samples, once for each cluster
#   type (raw and refined)
for (method in c("raw", "refined")) {
    vis_grid_clus(
        spe = spe,
        clustervar = paste0("imported_", method, "_cluster"),
        pdf_file = file.path(
            plot_dir, paste0(method, "_cluster_sample_grid.pdf")
        ),
        sort_clust = FALSE,
        colors = cols,
        spatial = FALSE,
        point_size = 2
    )
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
