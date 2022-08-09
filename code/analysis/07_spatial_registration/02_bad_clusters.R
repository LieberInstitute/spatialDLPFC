library(SpatialExperiment)
library(here)
library(pheatmap)
library(spatialLIBD)

# #load spe object
# load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)
#
# #load clusters
# spe <- cluster_import(
#   spe,
#   cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"),
#   prefix = ""
# )

# loaded this data object to make percent plots
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")


table_percent <- function(input.table) {
    list(
        "counts" = addmargins(input.table),
        "percent_all" = addmargins(round(input.table / sum(input.table) * 100, 2)), # divide my all spots

        "percent_col" = addmargins(round(sweep(input.table, 2, colSums(input.table), "/") * 100, 2)), # divide by spots in bayesSpace cluster

        "percent_row" = addmargins(round(sweep(input.table, 1, rowSums(input.table), "/") * 100, 2)) # divide by spots in Kristen's annotation
    )
}

# make table of bayesSpace k = 9 vs k = 28
my.table.28 <- with(colData(spe), table(bayesSpace_harmony_9, bayesSpace_harmony_28))
my.table.28[my.table.28 == 0] <- NA

pdf(file = here::here("plots", "07a_bad_clusters", "clusters_k9_k28.pdf"))
pheatmap(
    my.table.28,
    show_rownames = TRUE,
    cluster_rows = FALSE,
    cluster_cols = TRUE
)
dev.off()

# plot log transformed
my.table.log <- with(colData(spe), log10(table(bayesSpace_harmony_9, bayesSpace_harmony_28) + 1))

pdf(file = here::here("plots", "07a_bad_clusters", "clusters_k9_k28_log.pdf"))
pheatmap(
    my.table.log,
    show_rownames = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE
)
dev.off()

# plot it with percent
my.table.28 <- with(colData(spe), table(bayesSpace_harmony_9, bayesSpace_harmony_28))
my.table.28.percent <- table_percent(my.table.28)$percent_row
my.table.28.percent[my.table.28.percent == 0] <- NA

pdf(file = here::here("plots", "07a_bad_clusters", "clusters_k9_k28_percent.pdf"))
pheatmap(
    my.table.28.percent[c(1:9), c(1:26)],
    show_rownames = TRUE,
    cluster_rows = FALSE,
    cluster_cols = TRUE
)
dev.off()



# make table of bayesSpace k = 9 vs k = 16
my.table.16 <- with(colData(spe), table(bayesSpace_harmony_9, bayesSpace_harmony_16))
my.table.16[my.table.16 == 0] <- NA

pdf(file = here::here("plots", "07a_bad_clusters", "clusters_k9_k16.pdf"))
pheatmap(
    my.table.16,
    show_rownames = TRUE,
    cluster_rows = FALSE,
    cluster_cols = TRUE
)
dev.off()

# plot it with percent
my.table.16 <- with(colData(spe), table(bayesSpace_harmony_9, bayesSpace_harmony_16))
my.table.16.percent <- table_percent(my.table.16)$percent_row
my.table.16.percent[my.table.16.percent == 0] <- NA

pdf(file = here::here("plots", "07a_bad_clusters", "clusters_k9_k16_percent.pdf"))
pheatmap(
    my.table.16.percent[c(1:9), c(1:16)],
    show_rownames = TRUE,
    cluster_rows = FALSE,
    cluster_cols = TRUE
)
dev.off()


# make orange and gray plots of cluster 1 at k = 9
spe$BayesSpace_k9_c1 <- factor(spe$bayesSpace_harmony_9 == 1,
    levels = c("TRUE", "FALSE")
)
vis_grid_clus(
    spe = spe,
    clustervar = "BayesSpace_k9_c1",
    pdf = here("plots", "07a_bad_clusters", paste0("vis_clus_BayesSpace_k9_c1.pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange"),
    point_size = 0.75,
    height = 24,
    width = 90
)

# make orange and gray plots of cluster 2 at k = 9
spe$BayesSpace_k9_c2 <- factor(spe$bayesSpace_harmony_9 == 2,
    levels = c("TRUE", "FALSE")
)
vis_grid_clus(
    spe = spe,
    clustervar = "BayesSpace_k9_c2",
    pdf = here("plots", "07a_bad_clusters", paste0("vis_clus_BayesSpace_k9_c2.pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange"),
    point_size = 0.75,
    height = 24,
    width = 90
)

# make orange and gray plots of cluster 23 at k = 28
spe$BayesSpace_k28_c23 <- factor(spe$bayesSpace_harmony_28 == 23,
    levels = c("TRUE", "FALSE")
)
vis_grid_clus(
    spe = spe,
    clustervar = "BayesSpace_k28_c23",
    pdf = here("plots", "07a_bad_clusters", paste0("vis_clus_BayesSpace_k28_c23.pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange"),
    point_size = 0.75,
    height = 24,
    width = 90
)

# make orange and gray plots of cluster 27 at k = 28
spe$BayesSpace_k28_c27 <- factor(spe$bayesSpace_harmony_28 == 27,
    levels = c("TRUE", "FALSE")
)
vis_grid_clus(
    spe = spe,
    clustervar = "BayesSpace_k28_c27",
    pdf = here("plots", "07a_bad_clusters", paste0("vis_clus_BayesSpace_k28_c27.pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange"),
    point_size = 0.75,
    height = 24,
    width = 90
)
