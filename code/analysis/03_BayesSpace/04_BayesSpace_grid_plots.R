library("here")
library("sessioninfo")
library("spatialLIBD")
library("Polychrome")
library("ggplot2")

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# load spe object and clusters
load(
    file = here::here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters.Rdata"
    ),
    verbose = TRUE
)

mycolors <- Polychrome::palette36.colors(k)
## Reverse colors only for k2 to avoid WM vs GM confusion
if (k == 2) {
    mycolors <- rev(mycolors)
}
names(mycolors) <-
    sort(unique(colData(spe)[[paste0("bayesSpace_harmony_", k)]]))
sample_order <- unique(spe$sample_id)


pk <- vis_grid_clus(
    spe = spe,
    clustervar = paste0("bayesSpace_harmony_", k),
    sort_clust = FALSE,
    colors = mycolors,
    spatial = FALSE,
    point_size = 2,
    return_plots = TRUE
)

pk <- lapply(sample_order, function(sampleid) {
    p <- pk[[sampleid]]
    p + theme(legend.position = "none",
        plot.title = ggplot2::element_text(size = 40))
})
names(pk) <- sample_order

pdf(
    file = here::here(
        "plots",
        "03_BayesSpace",
        paste0("polychrome_vis_grid_clus_sfigu_BayesSpace_k", k, ".pdf")
    ),
    height = 5 * 8,
    width = 6 * 8
)
print(cowplot::plot_grid(plotlist = pk))
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
