#   Find marker genes that will be shared for both the IF and non-IF analyses.

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("HDF5Array"))

cell_group = "layer" # "broad" or "layer"

#   Number of marker genes to use per cell type. Note that cell2location seems
#   to need more markers than the other tools, motivating the exception below
n_markers_per_type <- 25
n_markers_per_type_c2l <- 25 # 100 for "broad"
stopifnot(n_markers_per_type_c2l >= n_markers_per_type)

#  Paths
sce_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("sce_", cell_group, ".rds")
)
marker_object_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("marker_stats_", cell_group, ".rds")
)

marker_out <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("markers_", cell_group, ".txt")
)
marker_c2l_out <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("markers_C2L_", cell_group, ".txt")
)

plot_dir <- here(
    "plots", "spot_deconvo", "05-shared_utilities", cell_group
)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
#  Load objects
###############################################################################

sce = readRDS(sce_in)
marker_stats <- readRDS(marker_object_in)
gc()

print(paste0("Running script at ", cell_group, "-resolution."))

###############################################################################
#  Functions
###############################################################################

write_markers = function(n_markers, out_path) {
    #   Take top N marker genes for each cell type
    marker_stats_temp <- marker_stats %>%
        filter(rank_ratio <= n_markers)
    
    markers_scratch <- marker_stats_temp$gene
    
    #   Verify all markers have more expression for the target cell type than
    #   next highest
    stopifnot(any(marker_stats_temp$ratio > 1))
    
    #   Write list of markers
    writeLines(markers_scratch, con = out_path)
}

plot_markers = function(
    sce, marker_stats, n_genes, ct, cell_column
) {
    p <- plot_marker_express(
        sce,
        marker_stats,
        cell_type = ct,
        n_genes = n_genes,
        rank_col = "rank_ratio",
        anno_col = "anno_ratio",
        cellType_col = cell_column
    )
    
    return(p)
}

#   Plot mean-ratio distribution by cell type/ layer
boxplot_mean_ratio = function(n_markers, plot_name) {
    p = marker_stats %>%
        filter(rank_ratio <= n_markers) %>%
        mutate(ratio, ratio = log(ratio)) %>%
        ggplot(aes(cellType.target, ratio)) +
        geom_boxplot() +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
        labs(y = "log(Mean Ratio)") +
        theme_bw()
    
    ggsave(
        p, filename = file.path(plot_dir, paste0(plot_name, ".png")),
        height = 10, width = 10
    )
}

###############################################################################
#  Subset and write markers
###############################################################################

print("Writing markers for tangram/SPOTlight...")
write_markers(n_markers_per_type, marker_out)

print("Writing markers for cell2location...")
write_markers(n_markers_per_type_c2l, marker_c2l_out)

###############################################################################
#  Visually check quality of markers
###############################################################################

if (cell_group == "broad") {
    cell_column = "cellType_broad_hc"
} else {
    cell_column = "layer_level"
}
#   Visually show how markers look for each cell type. In particular, look at
#   the best 5 markers, the worst 5 markers used for C2L, and the worst 5
#   markers used for tangram/SPOTlight
plot_list = lapply(
    unique(marker_stats$cellType.target),
    function(ct){
        plot_markers(
            sce, marker_stats, n_markers_per_type, ct, cell_column
        )
    }
)

#   Write a multi-page PDF with violin plots for each cell group and all
#   markers
# pdf(file.path(plot_dir, paste0("marker_gene_violin.pdf")))
png(
    file.path(plot_dir, paste0("marker_gene_violin.png")),
    width = 2000,
    height = 480 * length(unique(marker_stats$cellType.target))
)
lapply(plot_list, print)
dev.off()

#   Plot mean ratio against log fold-change for all genes, split by target cell
#   type and colored by whether each gene will be used as a marker
p = marker_stats %>% 
    mutate(
        Marker = case_when(
            rank_ratio <= n_markers_per_type ~ paste0('Marker top', n_markers_per_type),
            # rank_ratio <= n_markers_per_type_c2l ~ paste0('Marker top', n_markers_per_type_c2l),
            TRUE ~ 'Not marker'
        )
    ) %>%
    ggplot(aes(ratio, std.logFC, color = Marker)) +
    geom_point(size = 0.5, alpha = 0.5) +
    facet_wrap(~cellType.target, scales = "free_x") +
    labs(x = "Mean Ratio") +
    theme_bw()

ggsave(
    p, filename = file.path(plot_dir, paste0("mean_ratio_vs_1vall.png")),
    height = 10, width = 10
)

#   Plot mean-ratio distibution by group (cell type or layer label)
boxplot_mean_ratio(n_markers_per_type, "mean_ratio_boxplot")
# boxplot_mean_ratio(n_markers_per_type_c2l, "mean_ratio_boxplot_C2L")

session_info()
