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
    
    #   It's technically possible to find a single gene that is used as a
    #   "marker" for two different cell types via this method. Verify this is
    #   not the case, because that would indicate the use of bad markers
    stopifnot(
        length(unique(markers_scratch)) == n_markers * length(unique(marker_stats_temp$cellType.target))
    )
    
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
    ggsave(
        p,
        filename = file.path(
            plot_dir,
            paste0(plot_name, "_", ct, ".png")
        ),
        height = 10, width = 10
    )
}

#   Plot mean-ratio distribution by cell type/ layer
boxplot_mean_ratio = function(n_markers, plot_name) {
    p = marker_stats %>%
        filter(rank_ratio <= n_markers) %>%
        ggplot(aes(cellType.target, ratio)) +
        geom_boxplot() +
        labs(y = "Mean Ratio") +
        theme_bw()
    
    return(p)
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
plots_top_5 = list()
plots_bottom_5 = list()
walk(
    unique(marker_stats$cellType.target),
    function(ct){
        #  Plot the top 5 markers
        plots_top_5[[ct]] = plot_markers(
            sce, marker_stats, 5, ct, cell_column)
        
        #   Plot the worst 5 markers
        marker_stats_temp = marker_stats %>%
            filter(rank_ratio <= n_markers_per_type) %>%
            mutate(
                rank_ratio, 
                rank_ratio = 1 + n_markers_per_type - rank_ratio
            )
        plots_bottom_5[[ct]] = plot_markers(
            sce, marker_stats_temp, 5, ct, cell_column
        )
    }
)

#   Write a multi-page PDF with violin plots for each cell group and the top
#   5 markers
pdf(file.path(plot_dir, paste0("marker_genes_1-5.pdf")))
lapply(plots_top_5, print)
dev.off()

#   Write a multi-page PDF with violin plots for each cell group and the worst
#   5 markers
rank_range = paste(
    as.character(n_markers_per_type - 4),
    as.character(n_markers_per_type),
    sep = "-"
)
pdf(file.path(plot_dir, paste0("marker_genes_", rank_range, ".pdf")))
lapply(plots_bottom_5, print)
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
boxplot_mean_ratio(n_markers_per_type, "mean_ratio_boxplot_tangram")
# boxplot_mean_ratio(n_markers_per_type_c2l, "mean_ratio_boxplot_C2L")

session_info()
