#   Find marker genes that will be shared for both the IF and non-IF analyses.

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))

#   Number of marker genes to use per cell type. Note that cell2location seems
#   to need more markers than the other tools, motivating the exception below
n_markers_per_type <- 25
n_markers_per_type_c2l <- 100
stopifnot(n_markers_per_type_c2l >= n_markers_per_type)

#  Paths
sce_in <- here(
    "processed-data", "spot_deconvo", "sce.rds"
)
spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)
spe_nonIF_in <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
)
marker_object_in <- here(
    "processed-data", "spot_deconvo", "marker_stats.rds"
)

marker_out <- here(
    "processed-data", "spot_deconvo", "markers.txt"
)
marker_c2l_out <- here(
    "processed-data", "spot_deconvo", "markers_C2L.txt"
)

plot_dir <- here(
    "plots", "spot_deconvo", "05-shared_utilities"
)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
#  Load objects
###############################################################################

sce = readRDS(sce_in)
spe_IF <- readRDS(spe_IF_in)
load(spe_nonIF_in, verbose = TRUE)
spe_nonIF <- spe
rm(spe)

marker_stats <- readRDS(marker_object_in)

gc()

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
    
    #   All the marker genes are present in the single-cell data (sanity check),
    #   but some markers might be in neither the IF nor non-IF spatial data
    stopifnot(all(markers_scratch %in% rowData(sce)$gene_id))
    
    start_len <- length(markers_scratch)
    markers_scratch <- markers_scratch[
        (markers_scratch %in% rowData(spe_IF)$gene_id) &
            (markers_scratch %in% rowData(spe_nonIF)$gene_id)
    ]
    print(
        paste(
            "Dropped", start_len - length(markers_scratch),
            "markers, which were not in the spatial data."
        )
    )
    
    #   Write list of markers
    writeLines(markers_scratch, con = out_path)
}

plot_markers = function(sce, marker_stats, n_genes, plot_name, ct) {
    p <- plot_marker_express(
        sce,
        marker_stats,
        cell_type = ct,
        n_genes = n_genes,
        rank_col = "rank_ratio",
        anno_col = "anno_ratio",
        cellType_col = "cellType_broad_hc"
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

#   Visually show how markers look for each cell type. In particular, look at
#   the best 5 markers, the worst 5 markers used for C2L, and the worst 5
#   markers used for tangram/SPOTlight
walk(
    unique(marker_stats$cellType.target),
    function(ct){
        #   First plot the top 5 markers
        plot_markers(sce, marker_stats, 5, "marker_genes_1-5", ct)
        
        #   Next, plot worst 5 C2L markers: those with ranks in
        #   [n_markers_per_type_c2l - 4, n_markers_per_type_c2l]
        rank_range = paste(
            as.character(n_markers_per_type_c2l - 4),
            as.character(n_markers_per_type_c2l),
            sep = "-"
        )
        marker_stats_temp = marker_stats %>%
            filter(rank_ratio <= n_markers_per_type_c2l) %>%
            mutate(
                rank_ratio, 
                rank_ratio = 1 + n_markers_per_type_c2l - rank_ratio
            )
        plot_markers(
            sce, marker_stats_temp, 5,
            paste0("marker_genes_", rank_range), ct
        )
        
        #   Finally, plot worst 5 non-C2L markers
        rank_range = paste(
            as.character(n_markers_per_type - 4),
            as.character(n_markers_per_type),
            sep = "-"
        )
        marker_stats_temp = marker_stats %>%
            filter(rank_ratio <= n_markers_per_type) %>%
            mutate(
                rank_ratio, 
                rank_ratio = 1 + n_markers_per_type - rank_ratio
            )
        plot_markers(
            sce, marker_stats_temp, 5,
            paste0("marker_genes_", rank_range), ct
        )
    }
)

#   Plot mean ratio against log fold-change for all genes, split by target cell
#   type and colored by whether each gene will be used as a marker
p = marker_stats %>% 
    mutate(
        Marker = case_when(
            rank_ratio <= n_markers_per_type ~ paste0('Marker top', n_markers_per_type),
            rank_ratio <= n_markers_per_type_c2l ~ paste0('Marker top', n_markers_per_type_c2l),
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
boxplot_mean_ratio(n_markers_per_type_c2l, "mean_ratio_boxplot_C2L")

session_info()
