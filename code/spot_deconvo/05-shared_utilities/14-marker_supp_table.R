suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))

#   Number of marker genes to use per cell type
n_markers_per_type <- 25

#  Paths
sce_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    "sce_{cell_group}.rds"
)
marker_object_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    "marker_stats_{cell_group}.rds"
)

cell_types_broad <- c(
    "Astro", "EndoMural", "Micro", "Oligo", "OPC", "Excit", "Inhib"
)
cell_types_layer <- c(
    "Astro", "EndoMural", "Micro", "Oligo", "OPC", "Excit_L2/3", "Excit_L3",
    "Excit_L3/4/5", "Excit_L4", "Excit_L5", "Excit_L5/6", "Excit_L6",
    "Inhib"
)

out_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    "marker_stats_supp_table.csv"
)

marker_stats_list <- list()

for (cell_group in c("broad", "layer")) {
    sce <- readRDS(sub("\\{cell_group\\}", cell_group, sce_in))
    marker_stats <- readRDS(
        sub("\\{cell_group\\}", cell_group, marker_object_in)
    )

    if (cell_group == "broad") {
        cell_types <- cell_types_broad
    } else {
        cell_types <- cell_types_layer
    }

    #---------------------------------------------------------------------------
    #   Filter out mitochondrial genes and re-rank 'rank_ratio' values. Add gene
    #   symbols to 'marker_stats' object
    #---------------------------------------------------------------------------

    #   Add gene symbol
    marker_stats$symbol <- rowData(sce)$gene_name[
        match(marker_stats$gene, rownames(sce))
    ]

    #   Filter out mitochondrial genes
    marker_stats <- marker_stats[!grepl("^MT-", marker_stats$symbol), ]

    stopifnot(
        identical(
            sort(as.character(unique(marker_stats$cellType.target))),
            sort(cell_types)
        )
    )

    #   "Re-rank" rank_ratio, since there may be missing ranks now
    for (ct in cell_types) {
        old_ranks <- marker_stats |>
            filter(cellType.target == ct) |>
            pull(rank_ratio) |>
            sort()

        for (i in 1:length(which((marker_stats$cellType.target == ct)))) {
            index <- which(
                (marker_stats$cellType.target == ct) &
                    (marker_stats$rank_ratio == old_ranks[i])
            )
            stopifnot(length(index) == 1)
            marker_stats[index, "rank_ratio"] <- i
        }
    }

    marker_stats <- marker_stats |>
        #   Take top N marker genes for each cell type
        filter(
            rank_ratio <= n_markers_per_type,
            ratio > 1
        ) |>
        #   Label with cell-type resolution
        mutate(cellTypeResolution = cell_group)

    #   Warn if less than the intended number of markers is used for any cell
    #   type
    num_markers_table <- marker_stats |>
        group_by(cellType.target) |>
        summarize(num_markers = n())

    if (any(num_markers_table$num_markers < n_markers_per_type)) {
        warning(
            paste(
                "Used less than", n_markers_per_type,
                "markers for at least one cell type."
            )
        )
        print("Number of markers per cell type:")
        print(num_markers_table)
    }

    stopifnot(all(num_markers_table$num_markers > 0))

    marker_stats_list[[cell_group]] <- marker_stats
}

#   Combine tables for both cell-type resolutions and write to CSV
do.call(rbind, marker_stats_list) |>
    mutate(cellTypeResolution = as.factor(cellTypeResolution)) |>
    write.csv(file = out_path, quote = FALSE, row.names = FALSE)

session_info()
