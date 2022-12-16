#   Find marker genes that will be shared for both the IF and non-IF analyses.

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("HDF5Array"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("cowplot"))

#   Adds the 'spot_plot' function, a wrapper for 'vis_gene' or 'vis_clus' with
#   consistent manuscript-appropriate settings
source(
    here("code", "spot_deconvo", "05-shared_utilities", "shared_functions.R")
)

cell_group <- "layer" # "broad" or "layer"

#   Number of marker genes to use per cell type
n_markers_per_type <- 25

#  Paths
sce_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("sce_", cell_group, ".rds")
)
spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)
spe_nonIF_in <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
)
marker_object_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("marker_stats_", cell_group, ".rds")
)

marker_out <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("markers_", cell_group, ".txt")
)

plot_dir <- here(
    "plots", "spot_deconvo", "05-shared_utilities", cell_group
)

#   Symbols for markers of Layer 1-6, GM, and WM, respectively, for reference in
#   some plots
classical_markers <- c(
    "AQP4", "HPCAL1", "CUX2", "RORB", "PCP4", "KRT17", "SNAP25", "MOBP"
)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
#  Load objects and preprocess 'marker_stats'
###############################################################################

sce <- readRDS(sce_in)
marker_stats <- readRDS(marker_object_in)
gc()

print(paste0("Running script at ", cell_group, "-resolution."))

#-------------------------------------------------------------------------------
#   Filter out mitochondrial genes and re-rank 'rank_ratio' values. Add gene
#   symbols to 'marker_stats' object
#-------------------------------------------------------------------------------

#   Add gene symbol
marker_stats$symbol <- rowData(sce)$gene_name[
    match(marker_stats$gene, rownames(sce))
]

#   Filter out mitochondrial genes
marker_stats <- marker_stats[!grepl("^MT-", marker_stats$symbol), ]

#   "Re-rank" rank_ratio, since there may be missing ranks now
for (ct in unique(marker_stats$cellType.target)) {
    old_ranks <- marker_stats %>%
        filter(cellType.target == ct) %>%
        pull(rank_ratio) %>%
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

###############################################################################
#  Functions
###############################################################################

write_markers <- function(n_markers, out_path) {
    #   Take top N marker genes for each cell type
    marker_stats_temp <- marker_stats %>%
        filter(
            rank_ratio <= n_markers,
            ratio > 1
        )

    #   Warn if less than the intended number of markers is used for any cell
    #   type
    num_markers_table <- marker_stats_temp %>%
        group_by(cellType.target) %>%
        summarize(num_markers = n())

    if (any(num_markers_table$num_markers < n_markers)) {
        warning(
            paste(
                "Used less than", n_markers,
                "markers for at least one cell type."
            )
        )
        print("Number of markers per cell type:")
        print(num_markers_table)
    }

    stopifnot(all(num_markers_table$num_markers > 0))

    #   Write list of markers
    writeLines(marker_stats_temp$gene, con = out_path)
}

my_plotExpression <- function(sce, genes, assay = "logcounts", ct = "cellType", fill_colors = NULL,
    title = NULL, marker_stats) {
    cat_df <- as.data.frame(colData(sce))[, ct, drop = FALSE]
    expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes, ]))

    cat <- cat_df[expression_long$Var2, ]
    expression_long <- cbind(expression_long, cat)

    #   Use gene symbols for labels, not Ensembl ID
    symbols <- rowData(sce)$gene_name[match(genes, rownames(sce))]
    names(symbols) <- genes

    #   Add a data frame for adding mean-ratio labels to each gene
    text_df <- marker_stats
    text_df$ratio <- paste0("Mean ratio: ", round(text_df$ratio, 2))
    text_df$Var1 <- factor(text_df$gene, levels = levels(expression_long$Var1))

    expression_violin <- ggplot(
        data = expression_long, aes(x = cat, y = value, fill = cat)
    ) +
        geom_violin(scale = "width") +
        geom_text(
            data = text_df,
            mapping = aes(
                x = length(unique(sce[[ct]])), y = Inf, fill = NULL,
                label = ratio
            ),
            size = 10, hjust = 1, vjust = 1
        ) +
        facet_wrap(
            ~Var1,
            ncol = 5, scales = "free_y",
            labeller = labeller(Var1 = symbols)
        ) +
        labs(
            y = paste0("Expression (", assay, ")"),
            title = title
        ) +
        theme_bw(base_size = 35) +
        theme(
            legend.position = "None", axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(face = "italic")
        ) +
        stat_summary(
            fun = median,
            # fun.min = median,
            # fun.max = median,
            geom = "crossbar",
            width = 0.3
        )

    if (!is.null(fill_colors)) expression_violin <- expression_violin + scale_fill_manual(values = fill_colors)

    # expression_violin
    return(expression_violin)
}

if (cell_group == "broad") {
    colors_col <- "cell_type_colors_broad"
    cell_column <- "cellType_broad_hc"
    cell_type_nrow <- 2
} else {
    colors_col <- "cell_type_colors_layer"
    cell_column <- "layer_level"
    cell_type_nrow <- 3
}

#   Plot mean-ratio distribution by cell type/ layer
boxplot_mean_ratio <- function(n_markers, plot_name) {
    p <- marker_stats %>%
        filter(rank_ratio <= n_markers) %>%
        mutate(ratio, ratio = log(ratio)) %>%
        ggplot(aes(cellType.target, ratio, color = cellType.target)) +
        geom_boxplot() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        scale_color_manual(values = metadata(sce)[[colors_col]]) +
        labs(y = "log(Mean Ratio)") +
        theme_bw(base_size = 25) +
        guides(color = "none")

    ggsave(
        p,
        filename = file.path(plot_dir, paste0(plot_name, ".png")),
        height = 10, width = 10
    )
}

###############################################################################
#  Subset and write markers
###############################################################################

print("Writing markers...")
# write_markers(n_markers_per_type, marker_out)

###############################################################################
#  Visually check quality of markers
###############################################################################

#   Visually show how markers look for each cell type
plot_list <- lapply(
    unique(marker_stats$cellType.target),
    function(ct) {
        genes <- marker_stats %>%
            filter(
                rank_ratio <= n_markers_per_type,
                cellType.target == ct,
                ratio > 1
            ) %>%
            pull(gene)
        my_plotExpression(
            sce, genes,
            ct = cell_column,
            fill_colors = metadata(sce)[[colors_col]],
            title = paste("Top", length(genes), "for", ct),
            marker_stats = marker_stats %>%
                filter(
                    rank_ratio <= n_markers_per_type,
                    cellType.target == ct,
                    ratio > 1
                )
        )
    }
)

#   Write a multi-page PDF with violin plots for each cell group and all
#   markers
pdf(
    file.path(plot_dir, paste0("marker_gene_violin.pdf")),
    width = 35, height = 35
)
print(plot_list)
dev.off()

#   Plot mean ratio against log fold-change for all genes, split by target cell
#   type and colored by whether each gene will be used as a marker
p <- marker_stats %>%
    mutate(
        Marker = case_when(
            rank_ratio <= n_markers_per_type ~ paste0(
                "Top-", n_markers_per_type, " Marker"
            ),
            TRUE ~ "Not Marker"
        )
    ) %>%
    ggplot(aes(ratio, std.logFC, color = Marker)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    facet_wrap(~cellType.target, scales = "free_x", nrow = cell_type_nrow) +
    labs(x = "Mean Ratio") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    guides(col = guide_legend(override.aes = list(size = 2)))

pdf(file.path(plot_dir, paste0("mean_ratio_vs_1vall.pdf")), width = 10)
print(p)
dev.off()

#   Plot mean-ratio distibution by group (cell type or layer label)
boxplot_mean_ratio(n_markers_per_type, "mean_ratio_boxplot")

#   Get Ensembl ID for classical markers
spe <- readRDS(spe_IF_in)
stopifnot(all(classical_markers %in% rowData(sce)$gene_name))
classical_markers_ens <- rownames(sce)[
    match(classical_markers, rowData(sce)$gene_name)
]
stopifnot(all(classical_markers_ens %in% rownames(spe)))

#-------------------------------------------------------------------------------
#   First, plot the spatial expression of classical markers as reference
#-------------------------------------------------------------------------------

plot_list <- list()
plot_list_paper <- list()
i <- 1

for (j in 1:length(classical_markers)) {
    for (sample_id in unique(spe$sample_id)) {
        #   Determine the title for this subplot
        if (j <= 6) {
            title <- paste0(
                classical_markers[j], ": marker for layer ", j, "\n(",
                sample_id, ")"
            )
        } else if (classical_markers[j] == "SNAP25") {
            title <- paste0(
                "SNAP25: marker for gray matter\n(", sample_id, ")"
            )
        } else if (classical_markers[j] == "MOBP") {
            title <- paste0(
                "MOBP: marker for white matter\n(", sample_id, ")"
            )
        }

        #   Produce the ggplot object (grid version)
        plot_list[[i]] <- spot_plot(
            spe, sample_id = sample_id, var_name = classical_markers_ens[j],
            include_legend = TRUE, is_discrete = FALSE,
            title = title, assayname = "counts", minCount = 0
        )
        
        #   Produce the ggplot object (manuscript version)
        plot_list_paper[[i]] <-  spot_plot(
            spe, sample_id = sample_id, var_name = classical_markers_ens[j],
            include_legend = FALSE, is_discrete = FALSE,
            title = sub('\n', ': ', title), assayname = "counts", minCount = 0
        )

        i <- i + 1
    }
}

#   Create a grid version as a single-page PDF
n_sample <- length(unique(spe$sample_id))
pdf(
    file.path(
        plot_dir, paste0("marker_spatial_sparsity_reference.pdf")
    ),
    width = 7 * n_sample,
    height = 7 * length(classical_markers)
)
print(plot_grid(plotlist = plot_list, ncol = n_sample))
dev.off()

#   Create a manuscript-compatible version with one plot per PDF page
#   For figures in the paper, create a PDF version with one plot per
#   page. Exactly match the shape (aspect ratio) and visual details
#   of other spatial plots in this script
pdf(
    file.path(
        plot_dir, paste0("marker_spatial_sparsity_reference_paper.pdf")
    )
)
print(plot_list_paper)
dev.off()

#-------------------------------------------------------------------------------
#   For IF, show a grid of plots summarizing how sparsely marker genes
#   for each cell type are expressed spatially. Repeat these plots for different
#   numbers of markers per cell type: 15, 25, 50
#-------------------------------------------------------------------------------

for (n_markers in c(15, 25, 50)) {
    plot_list <- list()
    plot_list_paper <- list()
    i <- 1

    #   Plot proportion of markers having nonzero expression for each cell type
    for (ct in unique(marker_stats$cellType.target)) {
        #   Get markers for this cell type
        markers <- marker_stats %>%
            filter(
                cellType.target == ct,
                rank_ratio <= n_markers,
                ratio > 1
            ) %>%
            pull(gene)

        for (sample_id in unique(spe$sample_id)) {
            spe_small <- spe[markers, spe$sample_id == sample_id]

            #   For each spot, compute proportion of marker genes with nonzero
            #   expression
            spe_small$prop_nonzero_marker <- colMeans(
                assays(spe_small)$counts > 0
            )
            
            plot_list[[i]] <- spot_plot(
                spe_small, sample_id = sample_id,
                var_name = "prop_nonzero_marker", include_legend = TRUE,
                is_discrete = FALSE, minCount = 0,
                title = paste0(
                    "Prop. markers w/ nonzero exp:\n", ct, " (", sample_id, ")"
                )
            )
            
            # Use a 1-line title and remove the legend for the paper version   
            plot_list_paper[[i]] <- spot_plot(
                spe_small, sample_id = sample_id,
                var_name = "prop_nonzero_marker", include_legend = FALSE,
                is_discrete = FALSE,  minCount = 0,
                title = paste0(ct, " (", sample_id, ")")
            )
            
            i <- i + 1
        }
    }
    n_sample <- length(unique(spe$sample_id))
    n_rows <- length(unique(marker_stats$cellType.target))

    pdf(
        file.path(
            plot_dir, paste0("marker_spatial_sparsity_n", n_markers, ".pdf")
        ),
        width = 7 * n_sample, height = 7 * n_rows
    )
    print(plot_grid(plotlist = plot_list, ncol = n_sample))
    dev.off()
    
    #   Create a manuscript-compatible version with one plot per PDF page
    #   For figures in the paper, create a PDF version with one plot per
    #   page. Exactly match the shape (aspect ratio) and visual details
    #   of other spatial plots in this script
    pdf(
        file.path(
            plot_dir,
            paste0("marker_spatial_sparsity_n", n_markers, "_paper.pdf")
        )
    )
    print(plot_list_paper)
    dev.off()
}

#-------------------------------------------------------------------------------
#   For the manuscript, make a 4x4 plot for IF layer-level: 4 columns (samples)
#   by 4 rows: expression of PCP4 then that of 'Excit_L5' for different numbers
#   of markers: 15,25,50
#-------------------------------------------------------------------------------

if (cell_group == "layer") {
    plot_list <- list()
    max_list <- list()
    i <- 1

    #   Plot expression of PCP4 for every sample
    for (sample_id in unique(spe$sample_id)) {
        #   Produce the ggplot object
        plot_list[[i]] <- spot_plot(
            spe, sample_id = sample_id,
            var_name = classical_markers_ens[classical_markers == "PCP4"],
            include_legend = TRUE, is_discrete = FALSE, title = NULL,
            assayname = "counts", minCount = 0
        )
        
        max_list[[i]] = 0

        i <- i + 1
    }

    for (n_markers in c(15, 25, 50)) {
        #   Get markers for this cell type
        markers <- marker_stats %>%
            filter(
                cellType.target == "Excit_L5",
                rank_ratio <= n_markers,
                ratio > 1
            ) %>%
            pull(gene)

        for (sample_id in unique(spe$sample_id)) {
            spe_small <- spe[markers, spe$sample_id == sample_id]

            #   For each spot, compute proportion of marker genes with nonzero
            #   expression
            spe_small$prop_nonzero_marker <- colMeans(
                assays(spe_small)$counts > 0
            )
            
            max_list[[i]] <- max(spe_small$prop_nonzero_marker)
            
            plot_list[[i]] <- spot_plot(
                spe_small, sample_id = sample_id,
                var_name = "prop_nonzero_marker", include_legend = TRUE,
                is_discrete = FALSE, title = NULL, minCount = 0
            )
            
            i <- i + 1
        }
    }
    
    max_mat <- matrix(
        unlist(max_list), ncol = length(unique(spe$sample_id)), byrow = TRUE
    )
    
    #   Now loop back through the plot list (which will be displayed in 2D)
    #   and overwrite the scale to go as high as the largest value in the
    #   column. This allows for easy comparison between number of markers for
    #   the same samples
    for (i_col in 1:length(unique(spe$sample_id))) {
        for (i_row in 2:4) {
            index <- (i_row - 1) * length(unique(spe$sample_id)) + i_col
            upper_limit <- max(max_mat[, i_col])
            
            plot_list[[index]] <- plot_list[[index]] +
                scale_fill_gradientn(
                    colors = viridisLite::plasma(21),
                    limits = c(0, upper_limit), na.value = c("#CCCCCC40")
                ) +
                scale_color_gradientn(
                    colors = viridisLite::plasma(21),
                    limits = c(0, upper_limit), na.value = c("#CCCCCC40")
                ) 
        }
    }

    pdf(
        file.path(
            plot_dir, paste0("sparsity_figure.pdf")
        ),
        width = 7 * n_sample, height = 28
    )
    print(plot_grid(plotlist = plot_list, ncol = n_sample))
    dev.off()
}

session_info()
