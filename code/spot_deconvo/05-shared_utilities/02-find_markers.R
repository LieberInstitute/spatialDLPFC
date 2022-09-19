#   Find marker genes that will be shared for both the IF and non-IF analyses.

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("HDF5Array"))

cell_group = "layer" # "broad" or "layer"

#   Number of marker genes to use per cell type
n_markers_per_type <- 25

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

plot_dir <- here(
    "plots", "spot_deconvo", "05-shared_utilities", cell_group
)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
#  Load objects and preprocess 'marker_stats'
###############################################################################

sce = readRDS(sce_in)
marker_stats <- readRDS(marker_object_in)
gc()

print(paste0("Running script at ", cell_group, "-resolution."))

#-------------------------------------------------------------------------------
#   Filter out mitochondrial genes and re-rank 'rank_ratio' values. Add gene
#   symbols to 'marker_stats' object
#-------------------------------------------------------------------------------

#   Add gene symbol
marker_stats$symbol = rowData(sce)$gene_name[
    match(marker_stats$gene, rownames(sce))
]

#   Filter out mitochondrial genes
marker_stats = marker_stats[!grepl('^MT-', marker_stats$symbol),]

#   "Re-rank" rank_ratio, since there may be missing ranks now
for (ct in unique(marker_stats$cellType.target)) {
    old_ranks = marker_stats %>%
        filter(cellType.target == ct) %>%
        pull(rank_ratio) %>%
        sort()
    
    for (i in 1:length(which((marker_stats$cellType.target == ct)))) {
        index = which(
            (marker_stats$cellType.target == ct) &
                (marker_stats$rank_ratio == old_ranks[i])
        )
        stopifnot(length(index) == 1)
        marker_stats[index, 'rank_ratio'] = i
    }
}

###############################################################################
#  Functions
###############################################################################

write_markers = function(n_markers, out_path) {
    #   Take top N marker genes for each cell type
    marker_stats_temp <- marker_stats %>%
        filter(
            rank_ratio <= n_markers,
            ratio > 1
        )
    
    #   Warn if less than the intended number of markers is used for any cell
    #   type
    num_markers_table = marker_stats_temp %>%
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

my_plotExpression <- function(
        sce, genes, assay = "logcounts", ct = "cellType", fill_colors = NULL,
        title = NULL, marker_stats
    ) {
    cat_df <- as.data.frame(colData(sce))[, ct, drop = FALSE]
    expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes, ]))
    
    cat <- cat_df[expression_long$Var2, ]
    expression_long <- cbind(expression_long, cat)

    #   Use gene symbols for labels, not Ensembl ID
    symbols = rowData(sce)$gene_name[match(genes, rownames(sce))]
    names(symbols) = genes
    
    #   Add a data frame for adding mean-ratio labels to each gene
    text_df = marker_stats
    text_df$ratio = paste0('Mean ratio: ', round(text_df$ratio, 2))
    text_df$Var1 = text_df$gene
    
    expression_violin <- ggplot(
        data = expression_long, aes(x = cat, y = value, fill = cat)
    ) +
        geom_violin(scale = "width") +
        geom_text(
            data = text_df,
            mapping = aes(
                x = length(unique(sce[[ct]])), y = Inf, fill = NULL,
                label = ratio, size = 10
            ),
            hjust = 1, vjust = 1
        ) +
        facet_wrap(
            ~Var1, ncol = 5, scales = "free_y",
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

#   Plot mean-ratio distribution by cell type/ layer
boxplot_mean_ratio = function(n_markers, plot_name) {
    p = marker_stats %>%
        filter(rank_ratio <= n_markers) %>%
        mutate(ratio, ratio = log(ratio)) %>%
        ggplot(aes(cellType.target, ratio, color = cellType.target)) +
            geom_boxplot() +
            geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
            scale_color_manual(values = metadata(sce)$cell_type_colors_layer) +
            labs(y = "log(Mean Ratio)") +
            theme_bw() +
            guides(color = "none")
    
    ggsave(
        p, filename = file.path(plot_dir, paste0(plot_name, ".png")),
        height = 10, width = 10
    )
}

###############################################################################
#  Subset and write markers
###############################################################################

print("Writing markers...")
write_markers(n_markers_per_type, marker_out)

###############################################################################
#  Visually check quality of markers
###############################################################################

if (cell_group == "broad") {
    cell_column = "cellType_broad_hc"
} else {
    cell_column = "layer_level"
}

#   Visually show how markers look for each cell type
plot_list = lapply(
    unique(marker_stats$cellType.target),
    function(ct){
        genes = marker_stats %>%
            filter(
                rank_ratio <= n_markers_per_type,
                cellType.target == ct,
                ratio > 1
            ) %>%
            pull(gene)
        my_plotExpression(
            sce, genes, ct = cell_column,
            fill_colors = metadata(sce)$cell_type_colors_layer,
            title = paste('Top', length(genes), 'for', ct),
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
p = marker_stats %>% 
    mutate(
        Marker = case_when(
            rank_ratio <= n_markers_per_type ~ paste0('Marker top', n_markers_per_type),
            TRUE ~ 'Not marker'
        )
    ) %>%
    ggplot(aes(ratio, std.logFC, color = Marker)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
    geom_vline(xintercept = 1, linetype = 'dashed', color = 'red') +
    facet_wrap(~cellType.target, scales = "free_x") +
    labs(x = "Mean Ratio") +
    theme_bw()

ggsave(
    p, filename = file.path(plot_dir, paste0("mean_ratio_vs_1vall.png")),
    height = 10, width = 10
)

#   Plot mean-ratio distibution by group (cell type or layer label)
boxplot_mean_ratio(n_markers_per_type, "mean_ratio_boxplot")

session_info()
