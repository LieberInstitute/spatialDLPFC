library("here")
library("ggplot2")
library("jaffelab")
library("tidyverse")
library("reshape2")
library("spatialLIBD")
library("cowplot")
library("sessioninfo")

raw_results_broad_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "results_raw_broad.csv"
)
raw_results_layer_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "results_raw_layer.csv"
)

sample_ids_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "sample_ids.txt"
)

deconvo_tools <- c("tangram", "cell2location", "SPOTlight")

plot_dir <- here(
    "plots", "spot_deconvo", "05-shared_utilities", "combined"
)

layer_ann_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "annotations_{sample_id}_spots.csv"
)

spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)

#   Only needed to get colors for software-estimated cell types
sce_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "sce_layer.rds"
)

cell_types_actual <- c("astro", "micro", "neuron", "oligo", "other")
cell_types_broad <- c(
    "Astro", "EndoMural", "Excit", "Inhib", "Micro", "Oligo", "OPC"
)
cell_types_layer <- c(
    "Astro", "EndoMural", "Excit_L2_3", "Excit_L3", "Excit_L3_4_5",
    "Excit_L4", "Excit_L5", "Excit_L5_6", "Excit_L6", "Inhib", "Micro",
    "Oligo", "OPC"
)


#   Make a list of which layers we expect each cell type to be most highly
#   expressed in
corresponding_layers <- list(
    "Astro" = "Layer 1",
    "EndoMural" = "Layer 1",
    "Excit" = paste("Layer", 2:6),
    "Excit_L2_3" = c("Layer 2", "Layer 3"),
    "Excit_L3" = "Layer 3",
    "Excit_L3_4_5" = c("Layer 3", "Layer 4", "Layer 5"),
    "Excit_L4" = "Layer 4",
    "Excit_L5" = "Layer 5",
    "Excit_L5_6" = c("Layer 5", "Layer 6"),
    "Excit_L6" = "Layer 6",
    "Inhib" = paste("Layer", 2:6),
    "Micro" = c("Layer 1", "WM"),
    "Oligo" = "WM",
    "OPC" = c("Layer 1", "WM")
)

#   Name spatialLIBD colors with the layer names used in this script
names(libd_layer_colors)[
    match(c(paste0("Layer", 1:6), "WM"), names(libd_layer_colors))
] <- c(paste("Layer", 1:6), "WM")

cell_type_labels <- c(
    "#3BB273", "#663894", "#E49AB0", "#E07000", "#95B8D1"
)
names(cell_type_labels) <- c(cell_types_actual)

set.seed(11282022)
dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Read in and format cell counts into a table apt for plotting with ggplot
################################################################################

sample_ids <- readLines(sample_ids_path)
added_colnames <- c("barcode", "sample_id", "deconvo_tool")

shape_scale <- c(16, 17, 15, 3)
names(shape_scale) <- sample_ids

#-------------------------------------------------------------------------------
#   Read in broad results (raw)
#-------------------------------------------------------------------------------

#   Convert to long format
observed_df_broad <- read.csv(raw_results_broad_path) |>
    melt(
        id.vars = added_colnames, variable.name = "cell_type",
        value.name = "count"
    ) |>
    mutate(resolution = "broad")

#-------------------------------------------------------------------------------
#   Read in layer-level results (raw)
#-------------------------------------------------------------------------------

observed_df_layer <- read.csv(raw_results_layer_path)

#   Collapse excitatory cells to match broad resolution
observed_df_layer$Excit <- observed_df_layer |>
    select(starts_with("Excit_L")) |>
    rowSums()

#   Convert to long resolution
observed_df_layer <- observed_df_layer |>
    select(all_of(c(added_colnames, cell_types_broad))) |>
    melt(
        id.vars = added_colnames, variable.name = "cell_type",
        value.name = "count"
    ) |>
    mutate(resolution = "layer") |>
    as_tibble()


#   Ensure data for both resolutions have the same columns and that the cell
#   types are now identical
stopifnot(identical(colnames(observed_df_layer), colnames(observed_df_broad)))
stopifnot(
    identical(
        unique(observed_df_layer$cell_type), unique(observed_df_broad$cell_type)
    )
)

################################################################################
#   Read in manual layer annotation and load the SPE object
################################################################################

#   Read in SCE just to get colors for software-estimated cell types
sce <- readRDS(sce_in)
estimated_cell_labels <- metadata(sce)$cell_type_colors_broad

spe <- readRDS(spe_IF_in)

layer_ann <- data.frame()
for (sample_id in sample_ids) {
    this_layer_path <- sub("\\{sample_id\\}", sample_id, layer_ann_path)
    layer_ann_small <- read.csv(this_layer_path)

    layer_ann_small$barcode <- colnames(spe[, spe$sample_id == sample_id])[
        layer_ann_small$id + 1
    ]
    layer_ann_small$sample_id <- sample_id

    layer_ann <- rbind(layer_ann, layer_ann_small)
}
layer_ann$id <- NULL

#   Add layer label to observed_df_long
observed_df_long <- rbind(observed_df_broad, observed_df_layer) |>
    left_join(layer_ann, by = c("barcode", "sample_id")) |>
    filter(!is.na(label)) |>
    as_tibble()

#   Clean up labels
observed_df_long$label <- tolower(observed_df_long$label)
observed_df_long$label <- sub("layer", "Layer ", observed_df_long$label)
observed_df_long$label[observed_df_long$label == "wm"] <- "WM"
stopifnot(
    all(unlist(corresponding_layers) %in% unique(observed_df_long$label))
)

################################################################################
#   Plots
################################################################################

counts_df <- observed_df_long |>
    #   For each manually annotated label, deconvo tool and sample_id, normalize
    #   by the total counts of all cell types
    group_by(deconvo_tool, label, sample_id, resolution) |>
    mutate(count = count / sum(count)) |>
    #   Now for each label, deconvo tool, sample_id and cell type, add up counts
    #   for all relevant spots
    group_by(deconvo_tool, label, cell_type, sample_id, resolution) |>
    summarize(count = sum(count)) |>
    #   Now average across samples
    group_by(deconvo_tool, label, cell_type, resolution) |>
    summarize(count = mean(count)) |>
    ungroup()

#-------------------------------------------------------------------------------
#   Barplots comparing distribution of cell types in layers between resolutions
#-------------------------------------------------------------------------------

p <- ggplot(counts_df, aes(x = resolution, y = count, fill = cell_type)) +
    facet_grid(rows = vars(deconvo_tool), cols = vars(label)) +
    geom_bar(stat = "identity") +
    labs(
        x = "Cell-Type Resolution", y = "Proportion of Counts",
        fill = "Cell Type"
    ) +
    scale_fill_manual(values = estimated_cell_labels) +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

pdf(file.path(plot_dir, "barplot.pdf"), width = 10, height = 5)
print(p)
dev.off()

#-------------------------------------------------------------------------------
#   Scatterplots comparing section-wide cell-type proportions between
#   resolutions
#-------------------------------------------------------------------------------

counts_df <- observed_df_long |>
    #   Sum counts across each section
    group_by(deconvo_tool, sample_id, cell_type, resolution) |>
    summarize(count = sum(count)) |>
    #   Now pivot wider by resolution
    pivot_wider(names_from = resolution, values_from = count) |>
    ungroup()

metrics_df <- counts_df |>
    group_by(deconvo_tool) |>
    summarize(
        corr = paste0("Cor = ", round(cor(broad, layer), 2)),
        rmse = paste0("RMSE = ", signif(mean((broad - layer)**2)**0.5, 3))
    ) |>
    ungroup()

p <- ggplot(
    counts_df,
    aes(x = broad, y = layer, color = cell_type, shape = sample_id)
) +
    geom_point() +
    geom_abline(
        intercept = 0, slope = 1, linetype = "dashed", color = "red"
    ) +
    facet_wrap(~deconvo_tool) +
    scale_color_manual(values = estimated_cell_labels) +
    #   Correlation label
    geom_text(
        data = metrics_df,
        mapping = aes(
            x = max(counts_df$broad),
            y = max(counts_df$layer),
            label = corr,
            color = NULL,
            shape = NULL
        ),
        hjust = 1, vjust = 1, show.legend = FALSE
    ) +
    #   RMSE label
    geom_text(
        data = metrics_df,
        mapping = aes(
            x = max(counts_df$broad),
            y = 0.9 * max(counts_df$layer) + 0.1 * min(counts_df$layer),
            label = rmse,
            color = NULL,
            shape = NULL
        ),
        hjust = 1, vjust = 1, show.legend = FALSE
    ) +
    labs(
        x = "Total Broad Counts", y = "Total Layer-Level Counts",
        color = "Cell Type", shape = "Sample ID"
    ) +
    theme_bw(base_size = 15)

pdf(
    file.path(plot_dir, "sample_proportions_scatter.pdf"),
    width = 9, height = 3
)
print(p)
dev.off()

session_info()
