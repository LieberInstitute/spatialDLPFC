library("here")
library("ggplot2")
library("jaffelab")
library("tidyverse")
library("reshape2")
library("spatialLIBD")
library("cowplot")
library("sessioninfo")

#   "IF" or "nonIF"
dataset = "IF"

raw_results_broad_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", dataset,
    "results_raw_broad.csv"
)
raw_results_layer_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", dataset,
    "results_raw_layer.csv"
)

if (dataset == "IF") {
    collapsed_results_broad_path <- here(
        "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
        "results_collapsed_broad.csv"
    )
    collapsed_results_layer_path <- here(
        "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
        "results_collapsed_layer.csv"
    )
    
    layer_ann_path <- here(
        "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
        "annotations_{sample_id}_spots.csv"
    )
} else {
    nonIF_id_path <- here(
        "processed-data", "spot_deconvo", "nonIF_ID_table.csv"
    )
    nonIF_counts_path <- here(
        "processed-data", "rerun_spaceranger", "{sample_id}", "outs", "spatial",
        "tissue_spot_counts.csv"
    )
}

sample_ids_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", dataset,
    "sample_ids.txt"
)

deconvo_tools <- c("Tangram", "Cell2location", "SPOTlight")

plot_dir <- here(
    "plots", "spot_deconvo", "05-shared_utilities", "combined", dataset
)

spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)
spe_nonIF_in <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
)

#   Only needed to get colors for software-estimated cell types
sce_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "sce_layer.rds"
)

cell_types_actual <- c("Astro", "Micro", "Neuron", "Oligo", "Other")
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

if (dataset == "IF") {
    shape_scale <- c(16, 17, 15, 3)
    names(shape_scale) <- sample_ids
}

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

if (dataset == "IF") {
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
    
    #   Use titlecase (note this is MUCH faster than mutate with 'toTitleCase')
    observed_df_long$deconvo_tool[
        observed_df_long$deconvo_tool == "tangram"
    ] <- "Tangram"
    observed_df_long$deconvo_tool[
        observed_df_long$deconvo_tool == "cell2location"
    ] <- "Cell2location"
    
    #   Clean up labels
    observed_df_long$label <- tolower(observed_df_long$label)
    observed_df_long$label <- sub("layer", "Layer ", observed_df_long$label)
    observed_df_long$label[observed_df_long$label == "wm"] <- "WM"
    stopifnot(
        all(unlist(corresponding_layers) %in% unique(observed_df_long$label))
    )
    
    ############################################################################
    #   Broad-Resolution Plots
    ############################################################################
    
    counts_df <- observed_df_long |>
        #   For each manually annotated label, deconvo tool and sample_id,
        #   normalize by the total counts of all cell types
        group_by(deconvo_tool, label, sample_id, resolution) |>
        mutate(count = count / sum(count)) |>
        #   Now for each label, deconvo tool, sample_id and cell type, add up
        #   counts for all relevant spots
        group_by(deconvo_tool, label, cell_type, sample_id, resolution) |>
        summarize(count = sum(count)) |>
        #   Now average across samples
        group_by(deconvo_tool, label, cell_type, resolution) |>
        summarize(count = mean(count)) |>
        ungroup()
    
    #---------------------------------------------------------------------------
    #   Barplots comparing distribution of cell types in layers between
    #   resolutions
    #---------------------------------------------------------------------------
    
    p <- ggplot(counts_df, aes(x = resolution, y = count, fill = cell_type)) +
        facet_grid(rows = vars(deconvo_tool), cols = vars(label)) +
        geom_bar(stat = "identity") +
        labs(
            x = "Cell-Type Resolution", y = "Proportion of Counts",
            fill = "Cell Type"
        ) +
        scale_fill_manual(values = estimated_cell_labels) +
        scale_x_discrete(labels = c("broad" = "Broad", "layer" = "Layer")) +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    pdf(file.path(plot_dir, "barplot.pdf"), width = 10, height = 5)
    print(p)
    dev.off()
} else {
    #   This type of plot isn't possible for the nonIF data: just load the SPE
    #   object and continue
    load(spe_nonIF_in, verbose = TRUE)
    
    #   A table mapping one version of sample IDs to another
    id_table <- read.csv(nonIF_id_path)
    
    #   Add in the up-to-date VistoSeg counts
    spe$count <- NA
    for (sample_id in unique(spe$sample_id)) {
        #   Correctly determine the path for the cell counts for this sample,
        #   then read in
        long_id <- id_table[match(sample_id, id_table$short_id), "long_id"]
        this_path <- sub(
            "{sample_id}", long_id, nonIF_counts_path, fixed = TRUE
        )
        cell_counts <- read.csv(this_path)
        
        #   All spots in the object should have counts
        stopifnot(
            all(
                colnames(spe[, spe$sample_id == sample_id]) %in%
                    cell_counts$barcode
            )
        )
        
        #   Line up the rows of 'cell_counts' with the sample-subsetted SPE
        #   object
        cell_counts <- cell_counts[
            match(
                colnames(spe[, spe$sample_id == sample_id]),
                cell_counts$barcode
            ),
        ]
        
        #   Add this sample's counts to the SPE object
        spe$count[spe$sample_id == sample_id] <- cell_counts$Nmask_dark_blue
    }
    
    #   Ensure counts were read in for all spots in the object
    if (any(is.na(spe$count))) {
        stop("Did not find cell counts for all non-IF spots.")
    }
    
    #   Merge results for both resolutions
    observed_df_long <- rbind(observed_df_broad, observed_df_layer) |>
        as_tibble()
    
    #   Use titlecase (note this is MUCH faster than mutate with 'toTitleCase')
    observed_df_long$deconvo_tool[
        observed_df_long$deconvo_tool == "tangram"
    ] <- "Tangram"
    observed_df_long$deconvo_tool[
        observed_df_long$deconvo_tool == "cell2location"
    ] <- "Cell2location"
}

#-------------------------------------------------------------------------------
#   Scatterplots comparing section-wide cell-type proportions between
#   resolutions
#-------------------------------------------------------------------------------

sample_prop_scatter = function(
        counts_df, dataset, color_scale, filename, x_angle = 0
    ) {
    metrics_df <- counts_df |>
        group_by(deconvo_tool) |>
        summarize(
            corr = paste0("Cor = ", round(cor(broad, layer), 2)),
            rmse = paste0("RMSE = ", signif(mean((broad - layer)**2)**0.5, 3))
        ) |>
        ungroup()
    
    counts_df <- counts_df |>
        mutate(broad = log(broad), layer = log(layer))
    
    #   We'll shape by sample for IF data. For nonIF, there are too many samples
    if (dataset == "IF") {
        p <- ggplot(
            counts_df,
            aes(x = broad, y = layer, color = cell_type, shape = sample_id)
        )
    } else {
        p <- ggplot(counts_df, aes(x = broad, y = layer, color = cell_type))
    }
    
    p <- p + geom_point() +
        geom_abline(
            intercept = 0, slope = 1, linetype = 3, color = "black"
        ) +
        facet_wrap(~deconvo_tool) +
        scale_color_manual(values = color_scale) +
        #   Correlation label
        geom_text(
            data = metrics_df,
            mapping = aes(
                x = max(counts_df$broad),
                y = 0.1 * max(counts_df$layer) + 0.9 * min(counts_df$layer),
                label = corr,
                color = NULL,
                shape = NULL
            ),
            hjust = 1, vjust = 0, show.legend = FALSE
        ) +
        #   RMSE label
        geom_text(
            data = metrics_df,
            mapping = aes(
                x = max(counts_df$broad),
                y = min(counts_df$layer),
                label = rmse,
                color = NULL,
                shape = NULL
            ),
            hjust = 1, vjust = 0, show.legend = FALSE
        ) +
        labs(
            x = "log(Total Broad Counts)", y = "log(Total Layer-Level Counts)",
            color = "Cell Type", shape = "Sample ID"
        ) +
        theme_bw(base_size = 15)
    
        #   Account for misalignment of labels that occurs for certain
        #   rotations
        if (x_angle == 90 || x_angle == 270) {
            p <- p +
                theme(axis.text.x = element_text(angle = x_angle, vjust = 0.5))
        } else {
            p <- p +
                theme(axis.text.x = element_text(angle = x_angle))
        }
        
    
    pdf(
        file.path(plot_dir, filename),
        width = 9, height = 3
    )
    print(p)
    dev.off()
}

#   Plot at broad resolution
counts_df <- observed_df_long |>
    #   Sum counts across each section
    group_by(deconvo_tool, sample_id, cell_type, resolution) |>
    summarize(count = sum(count)) |>
    #   Now pivot wider by resolution
    pivot_wider(names_from = resolution, values_from = count) |>
    ungroup()

sample_prop_scatter(
    counts_df, dataset, estimated_cell_labels,
    "sample_proportions_scatter_broad.pdf", x_angle = 90
)

#   Plot at collapsed resolution

if (dataset == "IF") {
    cols_to_select = c(
        added_colnames, "resolution", "label",
        cell_types_actual[- match('Other', cell_types_actual)]
    )
} else {
    cols_to_select = c(
        added_colnames, "resolution",
        cell_types_actual[- match('Other', cell_types_actual)]
    )
}

counts_df = observed_df_long |>
    pivot_wider(names_from = cell_type, values_from = count) |>
    #   Collapse cell types
    mutate(Oligo = Oligo + OPC, Neuron = Excit + Inhib) |>
    select(all_of(cols_to_select)) |>
    #   Change back to one cell type per row
    pivot_longer(
        cols = cell_types_actual[- match('Other', cell_types_actual)],
        names_to = "cell_type", values_to = "count"
    ) |>
    #   Sum counts across each section
    group_by(deconvo_tool, sample_id, cell_type, resolution) |>
    summarize(count = sum(count)) |>
    #   Now pivot wider by resolution
    pivot_wider(names_from = resolution, values_from = count) |>
    ungroup()

sample_prop_scatter(
    counts_df, dataset, cell_type_labels, 
    "sample_proportions_scatter_collapsed.pdf", x_angle = 90
)

################################################################################
#   Total-Cells Plots
################################################################################

if (dataset == "IF") {
    collapsed_broad = read.csv(collapsed_results_broad_path) |>
        mutate(resolution = "broad")
    
    collapsed_results = read.csv(collapsed_results_layer_path) |>
        mutate(resolution = "layer") |>
        rbind(collapsed_broad) |>
        #   Make one cell type per row
        pivot_longer(
            cols = all_of(tolower(cell_types_actual)), names_to = "cell_type",
            values_to = "count"
        ) |>
        mutate(cell_type = str_to_title(cell_type)) |>
        #   Have an observed and actual column
        pivot_wider(
            names_from = obs_type, values_from = count,
        )
    
    #   Use titlecase (note this is MUCH faster than mutate with 'toTitleCase')
    collapsed_results$deconvo_tool[
        collapsed_results$deconvo_tool == "tangram"
    ] <- "Tangram"
    collapsed_results$deconvo_tool[
        collapsed_results$deconvo_tool == "cell2location"
    ] <- "Cell2location"
    
    #   Calculate total cells per spot and prepare for plotting
    count_df = collapsed_results |>
        filter(deconvo_tool %in% c("Tangram", "Cell2location")) |>
        #   Add counts of any cell type for each spot
        group_by(deconvo_tool, resolution, barcode, sample_id) |>
        summarize(observed = sum(observed), actual = sum(actual))

} else {
    #   Data frame of VistoSeg counts
    visto_df = as_tibble(
        cbind(colnames(spe), colData(spe)[, c("sample_id", "count")])
    )
    colnames(visto_df) = c("barcode", "sample_id", "actual")
    
    #   Add up counts of all cell types per spot
    count_df <- observed_df_long |>
        filter(deconvo_tool %in% c("Tangram", "Cell2location")) |>
        group_by(deconvo_tool, sample_id, barcode, resolution) |>
        summarize(observed = sum(count)) |>
        left_join(visto_df)
}

#   Compute metrics for each deconvolution tool and resolution: correlation
#   between observed and actual values as well as RMSE
metrics_df <- count_df |>
    group_by(deconvo_tool, resolution) |>
    summarize(
        corr = paste("Cor =", round(cor(observed, actual), 2)),
        rmse = paste("RMSE =", signif(mean((observed - actual)**2)**0.5, 3))
    ) |>
    ungroup()

p = ggplot(count_df, aes(x = observed, y = actual)) +
    geom_point(alpha = 0.01) +
    facet_grid(
        rows = vars(resolution), cols = vars(deconvo_tool),
        labeller = labeller(
            resolution = c("broad" = "Broad", "layer" = "Layer")
        )
    ) +
    geom_abline(
        intercept = 0, slope = 1, linetype = "dashed", color = "red"
    ) +
    geom_text(
        data = metrics_df,
        mapping = aes(
            x = max(count_df$observed), y = max(count_df$actual) / 7,
            label = corr
        ),
        hjust = 1, size = 5
    ) +
    geom_text(
        data = metrics_df,
        mapping = aes(x = max(count_df$observed), y = 0, label = rmse),
        hjust = 1, vjust = 0, size = 5
    ) +
    labs(
        x = "Calculated Cell Count",
        y = "Provided Cell Count (Cellpose)",
        title = "Provided vs. Calculated Total Cells Per Spot"
    ) +
    theme_bw(base_size = 15)

pdf(file.path(plot_dir, "total_cells_spot_paper.pdf"))
print(p)
dev.off()

session_info()
