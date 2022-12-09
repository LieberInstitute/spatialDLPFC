library("here")
library("ggplot2")
library("jaffelab")
library("tidyverse")
library("reshape2")
library("spatialLIBD")
library("cowplot")

#   Adds the 'spot_plot' function, a wrapper for 'vis_gene' or 'vis_clus' with
#   consistent manuscript-appropriate settings
source('shared_functions.R')

cell_group <- "layer" # "broad" or "layer"

sample_ids_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "sample_ids.txt"
)

deconvo_tools <- c("Tangram", "Cell2location", "SPOTlight")

plot_dir <- here(
    "plots", "spot_deconvo", "05-shared_utilities", cell_group
)

processed_dir <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF"
)

raw_results_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    paste0("results_raw_", cell_group, ".csv")
)

collapsed_results_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    paste0("results_collapsed_", cell_group, ".csv")
)

spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)

#   Only needed to get colors for software-estimated cell types
sce_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("sce_", cell_group, ".rds")
)

marker_object_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("marker_stats_", cell_group, ".rds")
)

marker_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("markers_", cell_group, ".txt")
)

layer_ann_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "annotations_{sample_id}_spots.csv"
)

cell_types_actual <- c("Astro", "Micro", "Neuron", "Oligo", "Other")
if (cell_group == "broad") {
    cell_types <- c(
        "Astro", "EndoMural", "Excit", "Inhib", "Micro", "Oligo", "OPC"
    )
} else {
    cell_types <- c(
        "Astro", "EndoMural", "Excit_L2_3", "Excit_L3", "Excit_L3_4_5",
        "Excit_L4", "Excit_L5", "Excit_L5_6", "Excit_L6", "Inhib", "Micro",
        "Oligo", "OPC"
    )
}

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
    "#3BB273", "#663894", "#E49AB0", "#E07000", "#95B8D1", "#000000"
)
names(cell_type_labels) <- c(cell_types_actual, "average")

set.seed(11282022)

################################################################################
#   Plotting functions
################################################################################

#   Scatterplots observed vs. actual cell counts for each cell type, faceted
#   by sample and deconvolution tool. Use all spots as points
all_spots <- function(count_df, plot_name) {
    #   Compute metrics for each deconvolution tool: correlation between
    #   observed and actual values as well as RMSE
    metrics_df <- count_df %>%
        group_by(deconvo_tool, sample_id, cell_type) %>%
        summarize(
            corr = round(cor(observed, actual), 2),
            rmse = signif(mean((observed - actual)**2)**0.5, 3)
        ) %>%
        ungroup()

    #   Improve labels for plotting
    metrics_df$corr <- paste("Cor =", metrics_df$corr)
    metrics_df$rmse <- paste("RMSE =", metrics_df$rmse)

    plot_list <- lapply(
        cell_types_actual,
        function(ct) {
            count_df_small <- count_df %>%
                filter(cell_type == ct)

            p <- ggplot(count_df_small) +
                geom_point(aes(x = observed, y = actual), alpha = 0.01) +
                geom_abline(
                    intercept = 0, slope = 1, linetype = "dashed", color = "red"
                ) +
                facet_grid(
                    rows = vars(sample_id), cols = vars(deconvo_tool)
                ) +
                guides(col = guide_legend(override.aes = list(alpha = 1))) +
                labs(
                    title = ct,
                    x = "Software-Estimated",
                    y = "CART-Calculated",
                ) +
                geom_text(
                    data = metrics_df %>% filter(cell_type == ct),
                    mapping = aes(
                        x = Inf, y = max(count_df_small$actual) / 5,
                        label = corr
                    ),
                    hjust = 1
                ) +
                geom_text(
                    data = metrics_df %>% filter(cell_type == ct),
                    mapping = aes(x = Inf, y = 0, label = rmse),
                    hjust = 1, vjust = 0
                ) +
                theme_bw(base_size = 10)

            return(p)
        }
    )

    pdf(file.path(plot_dir, plot_name))
    print(plot_list)
    dev.off()
}

#   Scatterplot of observed vs. actual total counts summed across spots,
#   faceted by deconvolution tool
across_spots <- function(count_df, plot_name, x_angle = 0) {
    #   Compute metrics for each deconvolution tool: correlation between
    #   observed and actual values as well as RMSE
    metrics_df <- count_df %>%
        group_by(deconvo_tool) %>%
        summarize(
            corr = round(cor(observed, actual), 2),
            rmse = signif(mean((observed - actual)**2)**0.5, 3)
        ) %>%
        ungroup()

    #   Add KL divergence from observed section-wide cell-type proportions to
    #   ground-truth proportions, averaged across sections
    metrics_df <- count_df |>
        #   Ensure counts or proportions are normalized to add to 1 across
        #   cell types
        group_by(sample_id, deconvo_tool) |>
        mutate(
            observed = observed / sum(observed),
            actual = actual / sum(actual),
        ) |>
        #   Compute each term in the sum for KL divergence
        group_by(sample_id, deconvo_tool, cell_type) |>
        summarize(kl_piece = actual * log(actual / observed)) |>
        #   Add all terms to form the sum for each sample
        group_by(sample_id, deconvo_tool) |>
        summarize(kl = sum(kl_piece)) |>
        #   Take the mean across samples to form one value per tool
        group_by(deconvo_tool) |>
        summarize(kl = round(mean(kl), 2)) |>
        left_join(metrics_df)

    #   Improve labels for plotting
    metrics_df$corr <- paste("Cor =", metrics_df$corr)
    metrics_df$rmse <- paste("RMSE =", metrics_df$rmse)
    metrics_df$kl <- paste("KL Div. =", metrics_df$kl)

    p <- ggplot(count_df) +
        geom_point(
            aes(x = observed, y = actual, shape = sample_id, color = cell_type)
        ) +
        facet_wrap(~deconvo_tool) +
        geom_abline(
            intercept = 0, slope = 1, linetype = "dashed", color = "red"
        ) +
        #   Correlation label
        geom_text(
            data = metrics_df,
            mapping = aes(
                x = max(count_df$observed),
                y = 0.05 * max(count_df$actual),
                label = corr
            ),
            hjust = 1, vjust = 0
        ) +
        #   RMSE label
        geom_text(
            data = metrics_df,
            mapping = aes(
                x = max(count_df$observed),
                y = 0.15 * max(count_df$actual),
                label = rmse
            ),
            hjust = 1, vjust = 0
        ) +
        #   KL divergence label
        geom_text(
            data = metrics_df,
            mapping = aes(
                x = max(count_df$observed),
                y = 0.25 * max(count_df$actual),
                label = kl
            ),
            hjust = 1, vjust = 0
        ) +
        scale_color_manual(values = cell_type_labels) +
        scale_x_continuous(
            limits = c(0, max(count_df$observed) * 1.05),
            expand = c(0, 0)
        ) +
        scale_y_continuous(
            limits = c(0, max(count_df$actual) * 1.05),
            expand = c(0, 0)
        ) +
        labs(
            x = "Software-Estimated", y = "CART-Calculated",
            color = "Cell Type", shape = "Sample ID"
        ) +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = x_angle))

    pdf(file.path(plot_dir, plot_name), height = 4, width = 10)
    print(p)
    dev.off()
}

#   Return a length-2 list containing a 'vis_grid_gene' plot and maximum value
#   in the plot given the SpatialExperiment, long-format
#   tibble of cell-type counts, the target sample ID, deconvo tool, and cell
#   type, column name of 'full_df' ('observed' or 'actual'), and a plot title
spatial_counts_plot <- function(spe_small, full_df, sample_id1, deconvo_tool1, cell_type1, c_name,
    title) {
    #   Grab counts for just this sample, deconvo tool, and cell type
    counts_df <- full_df %>%
        filter(
            sample_id == sample_id1,
            deconvo_tool == deconvo_tool1,
            cell_type == cell_type1
        )

    #   Add counts as a column in colData(spe_small)
    spe_small$temp_ct_counts <- counts_df[[c_name]][
        match(colnames(spe_small), counts_df$barcode)
    ]

    #   Plot spatial distribution
    p <- spot_plot(
        spe_small, sample_id = sample_id1, title = title,
        var_name = "temp_ct_counts", include_legend = TRUE,
        is_discrete = FALSE
    )

    return(list(p, max(spe_small$temp_ct_counts)))
}

#   For all samples individually, write PDF plots of each cell type spatially
#   for each deconvolution tool.
#
#   spe: SpatialExperiment including all spatial samples
#   full_df: tibble in long format containing columns "barcode", "sample_id",
#       "deconvo_tool", "cell_type", "observed", and, if include_actual,
#       "actual"
#   cell_type_vec: character vector of unique cell types expected in the full_df
#       column "cell_type"
#   include_actual: length-1 logical: add a row of plots for the ground-truth
#       counts?
#   pdf_prefix: length-1 character to prepend to PDF filenames
spatial_counts_plot_full <- function(spe, full_df, cell_type_vec, include_actual, pdf_prefix) {
    for (sample_id in sample_ids) {
        spe_small <- spe[, spe$sample_id == sample_id]

        i <- 1
        plot_list <- list()
        max_list <- list()

        #   For each deconvo tool, make a row of plots (each including all cell
        #   types) showed the observed distribution
        for (deconvo_tool in deconvo_tools) {
            for (cell_type in cell_type_vec) {
                temp <- spatial_counts_plot(
                    spe_small, full_df, sample_id, deconvo_tool, cell_type,
                    "observed",
                    paste0(cell_type, " Counts (", deconvo_tool, ")")
                )
                plot_list[[i]] <- temp[[1]]
                max_list[[i]] <- temp[[2]]
                i <- i + 1
            }
        }

        #   Add a row showing the ground-truth counts for each cell type
        if (include_actual) {
            for (cell_type in cell_types_actual) {
                temp <- spatial_counts_plot(
                    spe_small, full_df, sample_id, deconvo_tool,
                    str_to_title(cell_type), "actual",
                    paste0(cell_type, " Counts (CART-Calculated)")
                )
                plot_list[[i]] <- temp[[1]]
                max_list[[i]] <- temp[[2]]
                i <- i + 1
            }
        }

        max_mat <- matrix(
            unlist(max_list),
            ncol = length(cell_type_vec), byrow = TRUE
        )

        #   Now loop back through the plot list (which will be displayed in 2D)
        #   and overwrite the scale to go as high as the largest value in the
        #   column. This allows for easy comparison between deconvo tools
        #   (and optionally the ground truth)
        for (i_col in 1:length(cell_type_vec)) {
            for (i_row in 1:(length(deconvo_tools) + include_actual)) {
                index <- (i_row - 1) * length(cell_type_vec) + i_col
                upper_limit <- max(max_mat[, i_col])

                plot_list[[index]] <- plot_list[[index]] +
                    scale_color_continuous(
                        type = "viridis", limits = c(0, upper_limit),
                        na.value = c("black" = "#0000002D")
                    ) +
                    scale_fill_continuous(
                        type = "viridis", limits = c(0, upper_limit),
                        na.value = c("black" = "#0000002D")
                    )
            }
        }

        #   Plot in a grid where cell types are columns and rows are
        #   deconvolution tools. One PDF per sample
        pdf(
            file.path(plot_dir, paste0(pdf_prefix, sample_id, ".pdf")),
            width = 7 * length(cell_type_vec),
            height = 7 * (length(deconvo_tools) + include_actual)
        )
        print(plot_grid(plotlist = plot_list, ncol = length(cell_type_vec)))
        dev.off()
        
        #   For figures in the paper, create a PDF version with one plot per
        #   page. Remove the legend to match with other plots in the paper,
        #   including those with a discrete scale (vis_clus) where the legend
        #   is outside the plot
        for (i in 1:length(plot_list)) {
            plot_list[[i]] = plot_list[[i]] +
                theme(legend.position = "none")
        }
        pdf(
            file.path(
                plot_dir, paste0(pdf_prefix, sample_id, "_individual.pdf")
            )
        )
        print(plot_list)
        dev.off()
    }
}

#   Write a PDF to 'filename' of a barplot of section-wide cell-type
#   proportions against the ground-truth
#
#   prop_df: tibble with columns 'sample_id', 'deconvo_tool' (which should
#       include the ground-truth), and 'prop'. One row per sample per deconvo
#       tool per cell type.
#   filename: character relative to 'plot_dir', with extension ".pdf"
prop_barplot <- function(prop_df, filename) {
    plot_list <- list()
    for (sample_id in sample_ids) {
        plot_list[[sample_id]] <- ggplot(
            prop_df |> filter(sample_id == {{ sample_id }}),
            aes(x = deconvo_tool, y = prop, fill = cell_type)
        ) +
            geom_bar(stat = "identity") +
            labs(
                x = "Method", y = "Sample-Wide Proportion", fill = "Cell Type",
                title = sample_id
            ) +
            scale_x_discrete(labels = c("actual" = "CART-Calculated")) +
            scale_fill_manual(values = cell_type_labels) +
            theme_bw(base_size = 16)
    }
    pdf(file.path(plot_dir, filename))
    print(plot_list)
    dev.off()
}

#   Write a PDF to 'filename' of a barplot of section-wide cell-type
#   proportions against the ground-truth. Intended as a supplementary figure
#   in the manuscript
#
#   prop_df: tibble with columns 'sample_id', 'deconvo_tool' (which should
#       include the ground-truth), and 'prop'. One row per sample per deconvo
#       tool per cell type.
#   filename: character relative to 'plot_dir', with extension ".pdf"
prop_barplot_paper <- function(prop_df, filename) {
    p = ggplot(prop_df, aes(x = deconvo_tool, y = prop, fill = cell_type)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ sample_id, nrow = 1) +
        labs(x = "Method", y = "Sample-Wide Proportion", fill = "Cell Type") +
        scale_x_discrete(labels = c("actual" = "CART")) +
        scale_fill_manual(values = cell_type_labels) +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    pdf(file.path(plot_dir, filename), width = 10, height = 5)
    print(p)
    dev.off()
}

#   Print a table of KL divergences between measured cell-type proportions in
#   each section (averaged across sections) against the ground truth
kl_table <- function(full_df) {
    full_df %>%
        #   Compute cell-type proportions in each spot
        group_by(sample_id, deconvo_tool, cell_type) %>%
        summarize(observed = sum(observed), actual = sum(actual)) %>%
        group_by(sample_id, deconvo_tool) %>%
        mutate(
            observed = observed / sum(observed),
            actual = actual / sum(actual),
        ) |>
        #   Compute each term in the sum for KL divergence
        group_by(sample_id, deconvo_tool, cell_type) |>
        summarize(kl_piece = actual * log(actual / observed)) |>
        #   Add all terms to form the sum for each sample
        group_by(sample_id, deconvo_tool) |>
        summarize(kl = sum(kl_piece)) |>
        #   Take the mean across samples to form one value per tool
        group_by(deconvo_tool) |>
        summarize(kl = mean(kl))
}

#   Given a tibble 'metrics_df', write a scatterplot to PDF at 'filename'
corr_rmse_plot <- function(metrics_df, filename) {
    #   A tibble of labels for the average points on the plot
    text_df <- metrics_df |>
        group_by(deconvo_tool) |>
        summarize(Correlation = mean(Correlation), RMSE = mean(RMSE)) |>
        ungroup() |>
        mutate(
            Correlation = paste0("Avg. Cor = ", round(Correlation, 2)),
            RMSE = paste0("Avg. RMSE = ", round(RMSE, 2))
        )

    p <- ggplot(
        metrics_df,
        aes(x = Correlation, y = RMSE, color = cell_type, shape = sample_id)
    ) +
        facet_wrap(~deconvo_tool) +
        geom_point() +
        geom_text(
            data = text_df,
            aes(
                x = min(metrics_df$Correlation),
                y = 0.95 * max(metrics_df$RMSE),
                label = Correlation, color = NULL, shape = NULL
            ),
            vjust = 1, hjust = 0, show.legend = FALSE
        ) +
        geom_text(
            data = text_df,
            aes(
                x = min(metrics_df$Correlation),
                y = 0.85 * max(metrics_df$RMSE),
                label = RMSE, color = NULL, shape = NULL
            ),
            vjust = 1, hjust = 0, show.legend = FALSE
        ) +
        scale_y_continuous(
            limits = c(0, max(metrics_df$RMSE) * 1.05),
            expand = c(0, 0)
        ) +
        scale_color_manual(
            values = cell_type_labels,
            breaks = names(cell_type_labels)
        ) +
        scale_shape_manual(
            values = shape_scale,
            breaks = names(shape_scale)
        ) +
        labs(color = "Cell Type", shape = "Sample ID") +
        theme_bw(base_size = 15)

    pdf(file.path(plot_dir, filename), height = 3, width = 9)
    print(p)
    dev.off()
}

#   Given a tibble with columns 'label' (manual layer label), 'deconvo_tool',
#   'cell_type', and 'count', write a set of barplots to PDF under [plot_dir]
#   with name [filename]. 'ylab' give the y-axis label; 'x_var' is the x-axis
#   variable (as a string); 'fill_var' is the fill variable as a string;
#   'fill_scale' is passed to 'scale_fill_manual(values = [fill_scale])';
#   'fill_lab' is the fill label; 'xlab' is the x-axis label
#
#   The barplots are faceted by deconvo_tool, with x-axis including each
#   manually annotated layer. Each barplot includes counts for each cell type
#   in each layer. Each cell type is expected to have a maximal value (across
#   all bars in the facet) at a particular layer; an "O" is placed at the layer
#   with maximal value for each cell type if the layer is "correct" (e.g.
#   'Excit_L3' has maximal value in layer 3), and an "X" is placed if the layer
#   is incorrect. Total counts of "O"s for each facet across all cell types is
#   tallied and reported in the facet titles.
layer_dist_barplot <- function(counts_df, filename, ylab, x_var, fill_var, fill_scale, fill_lab, xlab) {
    #   Add a column 'layer_match' to indicate rows where each cell type
    #   has a maximal value across layers. We'll mark these with an "X"
    #   on the barplots
    counts_df <- counts_df |>
        group_by(deconvo_tool, cell_type) |>
        mutate(layer_match = count == max(count)) |>
        ungroup()

    #   Add a column 'correct_layer' indicating whether for a cell type
    #   and deconvo tool, the cell_type has maximal value in the correct/
    #   expected layer
    counts_df$correct_layer <- sapply(
        1:nrow(counts_df),
        function(i) {
            counts_df$layer_match[i] &&
                (counts_df$label[i] %in%
                    corresponding_layers[[
                        as.character(counts_df$cell_type)[i]
                    ]]
                )
        }
    )

    #   For each deconvo tool, add up how many times cell types have maximal
    #   value in the correct layers
    correct_df <- counts_df |>
        group_by(deconvo_tool) |>
        summarize(num_matches = sum(correct_layer)) |>
        ungroup()
    print("Number of times cell types have maximal value in the correct layer:")
    print(correct_df)

    print("Full list of which cell types matched the expected layer, by method:")
    counts_df |>
        group_by(deconvo_tool) |>
        filter(correct_layer) |>
        select(cell_type) |>
        ungroup() |>
        print(n = nrow(counts_df))

    #   Add the "layer accuracy" in the facet titles in the upcoming plot
    correct_labeller <- paste0(
        correct_df$deconvo_tool, ": ", correct_df$num_matches, "/",
        length(cell_types)
    )
    names(correct_labeller) <- correct_df |> pull(deconvo_tool)
    correct_labeller <- labeller(deconvo_tool = correct_labeller)

    p <- ggplot(
        counts_df,
        aes_string(x = x_var, y = "count", fill = fill_var)
    ) +
        facet_wrap(~deconvo_tool, labeller = correct_labeller) +
        geom_bar(stat = "identity") +
        labs(x = xlab, y = ylab, fill = fill_lab) +
        scale_fill_manual(values = fill_scale) +
        geom_text(
            aes(label = ifelse(correct_layer, "O", ifelse(layer_match, "X", ""))),
            position = position_stack(vjust = 0.5)
        ) +
        theme_bw(base_size = 16) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

    pdf(file.path(plot_dir, filename), width = 10, height = 5)
    print(p)
    dev.off()
}

################################################################################
#   Read in and format cell counts into a table apt for plotting with ggplot
################################################################################

sample_ids <- readLines(sample_ids_path)
added_colnames <- c("barcode", "sample_id", "deconvo_tool", "obs_type")

shape_scale <- c(16, 17, 15, 3, 8)
names(shape_scale) <- c(sample_ids, "average")

observed_df <- as_tibble(read.csv(raw_results_path))
observed_df$obs_type <- "observed"

#   Use titlecase (note this is MUCH faster than 'toTitleCase')
observed_df$deconvo_tool[observed_df$deconvo_tool == "tangram"] <- "Tangram"
observed_df$deconvo_tool[observed_df$deconvo_tool == "cell2location"] <-
    "Cell2location"

#   Plot counts for each cell type without collapsing cell categories
observed_df_long <- observed_df %>%
    melt(
        id.vars = added_colnames, variable.name = "cell_type",
        value.name = "count"
    ) %>%
    pivot_wider(
        names_from = obs_type, values_from = count,
    )

spe <- readRDS(spe_IF_in)

spatial_counts_plot_full(
    spe, observed_df_long, cell_types, FALSE, "spatial_counts_fullres_"
)

#   Gather collapsed cell counts so that each row is a unique
#   cell type, spot, sample, and deconvolution method with two values: measured
#   and ground-truth
full_df <- read.csv(collapsed_results_path) |>
    melt(
        id.vars = added_colnames, variable.name = "cell_type",
        value.name = "count"
    ) |>
    mutate(cell_type = str_to_title(cell_type)) |>
    pivot_wider(
        names_from = obs_type, values_from = count,
    )

#   Use titlecase (note this is MUCH faster than 'toTitleCase')
full_df$deconvo_tool[full_df$deconvo_tool == "tangram"] <- "Tangram"
full_df$deconvo_tool[full_df$deconvo_tool == "cell2location"] <-
    "Cell2location"

################################################################################
#   Exploratory plots
################################################################################

#   As a figure for the paper, plot the "all_spots" plot of counts for just one
#   sample and facet by cell type instead of sample
count_df <- full_df |>
    filter(sample_id == "Br6522_Ant_IF")

metrics_df <- count_df %>%
    group_by(deconvo_tool, cell_type) %>%
    summarize(
        corr = round(cor(observed, actual), 2),
        rmse = signif(mean((observed - actual)**2)**0.5, 3)
    ) %>%
    ungroup()

#   Improve labels for plotting
metrics_df$corr <- paste("Cor =", metrics_df$corr)
metrics_df$rmse <- paste("RMSE =", metrics_df$rmse)

p <- ggplot(count_df) +
    geom_point(
        aes(x = observed, y = actual, color = cell_type),
        alpha = 0.01
    ) +
    geom_abline(
        intercept = 0, slope = 1, linetype = "dashed", color = "red"
    ) +
    facet_grid(
        rows = vars(cell_type), cols = vars(deconvo_tool)
    ) +
    guides(col = guide_legend(override.aes = list(alpha = 1))) +
    labs(
        x = "Software-Estimated Counts", y = "CART-Calculated Counts",
        color = "Cell Type"
    ) +
    geom_text(
        data = metrics_df,
        mapping = aes(x =  0, y = max(count_df$actual), label = corr),
        hjust = 0, vjust = 1
    ) +
    geom_text(
        data = metrics_df,
        mapping = aes(x = 0, y = 0.85 * max(count_df$actual), label = rmse),
        hjust = 0, vjust = 1
    ) +
    scale_color_manual(values = cell_type_labels) +
    theme_bw(base_size = 15)

#   Data-specific adjustment for the sake of a quick figure
if (cell_group == "broad") {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

pdf(file.path(plot_dir, "counts_all_spots_figure.pdf"))
print(p)
dev.off()

#-------------------------------------------------------------------------------
#   Report overall performance to log
#-------------------------------------------------------------------------------

#   Consider the count of each cell type in each spot a single data point, and
#   compute correlation and RMSE between ground-truth counts and those measured
#   by each deconvolution software
print("Overall performance of deconvolution methods (accuracy of spatial variation):")
full_df |>
    #   For each sample ID and deconvo tool, compute correlation of all counts
    #   for all cell types
    group_by(deconvo_tool, sample_id) |>
    summarize(
        corr = cor(observed, actual),
        rmse = mean((observed - actual)**2)**0.5
    ) |>
    #   Now average across samples
    group_by(deconvo_tool) |>
    summarize(
        corr = round(mean(corr), 2), rmse = signif(mean(rmse), 3)
    )

print("Overall performance of deconvolution methods w/o 'other' (accuracy of spatial variation):")
full_df |>
    filter(cell_type != "Other") |>
    #   For each sample ID and deconvo tool, compute correlation of all counts
    #   for all cell types
    group_by(deconvo_tool, sample_id) |>
    summarize(
        corr = cor(observed, actual),
        rmse = mean((observed - actual)**2)**0.5
    ) |>
    #   Now average across samples
    group_by(deconvo_tool) |>
    summarize(
        corr = round(mean(corr), 2), rmse = signif(mean(rmse), 3)
    )

#-------------------------------------------------------------------------------
#   Stacked barplot of cell-type proportions
#-------------------------------------------------------------------------------

prop_df <- full_df %>%
    group_by(sample_id, deconvo_tool, cell_type) %>%
    summarize(observed = sum(observed), actual = sum(actual)) %>%
    group_by(sample_id, deconvo_tool) %>%
    mutate(
        observed = observed / sum(observed),
        actual = actual / sum(actual),
    ) %>%
    ungroup() |>
    pivot_longer(
        cols = c("observed", "actual"), values_to = "prop", names_to = "source"
    )

temp_actual <- prop_df |>
    filter(source == "actual", deconvo_tool == "Tangram") |>
    mutate(deconvo_tool = "actual")
temp_observed <- prop_df |>
    filter(source == "observed")
prop_df <- rbind(temp_actual, temp_observed)

prop_barplot(prop_df, "prop_barplots.pdf")
prop_barplot(
    prop_df |> filter(cell_type != "Other"), "prop_barplots_no_other.pdf"
)
prop_barplot_paper(
    prop_df |> filter(cell_type != "Other"), "prop_barplots_no_other_paper.pdf"
)

print("Accuracy of overall cell-type proportions per section (Mean KL divergence from ground-truth):")
print(kl_table(full_df))
print("Accuracy of overall cell-type proportions per section w/o 'other' (Mean KL divergence from ground-truth):")
print(kl_table(full_df |> filter(cell_type != "Other")))

#-------------------------------------------------------------------------------
#   Plot distribution of correlation & RMSE by sample and deconvo tool
#-------------------------------------------------------------------------------

metrics_df <- full_df |>
    group_by(deconvo_tool, sample_id, cell_type) |>
    summarize(
        Correlation = cor(observed, actual),
        RMSE = mean((observed - actual)**2)**0.5
    ) |>
    ungroup()

corr_rmse_plot(metrics_df, "corr_RMSE_scatter.pdf")
corr_rmse_plot(
    metrics_df |> filter(cell_type != "Other"),
    "corr_RMSE_scatter_no_other.pdf"
)

#-------------------------------------------------------------------------------
#   Plot total counts per sample for tangram and cell2location
#-------------------------------------------------------------------------------

count_df <- full_df |>
    filter(deconvo_tool %in% c("Tangram", "Cell2location")) |>
    group_by(sample_id, deconvo_tool) |>
    summarize(observed = sum(observed), actual = sum(actual)) |>
    mutate(diff = abs((observed - actual) / (observed + actual)))

pdf(file.path(plot_dir, "total_cells_sample.pdf"))
ggplot(count_df, aes(x = deconvo_tool, y = diff, color = sample_id)) +
    geom_point() +
    theme_bw(base_size = 10) +
    scale_y_continuous(limits = c(0, max(count_df$diff) * 1.1), expand = c(0, 0)) +
    labs(x = "Deconvolution tool", y = "Absolute proportion difference")
dev.off()

#-------------------------------------------------------------------------------
#   Plot total counts per spot for tangram and cell2location
#-------------------------------------------------------------------------------

#   For each spot, plot the provided vs. computed total number of cells for the
#   deconvolution methods that aren't constrained to use the same totals as
#   provided (tangram and cell2location)
count_df <- full_df %>%
    filter(deconvo_tool %in% c("Tangram", "Cell2location")) %>%
    group_by(barcode, sample_id, deconvo_tool) %>%
    summarize(observed = sum(observed), actual = sum(actual)) %>%
    ungroup()

#   Compute metrics for each deconvolution tool: correlation between
#   observed and actual values as well as RMSE
metrics_df <- count_df %>%
    group_by(deconvo_tool) %>%
    summarize(
        corr = round(cor(observed, actual), 2),
        rmse = signif(mean((observed - actual)**2)**0.5, 3)
    ) %>%
    ungroup()

#   Improve labels for plotting
metrics_df$corr <- paste("Cor =", metrics_df$corr)
metrics_df$rmse <- paste("RMSE =", metrics_df$rmse)

pdf(file.path(plot_dir, "total_cells_spot.pdf"))
ggplot(count_df) +
    geom_point(aes(x = observed, y = actual), alpha = 0.01) +
    facet_wrap(~deconvo_tool) +
    coord_fixed() +
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
dev.off()

#-------------------------------------------------------------------------------
#   Counts: "all" and "across"
#-------------------------------------------------------------------------------

#   For the "all_spots" plot, also write a CSV (to be used as a suppl. table)
#   of metrics
metrics_df <- full_df %>%
    group_by(deconvo_tool, sample_id, cell_type) %>%
    summarize(
        corr = cor(observed, actual),
        rmse = mean((observed - actual)**2)**0.5
    ) %>%
    ungroup()

write.csv(
    metrics_df,
    file.path(
        processed_dir, paste0("spatial_cell_type_metrics_", cell_group, ".csv")
    ),
    row.names = FALSE, quote = FALSE
)

#   Plot cell-type counts for all spots
all_spots(full_df, "counts_all_spots_scatter.pdf")

#   Plot cell-type counts summed across spots
count_df <- full_df %>%
    group_by(sample_id, deconvo_tool, cell_type) %>%
    summarize(observed = sum(observed), actual = sum(actual)) %>%
    ungroup()

across_spots(count_df, "counts_across_spots_scatter.pdf", x_angle = 90)
across_spots(
    count_df |> filter(cell_type != "Other"),
    "counts_across_spots_scatter_no_other.pdf",
    x_angle = 90
)

#-------------------------------------------------------------------------------
#   Proportions: "all" and "across"
#-------------------------------------------------------------------------------

#   Plot cell-type proportions for all spots

prop_df <- full_df %>%
    group_by(barcode, sample_id, deconvo_tool) %>%
    mutate(
        observed = observed / sum(observed),
        actual = actual / sum(actual),
    ) %>%
    filter(
        !is.na(observed) & !is.na(actual)
    ) %>%
    ungroup()
all_spots(prop_df, "props_all_spots_scatter.pdf")

#   First, sum counts for each (cell type, sample, deconvo tool) across
#   spots. Then assign a cell-type proportion by dividing by total cells
#   for each (sample, deconvo tool). Note this method yields many more non-NA
#   values since taking proportions in each spot and averaging across spots
#   introduces NAs whenever either the observed or actual total cell count is 0
#   for a spot.
prop_df <- full_df %>%
    group_by(sample_id, deconvo_tool, cell_type) %>%
    summarize(observed = sum(observed), actual = sum(actual)) %>%
    group_by(sample_id, deconvo_tool) %>%
    mutate(
        observed = observed / sum(observed),
        actual = actual / sum(actual),
    ) %>%
    ungroup()

#   Plot versions with and without the "Other" cell type
across_spots(prop_df, "props_across_spots_scatter.pdf")
across_spots(
    prop_df |> filter(cell_type != "Other"),
    "props_across_spots_scatter_no_other.pdf"
)

#-------------------------------------------------------------------------------
#   Adjusted counts: "all" and "across"
#-------------------------------------------------------------------------------

#   Counts, where we take the estimated proportion and multiply by the
#   "ground-truth" total count
count_df <- full_df %>%
    group_by(barcode, sample_id, deconvo_tool) %>%
    mutate(
        observed = sum(actual) * observed / sum(observed)
    ) %>%
    ungroup()
count_df$observed[is.na(count_df$observed)] <- 0
all_spots(count_df, "adjusted_counts_all_spots_scatter.pdf")

count_df <- full_df %>%
    group_by(sample_id, deconvo_tool, cell_type) %>%
    summarize(observed = sum(observed), actual = sum(actual)) %>%
    group_by(sample_id, deconvo_tool) %>%
    mutate(observed = sum(actual) * observed / sum(observed)) %>%
    ungroup()

#   Plot versions with and without the "Other" cell type
across_spots(count_df, "adjusted_counts_across_spots_scatter.pdf", x_angle = 90)
across_spots(
    count_df |> filter(cell_type != "Other"),
    "adjusted_counts_across_spots_scatter_no_other.pdf",
    x_angle = 90
)

#-------------------------------------------------------------------------------
#   Spatial distribution of counts for each sample, cell_type, deconvo method
#-------------------------------------------------------------------------------

spatial_counts_plot_full(
    spe, full_df, cell_types_actual, TRUE, "spatial_counts_"
)

#-------------------------------------------------------------------------------
#   Spatial distribution of marker gene expression vs. cell-type counts
#-------------------------------------------------------------------------------

#   Read in markers and marker stats
marker_stats <- readRDS(marker_object_in)
markers <- readLines(marker_in)

#   Subset marker_stats to include just the rows where each gene is a marker for
#   the given cell type
marker_stats <- marker_stats %>%
    group_by(gene) %>%
    filter(
        ratio == max(ratio),
        gene %in% markers
    ) %>%
    ungroup()

#   Collapse cell types to the resolution detectable on the IF images
#   ("ground-truth"). Note there are some subtle flaws that emerge from this
#   approach, which we'll ignore for the sake of an easy approximate
#   visualization: a marker might start with a good mean ratio between 'Excit'
#   and 'Inhib', but we collapse these categories, changing the real mean ratio
#   and thus rank of the marker (not accounted for here).
marker_stats$cellType.target <-as.character(marker_stats$cellType.target)
marker_stats$cellType.target[
    grep("^(Excit|Inhib)", marker_stats$cellType.target)
] <- "Neuron"
marker_stats$cellType.target[marker_stats$cellType.target == "OPC"] <- "Oligo"
marker_stats$cellType.target[
    marker_stats$cellType.target == "EndoMural"
] <- "Other"
marker_stats$cellType.target <- as.factor(marker_stats$cellType.target)

stopifnot(
    all(sort(unique(marker_stats$cellType.target)) == sort(cell_types_actual))
)

#   Add mean marker-gene expression for each associated cell type to the main
#   data frame with cell-type counts
full_df$marker_express <- NA
for (cell_type in cell_types_actual) {
    #   Grab markers for this cell type
    these_markers <- marker_stats %>%
        filter(cellType.target == cell_type) %>%
        arrange(desc(ratio)) %>%
        head(n = 25) %>%
        pull(gene)

    stopifnot(all(these_markers %in% rownames(spe)))

    for (sample_id in sample_ids) {
        spe_small <- spe[these_markers, spe$sample_id == sample_id]

        for (deconvo_tool in deconvo_tools) {
            #   Form a temporary tibble with mean expression across markers for
            #   each spot
            temp_df <- as_tibble(
                data.frame(
                    "marker_express_temp" = colMeans(assays(spe_small)$counts),
                    "barcode" = colnames(spe_small),
                    "sample_id" = sample_id,
                    "cell_type" = cell_type,
                    "deconvo_tool" = deconvo_tool
                )
            )

            #   Add this quantity to the full tibble that has cell-type counts
            temp_df <- left_join(
                full_df, temp_df,
                by = c("barcode", "sample_id", "cell_type", "deconvo_tool")
            )
            new_indices <- !is.na(temp_df$marker_express_temp)
            full_df[new_indices, "marker_express"] <- temp_df[
                new_indices, "marker_express_temp"
            ]
        }
    }
}

#   All rows in full_df should've been filled with an expression value
stopifnot(all(!is.na(full_df$marker_express)))

#   Compute correlation for each deconvolution tool
metrics_df <- full_df %>%
    group_by(deconvo_tool, sample_id, cell_type) %>%
    summarize(
        corr = round(cor(observed, marker_express), 2),
    ) %>%
    ungroup()

metrics_df$corr <- paste("Cor =", metrics_df$corr)

#   For each cell type, plot counts of this cell type, measured by each
#   deconvolution tool, against the mean marker-gene expression for that cell
#   type (scatterplot)
plot_list <- lapply(
    cell_types_actual,
    function(ct) {
        full_df_small <- full_df %>%
            filter(cell_type == ct)

        p <- ggplot(full_df_small) +
            geom_point(
                aes(x = observed, y = marker_express, color = sample_id),
                alpha = 0.01
            ) +
            facet_grid(
                rows = vars(sample_id), cols = vars(deconvo_tool)
            ) +
            guides(col = guide_legend(override.aes = list(alpha = 1))) +
            geom_text(
                data = metrics_df %>% filter(cell_type == ct),
                mapping = aes(
                    x = max(full_df_small$observed),
                    y = max(full_df_small$marker_express) / 7,
                    label = corr
                ),
                hjust = 1
            ) +
            labs(
                title = ct,
                x = paste0(
                    "Software-Estimated ", ct, " Counts"
                ),
                y = "Mean Marker-Gene Counts",
                color = "Sample ID"
            ) +
            theme_bw(base_size = 10)

        return(p)
    }
)

pdf(file.path(plot_dir, "marker_vs_ct_counts.pdf"))
print(plot_list)
dev.off()

#   Now plot a version suitable as a supplementary figure: this time, one page
#   per sample and facet by cell type
plot_list <- lapply(
    sample_ids,
    function(this_sample_id) {
        full_df_small <- full_df %>%
            filter(sample_id == this_sample_id)
        
        p <- ggplot(full_df_small) +
            geom_point(
                aes(x = observed, y = marker_express, color = cell_type),
                alpha = 0.01
            ) +
            facet_grid( rows = vars(cell_type), cols = vars(deconvo_tool)) +
            guides(col = guide_legend(override.aes = list(alpha = 1))) +
            geom_text(
                data = metrics_df %>% filter(sample_id == this_sample_id),
                mapping = aes(
                    x = max(full_df_small$observed),
                    y = max(full_df_small$marker_express) / 7,
                    label = corr
                ),
                hjust = 1
            ) +
            labs(
                title = this_sample_id,
                x = paste0(
                    "Software-Estimated Counts"
                ),
                y = "Mean Marker-Gene Counts",
                color = "Cell Type"
            ) +
            theme_bw(base_size = 15)
        
        return(p)
    }
)

pdf(file.path(plot_dir, "marker_vs_ct_counts_paper.pdf"))
print(plot_list)
dev.off()

#-------------------------------------------------------------------------------
#   Spatial distribution of cell-types compared against manual layer annotation
#   (boxplots)
#-------------------------------------------------------------------------------

#   Read layer annotation in for each sample and match to a barcode
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
observed_df_long <- left_join(
    observed_df_long, layer_ann,
    by = c("barcode", "sample_id")
)

#   Clean up labels
observed_df_long$label <- tolower(observed_df_long$label)
observed_df_long$label <- sub("layer", "Layer ", observed_df_long$label)
observed_df_long$label[observed_df_long$label == "wm"] <- "WM"
stopifnot(
    all(unlist(corresponding_layers) %in% unique(observed_df_long$label))
)

#   Average counts of each cell type in each layer as annotated; filter NA
#   labels (there are 2 intentional NAs where spots should be dropped)
counts_df <- observed_df_long |>
    filter(!is.na(label)) |>
    group_by(label, deconvo_tool, sample_id, cell_type) |>
    summarize(count = mean(observed)) |>
    ungroup()

#   Create plots for each cell type
plot_list <- list()
for (cell_type in cell_types) {
    y_max <- counts_df |>
        filter(cell_type == {{ cell_type }}) |>
        summarize(y_max = max(count)) |>
        pull(y_max)

    plot_list[[cell_type]] <- ggplot(
        counts_df |> filter(cell_type == {{ cell_type }}),
        aes(x = label, y = count, color = deconvo_tool)
    ) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.05) +
        labs(
            x = "Annotated Layer",
            y = "Average Predicted Count",
            color = "Deconvolution Tool"
        ) +
        #   Facet purely for aesthetic purposes: there is only one cell type
        facet_wrap(~cell_type) +
        theme_bw(base_size = 23) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        coord_cartesian(ylim = c(0, y_max)) +
        scale_y_continuous(expand = c(0, 0, 0, 0.05))
}

pdf(
    file.path(plot_dir, "layer_distribution.pdf"),
    width = 10
)
print(plot_list)
dev.off()

#-------------------------------------------------------------------------------
#   Spatial distribution of cell-types compared against manual layer annotation
#   (barplots- raw counts)
#-------------------------------------------------------------------------------

#   Read in SCE just to get colors for software-estimated cell types
sce <- readRDS(sce_in)
estimated_cell_labels <- metadata(sce)[[paste0("cell_type_colors_", cell_group)]]
names(estimated_cell_labels) <- gsub("/", "_", names(estimated_cell_labels))

counts_df <- observed_df_long |>
    filter(!is.na(label)) |>
    #   Average counts within sample for a given layer/cell type/ deconvo tool
    group_by(label, deconvo_tool, sample_id, cell_type) |>
    summarize(count = mean(observed)) |>
    #   Average these averages across sample
    group_by(label, deconvo_tool, cell_type) |>
    summarize(count = sum(count) / length(sample_ids)) |>
    ungroup()

#   Meaning of an example section of one bar in one facet of this plot:
#   Orange bar at white matter for cell2location: of all spots manually
#   annotated as white matter, the size of the bar section is the across-sample
#   mean of average oligo counts in those spots for each sample. Total bar
#   heights vary within deconvo tool because cell-count density varies between
#   layers. Total bar heights may vary between deconvo tools when total counts
#   per spot aren't preserved (e.g. C2L doesn't match the cell counts you
#   provide it)
layer_dist_barplot(
    counts_df,
    filename = "layer_distribution_barplot_raw.pdf",
    xlab = "Annotated Layer", ylab = "Average Predicted Count", x_var = "label",
    fill_lab = "Cell Type", fill_var = "cell_type",
    fill_scale = estimated_cell_labels
)

#   Also write a copy where will switch the x-axis and fill
layer_dist_barplot(
    counts_df,
    filename = "layer_distribution_barplot_raw_inverted.pdf",
    xlab = "Cell Type", ylab = "Average Predicted Count", x_var = "cell_type",
    fill_lab = "Annotated Layer", fill_var = "label",
    fill_scale = libd_layer_colors
)

#-------------------------------------------------------------------------------
#   Spatial distribution of cell-types compared against manual layer annotation
#   (barplots- proportions by layer)
#-------------------------------------------------------------------------------

#   Add all spots, regardless of sample
counts_df <- observed_df_long |>
    filter(!is.na(label)) |>
    #   For each manually annotated label and deconvo tool, normalize by the
    #   total counts of all cell types and samples
    group_by(deconvo_tool, label) |>
    mutate(observed = observed / sum(observed)) |>
    #   Now for each label, deconvo tool and cell type, add up counts for all
    #   samples and relevant spots
    group_by(deconvo_tool, label, cell_type) |>
    summarize(count = sum(observed))

layer_dist_barplot(
    counts_df,
    filename = "layer_distribution_barplot_prop_uneven.pdf",
    xlab = "Annotated Layer", ylab = "Proportion of Counts", x_var = "label",
    fill_lab = "Cell Type", fill_var = "cell_type",
    fill_scale = estimated_cell_labels
)

#   Also write a copy where will switch the x-axis and fill
layer_dist_barplot(
    counts_df,
    filename = "layer_distribution_barplot_prop_uneven_inverted.pdf",
    xlab = "Cell Type", ylab = "Proportion of Counts", x_var = "cell_type",
    fill_lab = "Annotated Layer", fill_var = "label",
    fill_scale = libd_layer_colors
)

#   Weight each sample equally
counts_df <- observed_df_long |>
    filter(!is.na(label)) |>
    #   For each manually annotated label, deconvo tool and sample_id, normalize
    #   by the total counts of all cell types
    group_by(deconvo_tool, label, sample_id) |>
    mutate(observed = observed / sum(observed)) |>
    #   Now for each label, deconvo tool, sample_id and cell type, add up counts
    #   for all relevant spots
    group_by(deconvo_tool, label, cell_type, sample_id) |>
    summarize(count = sum(observed)) |>
    #   Now average across samples
    group_by(deconvo_tool, label, cell_type) |>
    summarize(count = mean(count)) |>
    ungroup()

layer_dist_barplot(
    counts_df,
    filename = "layer_distribution_barplot_prop_even.pdf",
    xlab = "Annotated Layer", ylab = "Proportion of Counts", x_var = "label",
    fill_lab = "Cell Type", fill_var = "cell_type",
    fill_scale = estimated_cell_labels
)

#   Also write a copy where will switch the x-axis and fill
layer_dist_barplot(
    counts_df,
    filename = "layer_distribution_barplot_prop_even_inverted.pdf",
    xlab = "Cell Type", ylab = "Proportion of Counts", x_var = "cell_type",
    fill_lab = "Annotated Layer", fill_var = "label",
    fill_scale = libd_layer_colors
)
#-------------------------------------------------------------------------------
#   Spatial distribution of cell-types compared against manual layer annotation
#   (barplots- proportions by layer). Inverted so that x-axis is cell type
#-------------------------------------------------------------------------------

counts_df <- observed_df_long |>
    filter(!is.na(label)) |>
    #   Normalize counts by the sum within sample ID, cell type, and deconvo
    #   tool (each sample contributes equally)
    group_by(sample_id, deconvo_tool, cell_type) |>
    mutate(observed = observed / sum(observed)) |>
    #   Add up counts for all spots in each layer for each sample ID, cell
    #   type, and deconvo tool
    group_by(deconvo_tool, label, sample_id, cell_type) |>
    summarize(count = sum(observed)) |>
    #   Now average counts across samples (Note 'sum(count) / 4' must be used
    #   and not 'mean(count)', since not all samples have all layers
    group_by(deconvo_tool, label, cell_type) |>
    summarize(count = sum(count) / length(sample_ids)) |>
    ungroup()

#   Interpretation of this plot:
#   For a given deconvo tool, examine one cell type. The size of each component
#   of the bar can be interpreted as the probability a randomly selected cell of
#   this cell type belongs in this particular layer. Note that this does not
#   normalize for layer size (i.e. if more spots belong to layer 3 than white
#   matter, the maximal layer for a given cell type is more likely to be layer
#   3)!
layer_dist_barplot(
    counts_df,
    filename = "layer_distribution_barplot_probability_inverted.pdf",
    xlab = "Cell Type", ylab = "Proportion of Counts", x_var = "cell_type",
    fill_lab = "Annotated Layer", fill_var = "label",
    fill_scale = libd_layer_colors
)

#-------------------------------------------------------------------------------
#   Spatial plots of manual layer annotation
#-------------------------------------------------------------------------------

plot_list <- list()
for (sample_id in sample_ids) {
    spe_small <- spe[, spe$sample_id == sample_id]
    counts_df <- observed_df_long |> filter(sample_id == {{ sample_id }})

    spe_small$manual_layer <- counts_df[
        match(colnames(spe_small), counts_df$barcode), "label"
    ] |> pull(label)

    #   Verify no NA labels are present, except for sample 'Br6432_Ant_IF',
    #   where 2 NA spots are expected: these should be dropped for these plots
    spe_small$manual_layer[spe_small$manual_layer == ""] <- NA
    if (sample_id != "Br6432_Ant_IF" & any(is.na(spe_small$manual_layer))) {
        stop("No spot labels should be NA")
    } else if (sample_id == "Br6432_Ant_IF") {
        stopifnot(length(which(is.na(spe_small$manual_layer))) == 2)
        spe_small <- spe_small[, !is.na(spe_small$manual_layer)]
    }
    
    plot_list[[sample_id]] <- spot_plot(
        spe_small, sample_id = sample_id, title = sample_id,
        var_name = "manual_layer", colors = libd_layer_colors,
        is_discrete = TRUE, include_legend = FALSE
    )
}

pdf(file.path(plot_dir, "spot_layer_labels.pdf"))
print(plot_list)
dev.off()
