library("here")
library("ggplot2")
library("jaffelab")
library("tidyverse")
library("reshape2")
library("spatialLIBD")
library("cowplot")

cell_group <- "broad" # "broad" or "layer"

sample_ids_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "sample_ids.txt"
)

deconvo_tools = c('tangram', 'cell2location', 'SPOTlight')

plot_dir <- here(
    "plots", "spot_deconvo", "05-shared_utilities", cell_group
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

cell_types_actual <- c("astro", "micro", "neuron", "oligo", "other")
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
    
    #   For excitatory layers, make a list of the layers contained
    corresponding_layers = list(
        "Excit_L2_3" = c("layer2", "layer3"),
        "Excit_L3" = c("layer3"),
        "Excit_L3_4_5" = c("layer3", "layer4", "layer5"),
        "Excit_L4" = c("layer4"),
        "Excit_L5" = c("layer5"),
        "Excit_L5_6" = c("layer5", "layer6"),
        "Excit_L6" = c("layer6")
    )
}

cell_type_labels = c("#3BB273", "#663894", "#E49AB0", "#E07000", "#95B8D1")
names(cell_type_labels) = cell_types_actual

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

            p = ggplot(count_df_small) +
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
                    x = "Software-estimated",
                    y = "CART-calculated",
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
    # return(plot_list)
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

    #   Improve labels for plotting
    metrics_df$corr <- paste("Cor =", metrics_df$corr)
    metrics_df$rmse <- paste("RMSE =", metrics_df$rmse)

    p <- ggplot(count_df) +
        geom_point(
            aes(x = observed, y = actual, shape = sample_id, color = cell_type)
        ) +
        facet_wrap(~deconvo_tool) +
        geom_abline(
            intercept = 0, slope = 1, linetype = "dashed", color = "red"
        ) +
        geom_text(
            data = metrics_df,
            mapping = aes(
                x = max(count_df$observed),
                y = min(count_df$actual),
                label = corr
            ),
            hjust = 1, vjust = 0
        ) +
        geom_text(
            data = metrics_df,
            mapping = aes(
                x = max(count_df$observed),
                y =  0.15 * max(count_df$actual),
                label = rmse
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
            x = "Software-estimated", y = "CART-calculated",
            color = "Cell Type", shape = "Sample ID"
        ) +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = x_angle))

    pdf(file.path(plot_dir, plot_name), height = 4, width = 10)
    print(p)
    dev.off()
    # return(p)
}

#   Return a length-2 list containing a 'vis_grid_gene' plot and maximum value
#   in the plot given the SpatialExperiment, long-format
#   tibble of cell-type counts, the target sample ID, deconvo tool, and cell
#   type, column name of 'full_df' ('observed' or 'actual'), and a plot title
spatial_counts_plot = function(
        spe_small, full_df, sample_id1, deconvo_tool1, cell_type1, c_name,
        title
) {
    #   Grab counts for just this sample, deconvo tool, and cell type
    counts_df = full_df %>%
        filter(
            sample_id == sample_id1,
            deconvo_tool == deconvo_tool1,
            cell_type == cell_type1
        )
    
    #   Add counts as a column in colData(spe_small)
    spe_small$temp_ct_counts = counts_df[[c_name]][
        match(colnames(spe_small), counts_df$barcode)
    ]
    
    #   Plot spatial distribution
    p = vis_grid_gene(
        spe_small, geneid = 'temp_ct_counts', return_plots = TRUE,
        spatial = FALSE
    )[[1]] +
        labs(title = title)
    
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
spatial_counts_plot_full = function(
        spe, full_df, cell_type_vec, include_actual, pdf_prefix
    ) {
    for (sample_id in sample_ids) {
        spe_small = spe[,spe$sample_id == sample_id]
        
        i = 1
        plot_list = list()
        max_list = list()
        
        #   For each deconvo tool, make a row of plots (each including all cell
        #   types) showed the observed distribution
        for (deconvo_tool in deconvo_tools) {
            for (cell_type in cell_type_vec) {
                temp = spatial_counts_plot(
                    spe_small, full_df, sample_id, deconvo_tool, cell_type,
                    'observed',
                    paste0(cell_type, ' counts\n(', deconvo_tool, ')')
                )
                plot_list[[i]] = temp[[1]]
                max_list[[i]] = temp[[2]]
                i = i + 1
            }
        }
        
        #   Add a row showing the ground-truth counts for each cell type
        if (include_actual) {
            for (cell_type in cell_types_actual) {
                temp = spatial_counts_plot(
                    spe_small, full_df, sample_id, deconvo_tool, cell_type, 'actual',
                    paste0(cell_type, ' counts\n(Ground-truth)')
                )
                plot_list[[i]] = temp[[1]]
                max_list[[i]] = temp[[2]]
                i = i + 1
            }
        }
        
        max_mat = matrix(
            unlist(max_list), ncol = length(cell_type_vec), byrow= TRUE
        )
        
        #   Now loop back through the plot list (which will be displayed in 2D)
        #   and overwrite the scale to go as high as the largest value in the
        #   column. This allows for easy comparison between deconvo tools
        #   (and optionally the ground truth)
        for (i_col in 1:length(cell_type_vec)) {
            for (i_row in 1:(length(deconvo_tools) + include_actual)) {
                index = (i_row - 1) * length(cell_type_vec) + i_col
                upper_limit = max(max_mat[,i_col])
                
                plot_list[[index]] = plot_list[[index]] +
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
            file.path(plot_dir, paste0(pdf_prefix, sample_id, '.pdf')),
            width = 7 * length(cell_type_vec),
            height = 7 * (length(deconvo_tools) + include_actual)
        )
        print(plot_grid(plotlist = plot_list, ncol = length(cell_type_vec)))
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
prop_barplot = function(prop_df, filename) {
    plot_list = list()
    for (sample_id in sample_ids) {
        plot_list[[sample_id]] = ggplot(
            prop_df |> filter(sample_id == {{ sample_id }}),
            aes(x = deconvo_tool, y = prop, fill = cell_type)
        ) +
            geom_bar(stat = "identity") +
            labs(
                x = "Method", y = "Sample-Wide Proportion", fill = "Cell Type",
                title = sample_id
            ) +
            scale_x_discrete(labels = c("actual" = "Ground-Truth")) +
            scale_fill_manual(values = cell_type_labels) +
            theme_bw(base_size = 16)
    }
    pdf(file.path(plot_dir, filename))
    print(plot_list)
    dev.off()
}

#   Print a table of KL divergences between measured cell-type proportions in
#   each section (averaged across sections) against the ground truth
kl_table = function(full_df) {
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
        summarize(kl_piece = observed * log(observed / actual)) |>
        #   Add all terms to form the sum for each sample
        group_by(sample_id, deconvo_tool) |>
        summarize(kl = sum(kl_piece)) |>
        #   Take the mean across samples to form one value per tool
        group_by(deconvo_tool) |>
        summarize(kl = mean(kl))
}

#   Given a tibble 'metrics_df', write a scatterplot to PDF at 'filename'
corr_rmse_plot = function(metrics_df, filename) {
    p = ggplot(
        metrics_df,
        aes(x = Correlation, y = RMSE, color = cell_type, shape = sample_id)
    ) +
        facet_wrap(~deconvo_tool) +
        geom_point() +
        scale_color_manual(values = cell_type_labels) +
        labs(color = "Cell Type", shape = "Sample ID") +
        theme_bw(base_size = 13)
    
    pdf(file.path(plot_dir, filename), height = 4, width = 9)
    print(p)
    dev.off()
}

################################################################################
#   Read in and format cell counts into a table apt for plotting with ggplot
################################################################################

sample_ids <- readLines(sample_ids_path)
added_colnames <- c("barcode", "sample_id", "deconvo_tool", "obs_type")

observed_df <- as_tibble(read.csv(raw_results_path))
observed_df$obs_type = "observed"

#   Plot counts for each cell type without collapsing cell categories
observed_df_long <- observed_df %>%
    melt(
        id.vars = added_colnames, variable.name = "cell_type",
        value.name = "count"
    ) %>%
    pivot_wider(
        names_from = obs_type, values_from = count,
    )

spe = readRDS(spe_IF_in)

spatial_counts_plot_full(
    spe, observed_df_long, cell_types, FALSE, 'spatial_counts_fullres_'
)

#   Gather collapsed cell counts so that each row is a unique
#   cell type, spot, sample, and deconvolution method with two values: measured
#   and ground-truth
full_df <- read.csv(collapsed_results_path) %>%
    melt(
        id.vars = added_colnames, variable.name = "cell_type",
        value.name = "count"
    ) %>%
    pivot_wider(
        names_from = obs_type, values_from = count,
    )

################################################################################
#   Exploratory plots
################################################################################

#   Consider the count of each cell type in each spot a single data point, and
#   compute correlation and RMSE between ground-truth counts and those measured
#   by each deconvolution software
print("Overall performance of deconvolution methods (accuracy of spatial variation):")
full_df |>
    group_by(deconvo_tool) |>
    summarize(
        corr = round(cor(observed, actual), 2),
        rmse = signif(mean((observed - actual)**2)**0.5, 3)
    )

print("Overall performance of deconvolution methods w/o 'other' (accuracy of spatial variation):")
full_df |>
    group_by(deconvo_tool) |>
    filter(cell_type != "other") |>
    summarize(
        corr = round(cor(observed, actual), 2),
        rmse = signif(mean((observed - actual)**2)**0.5, 3)
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

temp_actual = prop_df |>
    filter(source == "actual", deconvo_tool == "tangram") |>
    mutate(deconvo_tool = "actual")
temp_observed = prop_df |>
    filter(source == "observed")
prop_df = rbind(temp_actual, temp_observed)

prop_barplot(prop_df, 'prop_barplots.pdf')
prop_barplot(
    prop_df |> filter(cell_type != "other"), 'prop_barplots_no_other.pdf'
)

print("Accuracy of overall cell-type proportions per section (Mean KL divergence from ground-truth):")
print(kl_table(full_df))
print("Accuracy of overall cell-type proportions per section w/o 'other' (Mean KL divergence from ground-truth):")
print(kl_table(full_df |> filter(cell_type != "other")))

#-------------------------------------------------------------------------------
#   Plot distribution of correlation & RMSE by sample and deconvo tool
#-------------------------------------------------------------------------------

metrics_df <- full_df |>
    group_by(deconvo_tool, sample_id, cell_type) |>
    summarize(
        Correlation = round(cor(observed, actual), 2),
        RMSE = signif(mean((observed - actual)**2)**0.5, 3)
    ) |>
    ungroup()

corr_rmse_plot(metrics_df, 'corr_RMSE_scatter.pdf')
corr_rmse_plot(
    metrics_df |> filter(cell_type != "other"),
    'corr_RMSE_scatter_no_other.pdf'
)

#-------------------------------------------------------------------------------
#   Plot total counts per sample for tangram and cell2location
#-------------------------------------------------------------------------------

count_df = full_df |>
    filter(deconvo_tool %in% c("tangram", "cell2location")) |>
    group_by(sample_id, deconvo_tool) |>
    summarize(observed = sum(observed), actual = sum(actual)) |>
    mutate(diff = abs((observed - actual) / (observed + actual)))

pdf(file.path(plot_dir, 'total_cells_sample.pdf'))
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
    filter(deconvo_tool %in% c("tangram", "cell2location")) %>%
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

pdf(file.path(plot_dir, 'total_cells_spot.pdf'))
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
        x = "Calculated cell count",
        y = "Provided cell count (cellpose)",
        title = "Provided vs. calculated total cells per spot"
    ) +
    theme_bw(base_size = 15)
dev.off()

#-------------------------------------------------------------------------------
#   Counts: "all" and "across"
#-------------------------------------------------------------------------------

#   Plot cell-type counts for all spots
all_spots(full_df, "counts_all_spots_scatter.pdf")

#   Plot cell-type counts summed across spots
count_df <- full_df %>%
    group_by(sample_id, deconvo_tool, cell_type) %>%
    summarize(observed = sum(observed), actual = sum(actual)) %>%
    ungroup()

across_spots(count_df, "counts_across_spots_scatter.pdf", x_angle = 90)
across_spots(
    count_df |> filter(cell_type != "other"),
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

#   Plot versions with and without the "other" cell type
across_spots(prop_df, "props_across_spots_scatter.pdf")
across_spots(
    prop_df |> filter(cell_type != "other"),
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

#   Plot versions with and without the "other" cell type
across_spots(count_df, "adjusted_counts_across_spots_scatter.pdf", x_angle = 90)
across_spots(
    count_df |> filter(cell_type != "other"),
    "adjusted_counts_across_spots_scatter_no_other.pdf",
    x_angle = 90
)

#-------------------------------------------------------------------------------
#   Spatial distribution of counts for each sample, cell_type, deconvo method
#-------------------------------------------------------------------------------

spatial_counts_plot_full(
    spe, full_df, cell_types_actual, TRUE, 'spatial_counts_'
)

#-------------------------------------------------------------------------------
#   Spatial distribution of marker gene expression vs. cell-type counts
#-------------------------------------------------------------------------------

#   Read in markers and marker stats
marker_stats = readRDS(marker_object_in)
markers = readLines(marker_in)

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
marker_stats$cellType.target = tolower(
    as.character(marker_stats$cellType.target)
)
marker_stats$cellType.target[
    grep('^(excit|inhib)', marker_stats$cellType.target)
] = "neuron"
marker_stats$cellType.target[marker_stats$cellType.target == 'opc'] = "oligo"
marker_stats$cellType.target[
    marker_stats$cellType.target == 'endomural'
] = "other"
marker_stats$cellType.target = as.factor(marker_stats$cellType.target)

stopifnot(
    all(sort(unique(marker_stats$cellType.target)) == sort(cell_types_actual))
)

#   Add mean marker-gene expression for each associated cell type to the main
#   data frame with cell-type counts
full_df$marker_express = NA
for (cell_type in cell_types_actual) {
    #   Grab markers for this cell type
    these_markers = marker_stats %>%
        filter(cellType.target == cell_type) %>%
        arrange(desc(ratio)) %>%
        head(n = 25) %>%
        pull(gene)
    
    stopifnot(all(these_markers %in% rownames(spe)))
    
    for (sample_id in sample_ids) {
        spe_small = spe[these_markers, spe$sample_id == sample_id]
        
        for (deconvo_tool in deconvo_tools) {
            #   Form a temporary tibble with mean expression across markers for
            #   each spot
            temp_df = as_tibble(
                data.frame(
                    'marker_express_temp' = colMeans(assays(spe_small)$counts),
                    'barcode' = colnames(spe_small),
                    'sample_id' = sample_id,
                    'cell_type' = cell_type,
                    'deconvo_tool' = deconvo_tool
                )
            )
            
            #   Add this quantity to the full tibble that has cell-type counts
            temp_df = left_join(
                full_df, temp_df,
                by = c("barcode", "sample_id", "cell_type", "deconvo_tool")
            )
            new_indices = !is.na(temp_df$marker_express_temp)
            full_df[new_indices, 'marker_express'] = temp_df[
                new_indices, 'marker_express_temp'
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

#   Improve label for plotting
metrics_df$corr <- paste("Cor =", metrics_df$corr)

#   For each cell type, plot counts of this cell type, measured by each
#   deconvolution tool, against the mean marker-gene expression for that cell
#   type (scatterplot)
plot_list <- lapply(
    cell_types_actual,
    function(ct) {
        full_df_small <- full_df %>%
            filter(cell_type == ct)
        
        p = ggplot(full_df_small) +
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
                    x = Inf, y = max(full_df_small$marker_express) / 7,
                    label = corr
                ),
                hjust = 1
            ) +
            labs(
                title = ct,
                x = paste0("Software-estimated ", ct, " counts"),
                y = "Mean marker-gene counts",
            ) +
            scale_fill_continuous(type = "viridis") +
            theme_bw(base_size = 10)
        
        return(p)
    }
)

pdf(file.path(plot_dir, 'marker_vs_ct_counts.pdf'))
print(plot_list)
dev.off()

#-------------------------------------------------------------------------------
#   Spatial distribution of cell-types compared against manual layer annotation
#   (boxplots)
#-------------------------------------------------------------------------------

#   Read layer annotation in for each sample and match to a barcode
layer_ann = data.frame()
for (sample_id in sample_ids) {
    this_layer_path <- sub("\\{sample_id\\}", sample_id, layer_ann_path)
    layer_ann_small <- read.csv(this_layer_path)
    
    layer_ann_small$barcode = colnames(spe[, spe$sample_id == sample_id])[
        layer_ann_small$id + 1
    ]
    layer_ann_small$sample_id = sample_id
    
    layer_ann = rbind(layer_ann, layer_ann_small)
}
layer_ann$id = NULL

#   Add layer label to observed_df_long
observed_df_long = left_join(
    observed_df_long, layer_ann, by = c('barcode', 'sample_id')
)

#   Clean up labels
observed_df_long$label = tolower(observed_df_long$label)

#   Average counts of each cell type in each layer as annotated; filter NA
#   labels (there are 2 inentional NAs where spots should be dropped)
counts_df = observed_df_long |> 
    filter(
        !is.na(label),
        label != ""
    ) |>
    group_by(label, deconvo_tool, sample_id, cell_type) |>
    summarize(count = mean(observed)) |>
    ungroup()

#   Create plots for each cell type
plot_list = list()
for (cell_type in cell_types) {
    y_max = counts_df |>
        filter(cell_type == {{ cell_type }}) |>
        summarize(y_max = max(count)) |>
        pull(y_max)
    
    plot_list[[cell_type]] = ggplot(
        counts_df |> filter(cell_type == {{ cell_type }}),
        aes(x = label, y = count, color = deconvo_tool)
    ) +
        geom_boxplot(outlier.shape = NA) +
        labs(
            x = "Annotated Layer",
            y = paste("Average Predicted", cell_type, "Count"),
            color = "Deconvolution Tool"
        ) +
        theme_bw(base_size = 20) +
        coord_cartesian(ylim = c(0, y_max)) +
        scale_y_continuous(expand = c(0, 0, 0, 0.05))
}

pdf(
    file.path(plot_dir, "layer_distribution.pdf"), width = 10
)
print(plot_list)
dev.off()

#-------------------------------------------------------------------------------
#   Spatial distribution of cell-types compared against manual layer annotation
#   (barplots- raw counts)
#-------------------------------------------------------------------------------

#   Read in SCE just to get colors for software-estimated cell types
sce = readRDS(sce_in)
estimated_cell_labels = metadata(sce)[[paste0('cell_type_colors_', cell_group)]]
names(estimated_cell_labels) = gsub('/', '_', names(estimated_cell_labels))

counts_df = counts_df |>
    group_by(label, deconvo_tool, cell_type) |>
    summarize(count = mean(count)) |>
    ungroup() |>
    #   "EndoMural" gets collapsed into "other", which we aren't considering
    filter(cell_type != "EndoMural")

#   Meaning of an example section of one bar in one facet of this plot:
#   Orange bar at white matter for cell2location: of all spots manually
#   annotated as white matter, the size of the bar section is the across-sample
#   mean of average oligo counts in those spots for each sample. Total bar
#   heights vary within deconvo tool because cell-count density varies between
#   layers. Total bar heights may vary between deconvo tools when total counts
#   per spot aren't preserved (e.g. C2L doesn't match the cell counts you
#   provide it)
pdf(
    file.path(plot_dir, 'layer_distribution_barplot_raw.pdf'), width = 10,
    height = 5
)
ggplot(
        counts_df,
        aes(x = label, y = count, fill = cell_type)
    ) +
    facet_wrap(~ deconvo_tool) +
    geom_bar(stat = "identity") +
    labs(
        x = "Annotated Layer", y = "Average Predicted Count",
        fill = "Cell Type"
    ) +
    scale_fill_manual(values = estimated_cell_labels) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90))
dev.off()

#-------------------------------------------------------------------------------
#   Spatial distribution of cell-types compared against manual layer annotation
#   (barplots- proportions by layer)
#-------------------------------------------------------------------------------

counts_df = observed_df_long |> 
    filter(!is.na(label), cell_type != "EndoMural") |>
    #   For each manually annotated label and deconvo tool, normalize by the
    #   total counts of all cell types and samples
    group_by(deconvo_tool, label) |>
    mutate(observed = observed / sum(observed)) |>
    #   Now for each label, deconvo tool and cell type, add up counts for all
    #   samples and relevant spots
    group_by(deconvo_tool, label, cell_type) |>
    summarize(observed = sum(observed)) |>
    #   Now add a column 'layer_match' to indicate rows where each cell type
    #   has a maximal proportion across layers. We'll mark these with an "X"
    #   on the barplots
    group_by(deconvo_tool, cell_type) |>
    mutate(layer_match = observed == max(observed)) |>
    ungroup()

if (cell_group == "layer") {
    #   Take just the excitatory cell types
    match_df = counts_df |>
        filter(layer_match, str_detect(cell_type, "^Excit"))
    
    #   Add a column 'is_match' indicating whether for an excitatory cell type
    #   and deconvo tool, the cell_type has maximal proportion in the correct/
    #   expected layer
    match_df$is_match = sapply(
        1:nrow(match_df),
        function(i) {
            match_df$label[i] %in%
                corresponding_layers[[as.character(match_df$cell_type)[i]]]
        }
    )
    
    #   For each deconvo tool, add up how many times excitatory cell types
    #   have maximal proportion in the correct layers
    print('Number of times excitatory cell types have maximal proportion in the correct layer:')
    match_df |>
        group_by(deconvo_tool) |>
        summarize(num_matches = sum(is_match))
}

pdf(
    file.path(plot_dir, 'layer_distribution_barplot_prop.pdf'), width = 10,
    height = 5
)
ggplot(
    counts_df,
    aes(x = label, y = observed, fill = cell_type)
) +
    facet_wrap(~ deconvo_tool) +
    geom_bar(stat = "identity") +
    labs(
        x = "Annotated Layer", y = "Proportion of Counts", fill = "Cell Type"
    ) +
    scale_fill_manual(values = estimated_cell_labels) +
    geom_text(
        aes(label = ifelse(layer_match, "X", "")),
        position = position_stack(vjust = 0.5)
    ) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90))
dev.off()

#-------------------------------------------------------------------------------
#   Spatial plots of manual layer annotation
#-------------------------------------------------------------------------------

plot_list = list()
for (sample_id in sample_ids) {
    spe_small = spe[, spe$sample_id == sample_id]
    counts_df = observed_df_long |> filter(sample_id == {{ sample_id }})
    
    spe_small$manual_layer = counts_df[
        match(colnames(spe_small), counts_df$barcode), 'label'
    ] |> pull(label)
    
    #   Verify no NA labels are present, except for sample 'Br6432_Ant_IF',
    #   where 2 NA spots are expected: these should be dropped for these plots
    spe_small$manual_layer[spe_small$manual_layer == ""] = NA
    if (sample_id != 'Br6432_Ant_IF' & any(is.na(spe_small$manual_layer))) {
        stop("No spot labels should be NA")
    } else if (sample_id == 'Br6432_Ant_IF') {
        stopifnot(length(which(is.na(spe_small$manual_layer))) == 2)
        spe_small = spe_small[, !is.na(spe_small$manual_layer)]
    }
    
    plot_list[[sample_id]] = vis_grid_gene(
        spe_small, geneid = 'manual_layer', return_plots = TRUE,
        spatial = FALSE
    )[[1]] +
        scale_color_discrete() +
        scale_fill_discrete() +
        theme_bw(base_size = 15) +
        coord_fixed()
}

pdf(file.path(plot_dir, 'spot_layer_labels.pdf'))
print(plot_list)
dev.off()
