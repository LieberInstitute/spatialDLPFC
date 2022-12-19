library("here")
library("ggplot2")
library("jaffelab")
library("tidyverse")
library("reshape2")
library("spatialLIBD")
library("cowplot")
library("sessioninfo")

#   Adds the 'spot_plot' function, a wrapper for 'vis_gene' or 'vis_clus' with
#   consistent manuscript-appropriate settings
source(
    here("code", "spot_deconvo", "05-shared_utilities", "shared_functions.R")
)

cell_group <- "layer" # "broad" or "layer"

sample_ids_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF",
    "sample_ids.txt"
)

raw_results_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF",
    paste0("results_raw_", cell_group, ".csv")
)

collapsed_results_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF",
    paste0("results_collapsed_", cell_group, ".csv")
)

plot_dir <- here(
    "plots", "spot_deconvo", "05-shared_utilities", "nonIF", cell_group
)

spe_nonIF_in <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
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

deconvo_tools <- c("Tangram", "Cell2location", "SPOTlight")

cell_types_actual <- c("Astro", "Micro", "Neuron", "Oligo", "Other")
if (cell_group == "broad") {
    cell_types <- c(
        "Astro", "EndoMural", "Micro", "Oligo", "OPC", "Excit", "Inhib"
    )
} else {
    cell_types <- c(
        "Astro", "EndoMural", "Micro", "Oligo", "OPC", "Excit_L2_3", "Excit_L3",
        "Excit_L3_4_5", "Excit_L4", "Excit_L5", "Excit_L5_6", "Excit_L6",
        "Inhib"
    )
}

cell_type_labels <- c("#3BB273", "#663894", "#E49AB0", "#E07000", "#95B8D1")
names(cell_type_labels) <- cell_types_actual

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
#   Functions
################################################################################

#   Return a length-2 list containing a 'vis_grid_gene' plot and maximum value
#   in the plot given the SpatialExperiment, long-format
#   tibble of cell-type counts, the target sample ID, deconvo tool, and cell
#   type, column name of 'full_df' ('observed' or 'actual'), and a plot title
spatial_counts_plot <- function(spe_small, full_df, sample_id1, deconvo_tool1, cell_type1, c_name,
                                title) {
    #   Grab counts for just this sample, deconvo tool, and cell type
    counts_df <- full_df |>
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
                    paste0(cell_type, " counts (", deconvo_tool, ")")
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
                    spe_small, full_df, sample_id, deconvo_tool, cell_type, "actual",
                    paste0(cell_type, " counts (CART-calculated)")
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
            ),
            width = 7 * length(cell_type_vec),
            height = 7 * (length(deconvo_tools) + include_actual)
        )
        print(plot_grid(plotlist = plot_list, ncol = length(cell_type_vec)))
        dev.off()
    }
}

#   Write a PDF to 'filename' of a barplot of section-wide cell-type
#   proportions
#
#   prop_df: tibble with columns 'sample_id', 'deconvo_tool', and 'prop'. One
#       row per sample per deconvo tool per cell type.
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
            scale_fill_manual(values = estimated_cell_labels) +
            theme_bw(base_size = 16)
    }
    pdf(file.path(plot_dir, filename))
    print(plot_list)
    dev.off()
}

#   Write a PDF to 'filename' of a barplot of section-wide cell-type
#   proportions. Intended as a supplementary figure in the manuscript
#
#   prop_df: tibble with columns 'sample_id', 'deconvo_tool', and 'prop'.
#       One row per sample per deconvo tool per cell type.
#   filename: character relative to 'plot_dir', with extension ".pdf"
prop_barplot_paper <- function(prop_df, filename) {
    p = ggplot(prop_df, aes(x = deconvo_tool, y = prop, fill = cell_type)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ sample_id, nrow = 5, ncol = 6) +
        labs(x = "Method", y = "Sample-Wide Proportion", fill = "Cell Type") +
        scale_x_discrete(labels = c("actual" = "CART")) +
        scale_fill_manual(values = estimated_cell_labels) +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    pdf(file.path(plot_dir, filename))
    print(p)
    dev.off()
}

#   Given a tibble with columns 'label' (manual layer label), 'deconvo_tool',
#   'cell_type', and 'count', write a set of barplots to PDF under [plot_dir]
#   with name [filename]. 'ylab' give the y-axis label; 'x_var' is the x-axis
#   variable (as a string); 'fill_var' is the fill variable as a string;
#   'fill_scale' is passed to 'scale_fill_manual(values = [fill_scale])';
#   'fill_lab' is the fill label; 'xlab' is the x-axis label
layer_dist_barplot <- function(
        counts_df, filename, ylab, x_var, fill_var, fill_scale, fill_lab, xlab
        ) {
    p <- ggplot(
        counts_df,
        aes_string(x = x_var, y = "count", fill = fill_var)
    ) +
        facet_wrap(~deconvo_tool) +
        geom_bar(stat = "identity") +
        labs(x = xlab, y = ylab, fill = fill_lab) +
        scale_fill_manual(values = fill_scale) +
        theme_bw(base_size = 16) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    pdf(file.path(plot_dir, filename), width = 10, height = 5)
    print(p)
    dev.off()
}

################################################################################
#   Read in and format cell counts into a table apt for plotting with ggplot
################################################################################

#   Read in SCE just to get colors for software-estimated cell types
sce <- readRDS(sce_in)
estimated_cell_labels <- metadata(sce)[[paste0("cell_type_colors_", cell_group)]]
names(estimated_cell_labels) <- gsub("/", "_", names(estimated_cell_labels))

sample_ids <- readLines(sample_ids_path)
added_colnames <- c("barcode", "sample_id", "deconvo_tool", "obs_type")

observed_df <- as_tibble(read.csv(raw_results_path))
observed_df$obs_type <- "observed"

#   Use titlecase (note this is MUCH faster than 'toTitleCase')
observed_df$deconvo_tool[observed_df$deconvo_tool == "tangram"] <- "Tangram"
observed_df$deconvo_tool[observed_df$deconvo_tool == "cell2location"] <-
    "Cell2location"

#   Plot counts for each cell type without collapsing cell categories
observed_df_long <- observed_df |>
    melt(
        id.vars = added_colnames, variable.name = "cell_type",
        value.name = "count"
    ) |>
    pivot_wider(
        names_from = obs_type, values_from = count,
    ) |>
    #   Use an ordered factor for plotting cell types
    mutate(cell_type = factor(cell_type, levels = cell_types, order = TRUE))

load(spe_nonIF_in, verbose = TRUE)

################################################################################
#   Exploratory plots
################################################################################

spatial_counts_plot_full(
    spe, observed_df_long, cell_types, FALSE, "spatial_counts_fullres_"
)

#-------------------------------------------------------------------------------
#   Stacked barplot of cell-type proportions
#-------------------------------------------------------------------------------

prop_df <- observed_df_long |>
    group_by(sample_id, deconvo_tool, cell_type) |>
    summarize(observed = sum(observed)) |>
    group_by(sample_id, deconvo_tool) |>
    mutate(observed = observed / sum(observed)) |>
    ungroup() |>
    pivot_longer(
        cols = c("observed"), values_to = "prop", names_to = "source"
    )

prop_barplot(prop_df, "prop_barplots.pdf")
prop_barplot(
    prop_df |> filter(cell_type != "EndoMural"), "prop_barplots_no_endo.pdf"
)
prop_barplot_paper(
    prop_df |> filter(cell_type != "EndoMural"),
    "prop_barplots_no_endo_paper.pdf"
)

#-------------------------------------------------------------------------------
#   Plot total counts per sample for tangram and cell2location
#-------------------------------------------------------------------------------

#   Grab provided cell counts per spot and sum across section. Then gather into
#   a tibble easily mergable with 'observed_df_long'
provided_total <- c()
for (sample_id in sample_ids) {
    provided_total <- c(
        provided_total, rep(sum(spe$count[spe$sample_id == sample_id]), 2)
    )
}
provided_df <- as_tibble(
    data.frame(
        sample_id = rep(sample_ids, each = 2), actual = provided_total,
        deconvo_tool = rep(c("Tangram", "Cell2location"), length(sample_ids))
    )
)

#   Merge with the estimated section-wide cell counts for tangram and C2L
count_df <- observed_df_long |>
    filter(deconvo_tool %in% c("Tangram", "Cell2location")) |>
    group_by(sample_id, deconvo_tool) |>
    summarize(observed = sum(observed)) |>
    ungroup() |>
    left_join(provided_df, by = c("sample_id", "deconvo_tool")) |>
    mutate(diff = abs((observed - actual) / (observed + actual)))

pdf(file.path(plot_dir, "total_cells_sample.pdf"))
ggplot(count_df, aes(x = deconvo_tool, y = diff, color = sample_id)) +
    geom_point() +
    theme_bw(base_size = 10) +
    scale_y_continuous(
        limits = c(0, max(count_df$diff) * 1.1), expand = c(0, 0)
    ) +
    labs(x = "Deconvolution Tool", y = "Absolute Proportion Difference")
dev.off()

#-------------------------------------------------------------------------------
#   Plot total counts per spot for tangram and cell2location
#-------------------------------------------------------------------------------

provided_df <- tibble(
    barcode = rep(colnames(spe), each = 2),
    actual = rep(spe$count, each = 2),
    sample_id = factor(rep(spe$sample_id, each = 2)),
    deconvo_tool = factor(rep(c("Tangram", "Cell2location"), ncol(spe)))
)

count_df <- observed_df_long |>
    filter(deconvo_tool %in% c("Tangram", "Cell2location")) |>
    group_by(barcode, sample_id, deconvo_tool) |>
    summarize(observed = sum(observed)) |>
    ungroup() |>
    left_join(provided_df, by = c("barcode", "sample_id", "deconvo_tool"))

#   Compute metrics for each deconvolution tool: correlation between
#   observed and actual values as well as RMSE
metrics_df <- count_df |>
    group_by(deconvo_tool) |>
    summarize(
        corr = round(cor(observed, actual), 2),
        rmse = signif(mean((observed - actual)**2)**0.5, 3)
    ) |>
    ungroup()

#   Improve labels for plotting
metrics_df$corr <- paste("Cor =", metrics_df$corr)
metrics_df$rmse <- paste("RMSE =", metrics_df$rmse)

pdf(file.path(plot_dir, "total_cells_spot.pdf"), height = 5, width = 10)
ggplot(count_df) +
    # geom_point(aes(x = observed, y = actual), alpha = 0.01) +
    geom_bin2d(aes(x = log(observed + 1), y = log(actual + 1))) +
    scale_fill_continuous(type = "viridis") +
    facet_wrap(~deconvo_tool) +
    geom_abline(
        intercept = 0, slope = 1, linetype = "dashed", color = "red"
    ) +
    geom_text(
        data = metrics_df,
        mapping = aes(
            x = max(log(count_df$observed + 1)),
            y = max(log(count_df$actual + 1)) / 8,
            label = corr
        ),
        hjust = 1, size = 6
    ) +
    geom_text(
        data = metrics_df,
        mapping = aes(x = max(log(count_df$observed + 1)), y = 0, label = rmse),
        hjust = 1, vjust = 0, size = 6
    ) +
    labs(
        x = "log(Calculated + 1)",
        y = "log(Provided:VistoSeg + 1)",
        title = "Provided vs. Calculated Total Cells Per Spot"
    ) +
    theme_bw(base_size = 16)
dev.off()

#-------------------------------------------------------------------------------
#   Spatial distribution of marker gene expression vs. cell-type counts
#-------------------------------------------------------------------------------

#   Read in markers and marker stats
marker_stats <- readRDS(marker_object_in)
markers <- readLines(marker_in)

#   Subset marker_stats to include just the rows where each gene is a marker for
#   the given cell type. Replace '/' with '_' in cell types
marker_stats <- marker_stats |>
    group_by(gene) |>
    filter(ratio == max(ratio), gene %in% markers) |>
    ungroup() |>
    mutate(cellType.target = gsub('/', '_', cellType.target))

stopifnot(all(marker_stats$cellType.target %in% cell_types))

#   Add mean marker-gene expression for each associated cell type to the main
#   data frame with cell-type counts
observed_df_long$marker_express <- NA
for (cell_type in cell_types) {
    #   Grab markers for this cell type
    these_markers <- marker_stats |>
        filter(cellType.target == cell_type) |>
        arrange(desc(ratio)) |>
        head(n = 25) |>
        pull(gene)

    stopifnot(all(these_markers %in% rownames(spe)))

    for (sample_id in sample_ids) {
        spe_small <- spe[these_markers, spe$sample_id == sample_id]

        for (deconvo_tool in deconvo_tools) {
            #   Form a temporary tibble with mean expression across markers for
            #   each spot
            temp_df <- tibble(
                "marker_express_temp" = colMeans(assays(spe_small)$counts),
                "barcode" = colnames(spe_small),
                "sample_id" = sample_id,
                "cell_type" = cell_type,
                "deconvo_tool" = deconvo_tool
            )

            #   Add this quantity to the full tibble that has cell-type counts
            temp_df <- left_join(
                observed_df_long, temp_df,
                by = c("barcode", "sample_id", "cell_type", "deconvo_tool")
            )
            new_indices <- !is.na(temp_df$marker_express_temp)
            observed_df_long[new_indices, "marker_express"] <- temp_df[
                new_indices, "marker_express_temp"
            ]
        }
    }
}

#   All rows in observed_df_long should've been filled with an expression value
stopifnot(all(!is.na(observed_df_long$marker_express)))

#   Compute correlation for each deconvolution tool
metrics_df <- observed_df_long |>
    group_by(deconvo_tool, sample_id, cell_type) |>
    summarize(
        corr = round(cor(observed, marker_express), 2),
    ) |>
    ungroup()

#   Improve label for plotting
metrics_df$corr <- paste("Cor =", metrics_df$corr)

#   Now plot a version suitable as a supplementary figure: this time, one page
#   per sample and facet by cell type
plot_list <- lapply(
    sample_ids,
    function(this_sample_id) {
        full_df_small <- observed_df_long |>
            filter(sample_id == this_sample_id)
        
        p <- ggplot(full_df_small) +
            geom_point(
                aes(x = observed, y = marker_express, color = cell_type),
                alpha = 0.01
            ) +
            facet_grid( rows = vars(cell_type), cols = vars(deconvo_tool)) +
            guides(col = guide_legend(override.aes = list(alpha = 1))) +
            geom_text(
                data = metrics_df |> filter(sample_id == this_sample_id),
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

#   Add k=9 spatial domains into the big tibble (observed_df_long)
temp_df = tibble(
    "barcode" = colnames(spe), "bs_k9" = spe$bayesSpace_pca_9,
    "sample_id" = spe$sample_id
) |>
    mutate(bs_k9 = paste0("Sp09D0", bs_k9))

observed_df_long <- left_join(
    observed_df_long, temp_df, by = c("barcode", "sample_id")
)
stopifnot(!any(is.na(observed_df_long$bs_k9)))

counts_df <- observed_df_long |>
    filter(deconvo_tool %in% c("Tangram", "Cell2location")) |>
    #   Average counts within sample for a given domain/ cell type/ deconvo tool
    group_by(bs_k9, deconvo_tool, sample_id, cell_type) |>
    summarize(count = mean(observed)) |>
    #   Average these averages across sample
    group_by(bs_k9, deconvo_tool, cell_type) |>
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
    xlab = "Annotated Layer", ylab = "Average Predicted Count", x_var = "bs_k9",
    fill_lab = "Cell Type", fill_var = "cell_type",
    fill_scale = estimated_cell_labels
)

#   Weight each sample equally
counts_df <- observed_df_long |>
    filter(deconvo_tool %in% c("Tangram", "Cell2location")) |>
    #   For each manually annotated domain, deconvo tool and sample_id,
    #   normalize by the total counts of all cell types
    group_by(deconvo_tool, bs_k9, sample_id) |>
    mutate(observed = observed / sum(observed)) |>
    #   Now for each domain, deconvo tool, sample_id and cell type, add up
    #   counts for all relevant spots
    group_by(deconvo_tool, bs_k9, cell_type, sample_id) |>
    summarize(count = sum(observed)) |>
    #   Now average across samples. Note that some samples have 0 counts for
    #   SPOTlight and Tangram in some domains, so we use 'na.rm = TRUE' here to
    #   just average over the samples that have nonzero counts
    group_by(deconvo_tool, bs_k9, cell_type) |>
    summarize(count = mean(count, na.rm = TRUE)) |>
    ungroup()

layer_dist_barplot(
    counts_df,
    filename = "layer_distribution_barplot_prop_even.pdf",
    xlab = "Annotated Layer", ylab = "Proportion of Counts", x_var = "bs_k9",
    fill_lab = "Cell Type", fill_var = "cell_type",
    fill_scale = estimated_cell_labels
)

#-------------------------------------------------------------------------------
#   Spatial distribution of cell-types compared against BayesSpace k = 9 domains
#   (barplots- proportions by domain). Inverted so that x-axis is cell type
#-------------------------------------------------------------------------------

counts_df <- observed_df_long |>
    filter(deconvo_tool %in% c("Tangram", "Cell2location")) |>
    #   Normalize counts by the sum within sample ID, cell type, and deconvo
    #   tool (each sample contributes equally)
    group_by(sample_id, deconvo_tool, cell_type) |>
    mutate(observed = observed / sum(observed)) |>
    #   Add up counts for all spots in each layer for each sample ID, cell
    #   type, and deconvo tool
    group_by(deconvo_tool, bs_k9, sample_id, cell_type) |>
    summarize(count = sum(observed)) |>
    #   Now average counts across samples (Note that
    #   'sum(count) / length(sample_ids)' must be used and not 'mean(count)',
    #   since not all samples have counts for all cell types in all domains
    group_by(deconvo_tool, bs_k9, cell_type) |>
    summarize(count = sum(count) / length(sample_ids)) |>
    ungroup()

#   Interpretation of this plot:
#   For a given deconvo tool, examine one cell type. The size of each component
#   of the bar can be interpreted as the probability a randomly selected cell of
#   this cell type belongs in this particular domain. Note that this does not
#   normalize for domain size (i.e. if more spots belong to domain 3 than white
#   matter, the maximal domain for a given cell type is more likely to be domain
#   3)!
layer_dist_barplot(
    counts_df,
    filename = "layer_distribution_barplot_probability_inverted.pdf",
    xlab = "Cell Type", ylab = "Proportion of Counts", x_var = "cell_type",
    fill_lab = "Annotated Layer", fill_var = "bs_k9",
    fill_scale = libd_layer_colors
)

session_info()
