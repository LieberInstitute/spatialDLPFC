library("here")
library("ggplot2")
library("jaffelab")
library("tidyverse")
library("reshape2")
library("spatialLIBD")
library("cowplot")

cell_group <- "layer" # "broad" or "layer"

sample_ids <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "sample_ids.txt"
)

deconvo_tools <- c("01-tangram", "03-cell2location", "04-spotlight")

#   "Ground-truth" cell counts from cellpose + trained classification tree
actual_paths <- here(
    "processed-data", "spot_deconvo", "02-cellpose", "{sample_id}", "clusters.csv"
)

#   Cell counts estimated by each deconvolution tool
observed_paths <- here(
    "processed-data", "spot_deconvo", "{deconvo_tool}", "IF", cell_group,
    "{sample_id}", "clusters.csv"
)

plot_dir <- here(
    "plots", "spot_deconvo", "05-shared_utilities", cell_group
)

spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)

marker_object_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("marker_stats_", cell_group, ".rds")
)

marker_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("markers_", cell_group, ".txt")
)

cell_types_actual <- c("micro", "neuron", "oligo", "other")
if (cell_group == "broad") {
    cell_types <- c("Astro", "Excit", "Inhib", "Micro", "Oligo", "OPC")
} else {
    cell_types <- c(
        "Astro", "Excit_L2_3", "Excit_L3", "Excit_L3_4_5", "Excit_L4",
        "Excit_L5", "Excit_L5_6", "Excit_L6", "Inhib", "Micro", "Oligo", "OPC"
    )
}

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
                geom_point(
                    aes(x = observed, y = actual, color = sample_id),
                    alpha = 0.01
                ) +
                geom_abline(
                    intercept = 0, slope = 1, linetype = "dashed", color = "red"
                ) +
                coord_fixed() +
                facet_grid(rows = vars(sample_id), cols = vars(deconvo_tool)) +
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
                scale_fill_continuous(type = "viridis") +
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
across_spots <- function(count_df, plot_name) {
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
            aes(x = observed, y = actual, color = cell_type, shape = sample_id)
        ) +
        coord_fixed() +
        facet_wrap(~deconvo_tool) +
        geom_abline(
            intercept = 0, slope = 1, linetype = "dashed", color = "red"
        ) +
        geom_text(
            data = metrics_df,
            mapping = aes(
                x = Inf, y = max(count_df$observed) / 7, label = corr
            ),
            hjust = 1
        ) +
        geom_text(
            data = metrics_df,
            mapping = aes(x = Inf, y = 0, label = rmse),
            hjust = 1, vjust = 0
        ) +
        labs(x = "Software-estimated", y = "CART-calculated") +
        theme_bw(base_size = 10)

    pdf(file.path(plot_dir, plot_name))
    print(p)
    dev.off()
    # return(p)
}

################################################################################
#   Read in and format cell counts into a table apt for plotting with ggplot
################################################################################

sample_ids <- readLines(sample_ids)

actual_list <- list()
observed_list <- list()
index <- 1

added_colnames <- c("barcode", "sample_id", "deconvo_tool", "obs_type")

#   Read in counts for all cell types for all samples and deconvolution tools
#   (as well as the "ground-truth")
for (sample_id in sample_ids) {
    #   Read in the "ground-truth" counts for this sample
    actual_path <- sub("\\{sample_id\\}", sample_id, actual_paths)
    actual_df_small <- read.csv(actual_path)

    #   Make sure we have just counts for each cell type, barcode, and
    #   sample ID variables-- nothing else. Order column names
    actual_df_small$barcode <- ss(actual_df_small$key, "_", 1)
    actual_df_small$sample_id <- sample_id
    actual_df_small$obs_type <- "actual"

    for (deconvo_tool in deconvo_tools) {
        #   Read in estimated cell counts for this deconvo tool and sample
        observed_path <- sub("\\{sample_id\\}", sample_id, observed_paths)
        observed_path <- sub("\\{deconvo_tool\\}", deconvo_tool, observed_path)
        observed_df_small <- read.csv(observed_path)
        colnames(observed_df_small) = gsub(
            '\\.', '_', colnames(observed_df_small)
        )

        #   Make sure column names are consistent and include only info about
        #   barcode, sample_id, deconvo tool, and cell-type counts
        observed_df_small$barcode <- ss(observed_df_small$key, "_", 1)
        observed_df_small$sample_id <- sample_id
        observed_df_small$deconvo_tool <- deconvo_tool
        observed_df_small$obs_type <- "observed"
        observed_df_small <- observed_df_small[
            , c(added_colnames, cell_types)
        ]

        actual_df_small$deconvo_tool <- deconvo_tool

        actual_list[[index]] <- actual_df_small[
            , c(added_colnames, cell_types_actual)
        ]
        observed_list[[index]] <- observed_df_small
        index <- index + 1
    }
}

#   Form data frames containing cell-type counts for all spots, samples, and
#   deconvolution tools
actual_df <- as_tibble(do.call(rbind, actual_list))
observed_df <- as_tibble(do.call(rbind, observed_list))

#   Combine cell types as appropriate for comparison against the relatively
#   narrow types in the ground-truth
colnames(observed_df) <- tolower(colnames(observed_df))

if (cell_group == "broad") {
    observed_df <- observed_df %>%
        mutate("neuron" = excit + inhib) %>%
        mutate("other" = astro + opc) %>%
        select(all_of(c(added_colnames, cell_types_actual)))
} else {
    observed_df = observed_df %>%
        rowwise() %>%
        mutate(
            neuron = sum(c_across(starts_with(c('excit_', 'inhib')))),
            other = astro + opc
        ) %>%
        ungroup() %>%
        select(all_of(c(added_colnames, cell_types_actual)))
}

stopifnot(all(colnames(observed_df) == colnames(actual_df)))

#   Combine observed and actual cell counts so that each row is a unique
#   cell type, spot, sample, and deconvolution method with two values: measured
#   and ground-truth
full_df <- rbind(observed_df, actual_df) %>%
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

#-------------------------------------------------------------------------------
#   Plot total counts per spot for tangram and cell2location
#-------------------------------------------------------------------------------

#   For each spot, plot the provided vs. computed total number of cells for the
#   deconvolution methods that aren't constrained to use the same totals as
#   provided (tangram and cell2location)
count_df <- full_df %>%
    filter(deconvo_tool %in% c("01-tangram", "03-cell2location")) %>%
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

pdf(file.path(plot_dir, 'total_cells.pdf'))
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
            x = Inf, y = max(count_df$observed) / 7, label = corr
        ),
        hjust = 1
    ) +
    geom_text(
        data = metrics_df,
        mapping = aes(x = Inf, y = 0, label = rmse),
        hjust = 1, vjust = 0
    ) +
    labs(
        x = "Calculated cell count",
        y = "Provided cell count (cellpose)",
        title = "Provided vs. calculated total cells per spot"
    ) +
    theme_bw(base_size = 10)
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

across_spots(count_df, "counts_across_spots_scatter.pdf")

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

across_spots(prop_df, "props_across_spots_scatter.pdf")

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

across_spots(count_df, "adjusted_counts_across_spots_scatter.pdf")

#-------------------------------------------------------------------------------
#   Spatial distribution of counts for each sample, cell_type, deconvo method
#-------------------------------------------------------------------------------

#   TODO: This is hardcoded for IF! Not sure if we'll do similar plots for nonIF

spe = readRDS(spe_IF_in)

for (sample_id in sample_ids) {
    spe_small = spe[,spe$sample_id == sample_id]
    
    i = 1
    plot_list = list()
    for (deconvo_tool in deconvo_tools) {
        for (cell_type in cell_types_actual) {
            #   Grab counts for just this sample, deconvo tool, and cell type
            counts_df = full_df %>%
                filter(
                    sample_id == {{ sample_id }},
                    deconvo_tool == {{ deconvo_tool }},
                    cell_type == {{ cell_type }}
                )
            
            #   Add counts as a column in colData(spe_small)
            spe_small$temp_ct_counts = counts_df$observed[
                match(colnames(spe_small), counts_df$barcode)
            ]
            
            #   Plot spatial distribution
            plot_list[[i]] = vis_grid_gene(
                spe_small, geneid = 'temp_ct_counts', return_plots = TRUE,
                spatial = FALSE
            )[[1]] +
                labs(
                    title = paste0(
                        cell_type, ' counts\n(', deconvo_tool, ')'
                    )
                )
            
            i = i + 1
        }
    }
    
    #   Plot in a grid where cell types are columns and rows are deconvolution
    #   tools. One PDF per sample
    pdf(
        file.path(plot_dir, paste0('spatial_counts_', sample_id, '.pdf')),
        width = 7 * length(cell_types_actual),
        height = 7 * length(deconvo_tools)
    )
    print(plot_grid(plotlist = plot_list, ncol = length(cell_types_actual)))
    dev.off()
}

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
marker_stats$cellType.target[
    marker_stats$cellType.target %in% c('opc', 'astro')
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
        corr = round(cor(observed, actual), 2),
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
            facet_grid(rows = vars(sample_id), cols = vars(deconvo_tool)) +
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
