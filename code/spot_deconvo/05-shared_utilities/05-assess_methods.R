library('here')
library('ggplot2')
library('jaffelab')
library('tidyverse')
library('reshape2')

cell_group = "broad" # "broad" or "layer"

sample_ids = here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "sample_ids.txt"
)

deconvo_tools = c('01-tangram', '03-cell2location', '04-spotlight')

#   "Ground-truth" cell counts from cellpose + trained classification tree
actual_paths = here(
    "processed-data", "spot_deconvo", "02-cellpose", "{sample_id}", "clusters.csv"
)

#   Cell counts estimated by each deconvolution tool
observed_paths <- here(
    "processed-data", "spot_deconvo", "{deconvo_tool}", "IF", cell_group,
    "{sample_id}", "clusters.csv"
)

plot_dir = here(
    "plots", "spot_deconvo", "05-shared_utilities", cell_group
)

cell_types_actual = c('micro', 'neuron', 'oligo', 'other')
if (cell_group == "broad") {
    cell_types = c('Astro', 'Excit', 'Inhib', 'Micro', 'Oligo', 'OPC')
} else {
    cell_types = c() # TODO
}

################################################################################
#   Plotting functions
################################################################################

#   Scatterplots observed vs. actual cell counts for each cell type, faceted
#   by sample and deconvolution tool. Use all spots as points
all_spots = function(count_df, plot_name) {
    #   Compute metrics for each deconvolution tool: correlation between
    #   observed and actual values as well as RMSE
    metrics_df = count_df %>%
        group_by(deconvo_tool, sample_id) %>%
        summarize(
            corr = round(cor(observed, actual), 2),
            rmse = signif(mean((observed - actual) ** 2) ** 0.5, 3)
        )
    
    #   Improve labels for plotting
    metrics_df$corr = paste('Cor =', metrics_df$corr)
    metrics_df$rmse = paste('RMSE =', metrics_df$rmse)
    
    plot_list = lapply(
        cell_types_actual,
        function(ct) {
            count_df_small = count_df %>%
                filter(cell_type == ct)
            
            ggplot(count_df_small) +
                geom_point(
                    aes(x = observed, y = actual, color = sample_id),
                    alpha = 0.01
                ) +
                geom_abline(
                    intercept = 0, slope = 1, linetype = 'dashed', color = 'red'
                ) +
                coord_fixed() +
                facet_grid(rows = vars(sample_id), cols = vars(deconvo_tool)) +
                guides(col = guide_legend(override.aes = list(alpha = 1))) +
                labs(title = ct) +
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
                )
        }
    )
    
    # pdf(file.path(plot_dir, plot_name))
    # print(plot_list)
    # dev.off()
    return(plot_list)
}

#   Scatterplot of observed vs. actual total counts summed across spots,
#   faceted by cell type and deconvolution tool
across_spots = function(count_df, plot_name) {
    #   Compute metrics for each deconvolution tool: correlation between
    #   observed and actual values as well as RMSE
    metrics_df = count_df %>%
        group_by(deconvo_tool) %>%
        summarize(
            corr = round(cor(observed, actual), 2),
            rmse = signif(mean((observed - actual) ** 2) ** 0.5, 3)
        )
    
    #   Improve labels for plotting
    metrics_df$corr = paste('Cor =', metrics_df$corr)
    metrics_df$rmse = paste('RMSE =', metrics_df$rmse)
    
    p = ggplot(count_df) +
        geom_point(
            aes(x = observed, y = actual, color = cell_type, shape = sample_id)
        ) +
        coord_fixed() +
        facet_wrap(~deconvo_tool) +
        geom_abline(
            intercept = 0, slope = 1, linetype = 'dashed', color = 'red'
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
        )
    
    pdf(file.path(plot_dir, plot_name))
    print(p)
    dev.off()
    # return(p)
}

################################################################################
#   Analysis/ data exploration
################################################################################

#-------------------------------------------------------------------------------
#   Read in and format cell counts into a table apt for plotting with ggplot
#-------------------------------------------------------------------------------

sample_ids = readLines(sample_ids)

actual_list = list()
observed_list = list()
index = 1

added_colnames = c("barcode", "sample_id", "deconvo_tool", "obs_type")

#   Read in counts for all cell types for all samples and deconvolution tools
#   (as well as the "ground-truth")
for (sample_id in sample_ids) {
    #   Read in the "ground-truth" counts for this sample
    actual_path = sub("\\{sample_id\\}", sample_id, actual_paths)
    actual_df_small = read.csv(actual_path)
    
    #   Make sure we have just counts for each cell type, barcode, and
    #   sample ID variables-- nothing else. Order column names
    actual_df_small$barcode = ss(actual_df_small$key, '_', 1)
    actual_df_small$sample_id = sample_id
    actual_df_small$obs_type = "actual"
    
    for (deconvo_tool in deconvo_tools) {
        #   Read in estimated cell counts for this deconvo tool and sample
        observed_path = sub("\\{sample_id\\}", sample_id, observed_paths)
        observed_path = sub("\\{deconvo_tool\\}", deconvo_tool, observed_path)
        observed_df_small = read.csv(observed_path)
        
        #   Make sure column names are consistent and include only info about
        #   barcode, sample_id, deconvo tool, and cell-type counts
        observed_df_small$barcode = ss(observed_df_small$key, '_', 1)
        observed_df_small$sample_id = sample_id
        observed_df_small$deconvo_tool = deconvo_tool
        observed_df_small$obs_type = "observed"
        observed_df_small = observed_df_small[
            , c(added_colnames, cell_types)
        ]
        
        actual_df_small$deconvo_tool = deconvo_tool
        
        actual_list[[index]] = actual_df_small[
            , c(added_colnames, cell_types_actual)
        ]
        observed_list[[index]] = observed_df_small
        index = index + 1
    }
}

#   Form data frames containing cell-type counts for all spots, samples, and
#   deconvolution tools
actual_df = as_tibble(do.call(rbind, actual_list))
observed_df = as_tibble(do.call(rbind, observed_list))

#   Combine cell types as appropriate for comparison against the relatively
#   narrow types in the ground-truth
if (cell_group == "broad") {
    colnames(observed_df) = tolower(colnames(observed_df))
    
    observed_df = observed_df %>%
        mutate("neuron" = excit + inhib) %>%
        mutate("other" = astro + opc) %>%
        select(all_of(c(added_colnames, cell_types_actual)))
    
    stopifnot(all(colnames(observed_df) == colnames(actual_df)))
} else {
    #   TODO
}

#   Combine observed and actual cell counts so that each row is a unique
#   cell type, spot, sample, and deconvolution method with two values: measured
#   and ground-truth
full_df = rbind(observed_df, actual_df) %>%
    melt(
        id.vars = added_colnames, variable.name = "cell_type",
        value.name = "count"
    ) %>%
    pivot_wider(
        names_from = obs_type, values_from = count,
    )

#-------------------------------------------------------------------------------
#   Exploratory plots
#-------------------------------------------------------------------------------

#   For each spot, plot the provided vs. computed total number of cells for the
#   deconvolution methods that aren't constrained to use the same totals as
#   provided (tangram and cell2location)
count_df = full_df %>%
    filter(deconvo_tool %in% c('01-tangram', '03-cell2location')) %>%
    group_by(barcode, sample_id, deconvo_tool) %>%
    summarize(observed = sum(observed), actual = sum(actual))

ggplot(count_df) +
    geom_point(aes(x = observed, y = actual), alpha = 0.01) +
    facet_wrap(~ deconvo_tool) +
    coord_fixed() +
    geom_abline(
        intercept = 0, slope = 1, linetype = 'dashed', color = 'red'
    ) +
    labs(title = "Provided vs. calculated total cells per spot")

all_spots(full_df, 'counts_all_spots_scatter.pdf')

count_df = full_df %>%
    group_by(sample_id, deconvo_tool, cell_type) %>%
    summarize(observed = sum(observed), actual = sum(actual))



across_spots(count_df, 'counts_across_spots_scatter.pdf')

#   Now convert from counts to proportions for each spot
prop_df = full_df %>%
    group_by(barcode, sample_id, deconvo_tool) %>%
    mutate(
        observed = observed / sum(observed),
        actual = actual / sum(actual),
    )
all_spots(prop_df, 'props_all_spots_scatter.pdf')

#   First, sum counts for each (cell type, sample, deconvo tool) across
#   spots. Then assign a cell-type proportion by dividing by total cells
#   for each (sample, deconvo tool). Note this method yields many more non-NA
#   values since taking proportions in each spot and averaging across spots
#   introduces NAs whenever either the observed or actual total cell count is 0
#   for a spot.
prop_df = full_df %>%
    group_by(sample_id, deconvo_tool, cell_type) %>%
    summarize(observed = sum(observed), actual = sum(actual)) %>%
    group_by(sample_id, deconvo_tool) %>%
    mutate(
        observed = observed / sum(observed),
        actual = actual / sum(actual),
    )
across_spots(prop_df, 'props_across_spots_scatter.pdf')
