library("here")
library("ggplot2")
library("jaffelab")
library("tidyverse")
library("reshape2")
library("SpatialExperiment")
library("sessioninfo")

#   Adds the 'layer_dist_barplot' function
source(
    here("code", "spot_deconvo", "05-shared_utilities", "shared_functions.R")
)

cell_group <- "layer" # "broad" or "layer"

#   Named vector where names are used in plots, and values give the directory
#   names where SPOTlight results are saved
subset_choices <- c(
    "Full Data" = "full_data",
    "100 Cells (Seeded)" = "subset_n100",
    "100 Cells (No Seed)" = "subset_n100_bad"
)

sample_ids_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "sample_ids.txt"
)

plot_dir <- here(
    "plots", "spot_deconvo", "04-spotlight", "IF", cell_group, "comparison"
)

processed_dir <- here(
    "processed-data", "spot_deconvo", "04-spotlight", "IF", cell_group,
    "comparison"
)

#   "clusters.csv" files of spot-deconvolution results
observed_paths <- here(
    "processed-data", "spot_deconvo", "04-spotlight", "IF",
    cell_group, "{subset_choice}", "{sample_id}", "clusters.csv"
)

spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)

#   Only needed to get colors for software-estimated cell types
sce_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("sce_", cell_group, ".rds")
)

layer_ann_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "annotations_{sample_id}_spots.csv"
)

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

set.seed(11282022)
dir.create(plot_dir, showWarnings = FALSE)
dir.create(processed_dir, showWarnings = FALSE)

################################################################################
#   Analysis
################################################################################

sample_ids <- readLines(sample_ids_path)

#-------------------------------------------------------------------------------
#   Read in SPE and SCE
#-------------------------------------------------------------------------------

#   Read in SCE just to get colors for software-estimated cell types
sce <- readRDS(sce_in)
estimated_cell_labels <- metadata(sce)[[paste0("cell_type_colors_", cell_group)]]
names(estimated_cell_labels) <- gsub("/", "_", names(estimated_cell_labels))

spe <- readRDS(spe_IF_in)

#-------------------------------------------------------------------------------
#   Read in deconvolution results for SPOTlight run with the different settings
#-------------------------------------------------------------------------------

#   Form a list of data frames containing the deconvolution results
observed_df_list = list()
i = 1

for (sample_id in sample_ids) {
    for (subset_choice in subset_choices) {
        #   Read in estimated cell counts
        observed_path <- sub("\\{sample_id\\}", sample_id, observed_paths)
        observed_path <- sub(
            "\\{subset_choice\\}", subset_choice, observed_path
        )
        observed_df_list[[i]] <- read.csv(observed_path) |>
            mutate(
                subset_choice = {{ subset_choice }},
                sample_id = {{ sample_id }},
                deconvo_tool = "SPOTlight",
                barcode = ss(key, "_", 1)
            )  |>
            select(- key)
        
        i = i + 1
    }
}

#-------------------------------------------------------------------------------
#   Read in layer annotation and merge with the deconvolution results
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

#   Collapse list of results into one tibble and add layer annotation
observed_df = do.call(rbind, observed_df_list) |>
    melt(
        id.vars = c("barcode", "sample_id", "deconvo_tool", "subset_choice"),
        variable.name = "cell_type", value.name = "count"
    ) |>
    #   Use an ordered factor for plotting cell types
    mutate(cell_type = factor(cell_type, levels = cell_types, order = TRUE)) |>
    left_join(layer_ann, by = c("barcode", "sample_id")) |>
    filter(!is.na(label)) |>
    as_tibble()

#   Clean up labels
observed_df$label <- tolower(observed_df$label)
observed_df$label <- sub("layer", "L", observed_df$label)
observed_df$label[observed_df$label == "wm"] <- "WM"
stopifnot(
    all(unlist(corresponding_layers) %in% unique(observed_df$label))
)

#-------------------------------------------------------------------------------
#   Format and plot data
#-------------------------------------------------------------------------------

#   Format data for plotting
count_df = observed_df |>
    #   For each manually annotated label, subset choice and sample_id,
    #   normalize by the total counts of all cell types
    group_by(subset_choice, label, sample_id) |>
    mutate(count = count / sum(count)) |>
    #   Now for each label, subset choice, sample_id and cell type, add up
    #   counts for all relevant spots
    group_by(subset_choice, label, cell_type, sample_id) |>
    summarize(count = sum(count)) |>
    #   Now average across samples
    group_by(subset_choice, label, cell_type) |>
    summarize(count = mean(count)) |>
    ungroup() |>
    #   Improve labels for plot facets; treat subset choice as the deconvolution
    #   tool as a trick to re-use plotting code ('layer_dist_barplot' function)
    mutate(
        deconvo_tool = names(subset_choices)[
            match(subset_choice, subset_choices)
        ]
    ) |>
    select(- subset_choice)

layer_dist_barplot(
    count_df, filename = "layer_distribution_barplot_prop_even.pdf",
    ylab = "Proportion of Counts", x_var = "label", fill_var = "cell_type",
    fill_scale = estimated_cell_labels, fill_lab = "Cell Type",
    xlab = "Annotated Layer"
)

session_info()
