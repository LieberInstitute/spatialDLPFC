library("here")
library("SpatialExperiment")
library("ggplot2")
library("sessioninfo")
library("tidyverse")

spe_path <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
)

sce_path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"

plot_dir <- here("plots", "spot_deconvo", "05-shared_utilities")

spe <- readRDS(spe_path)
load(sce_path, verbose = TRUE)

#   Form a tibble with sample ID, cell type, and the proportion of cells of that
#   type, for the snRNA-seq data
prop_df <- colData(sce) |>
    as_tibble() |>
    rename(sample_id = Sample) |>
    #   For each sample, sum up counts of each cell type
    group_by(sample_id, cellType_broad_hc) |>
    summarize(n = n()) |>
    #   Now convert to cell-type proportions
    group_by(sample_id) |>
    mutate(
        cellType_broad_hc_prop = n / sum(n),
        position = colData(sce)[match(sample_id, sce$Sample), "pos"]
    ) |>
    select(-n) |>
    ungroup()

#   Now add the proportion of white matter for each spatial sample to 'prop_df'
prop_df <- colData(spe) |>
    as_tibble() |>
    #   Grab only the spatial samples with matching snRNA-seq data
    filter(sample_id %in% unique(prop_df$sample_id)) |>
    group_by(sample_id) |>
    summarize(
        prop_WM = sum(BayesSpace_harmony_28 %in% c(6, 16, 17, 20, 28)) / n()
    ) |>
    ungroup() |>
    right_join(prop_df, multiple = "all")

pdf(file.path(plot_dir, "prop_WM_position_scatter.pdf"))
ggplot(
    prop_df, aes(x = cellType_broad_hc_prop, y = prop_WM, color = position)
) +
    facet_wrap(~cellType_broad_hc) +
    geom_point() +
    labs(x = "Cell Type Prop.", y = "Prop. WM Spots", color = "Position") +
    theme_bw(base_size = 10)
dev.off()

session_info()
