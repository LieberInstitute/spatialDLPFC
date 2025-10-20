#   Produce custom heatmaps showing a subset of views in the column annotation

library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

dataset = 'combined'
num_factors = 6
specific_factor = 'Factor3'
num_views_per_type = 5

model_path = here(
    'processed-data', 'MFA', 'models',
    sprintf('%s_%d.rds', dataset, num_factors)
)
plot_path = here(
    'plots', 'MFA', 'heatmap',
    sprintf('custom_views_%s_n%d_%s.pdf', dataset, num_factors, specific_factor)
)
out_path = here(
    'processed-data', 'MFA',
    sprintf('variance_explained_%s_n%d.csv', dataset, num_factors)
)

################################################################################
#   R2 column annotation
################################################################################

model = load_model(model_path)
factor_weights = model@cache$variance_explained$r2_per_factor$single_group

factor_weights |>
    as.data.frame() |>
    rownames_to_column('factor_num') |>
    write_csv(out_path)

#   For the specific factor of interest, show just the top few views per dataset
#   (the 4 combos of DLPFC/HPC Visium/snRNA-seq) by R^2
factor_weights_grouped = factor_weights |>
    as.data.frame() |>
    rownames_to_column('factor_num') |>
    as_tibble() |>
    pivot_longer(
        cols = -factor_num,
        names_to = 'view',
        values_to = 'r2'
    ) |>
    mutate(view_type = str_extract(view, '^(DLPFC|HPC)_(visium|sn)')) |>
    filter(factor_num == specific_factor) |>
    group_by(view_type) |>
    arrange(desc(r2)) |>
    slice_head(n = num_views_per_type)
stopifnot(all(factor_weights_grouped$view %in% colnames(factor_weights)))

#   Create a separate color scale for each view type
col_list = list()
weight_list = list()
for (this_view in unique(factor_weights_grouped$view_type)) {
    #   Grab just the weights for this view type and remove the view-type prefix
    weight_list[[this_view]] = factor_weights[
        specific_factor,
        (str_detect(colnames(factor_weights), paste0('^', this_view))) &
        (colnames(factor_weights) %in% factor_weights_grouped$view),
        drop = FALSE
    ]
    colnames(weight_list[[this_view]]) = str_replace(
        colnames(weight_list[[this_view]]), paste0('^', this_view, '_'), ''
    )

    #   Now order by descending R^2
    weight_list[[this_view]] = weight_list[[this_view]][
        ,
        order(weight_list[[this_view]][specific_factor,], decreasing = TRUE),
        drop = FALSE
    ]

    max_val = max(weight_list[[this_view]], na.rm = TRUE)
    col_list[[paste0('R2_', this_view)]] = colorRamp2(
        seq(0, min(100, max_val), length = 50),
        hcl.colors(50, "Oranges", rev = TRUE)
    )
}

column_ha = HeatmapAnnotation(
    R2_DLPFC_visium = weight_list[["DLPFC_visium"]],
    R2_DLPFC_sn = weight_list[["DLPFC_sn"]],
    R2_HPC_visium = weight_list[["HPC_visium"]],
    R2_HPC_sn = weight_list[["HPC_sn"]],
    gap = unit(2.5, "mm"),
    border = TRUE,
    col = col_list
)

################################################################################
#   Main heatmap body
################################################################################

factor_matrix = get_factors(model, factors = "all")$single_group[
    , specific_factor, drop = FALSE
]

max_fact = abs(max(factor_matrix))
col_fun_fact = colorRamp2(
    seq((max_fact + 0.5) * -1, max_fact + 0.5, length = 50),
    hcl.colors(50, "Green-Brown", rev = TRUE)
)

################################################################################
#   Putting it all together
################################################################################

scores_hmap = Heatmap(
    factor_matrix,
    name = "factor_scores",
    top_annotation = column_ha,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    border = TRUE,
    gap = unit(2.5, "mm"),
    col = col_fun_fact
)

pdf(
    plot_path,
    width = 2.7,
    height = as.integer(round(num_views_per_type * 0.75) + 5)
)
ComplexHeatmap::draw(scores_hmap, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

session_info()
