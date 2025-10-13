#   Produce custom heatmaps showing a subset of views in the column annotation

library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

task_id = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (task_id == 1) {
    specific_factor = 'Factor4'
} else {
    specific_factor = 'Factor2'
}

dataset = 'combined'
num_factors = 6
num_views_per_type = 5

model_path = here(
    'processed-data', 'MFA', 'models',
    sprintf('%s_%d.rds', dataset, num_factors)
)
plot_path = here(
    'plots', 'MFA', 'heatmap',
    sprintf('custom_views_%s_n%d_%s.pdf', dataset, num_factors, specific_factor)
)

################################################################################
#   R2 column annotation
################################################################################

model = load_model(model_path)
factor_weights = model@cache$variance_explained$r2_per_factor$single_group

#   For the specific factor of interest, show just the top few views per dataset
#   (the 4 combos of DLPFC/HPC Visium/snRNA-seq) by R^2
highlighted_views = factor_weights |>
    as.data.frame() |>
    rownames_to_column('factor_num') |>
    as_tibble() |>
    pivot_longer(
        cols = -factor_num,
        names_to = 'view',
        values_to = 'r2'
    ) |>
    mutate(
        view_type = case_when(
            grepl('^HPC_sn_', view) ~ 'HPC_sn',
            grepl('^HPC_visium', view) ~ 'HPC_visium',
            grepl('^DLPFC_Sp09', view) ~ 'DLPFC_visium',
            TRUE ~ 'DLPFC_sn'
        )
    ) |>
    filter(factor_num == specific_factor) |>
    group_by(view_type) |>
    arrange(desc(r2)) |>
    slice_head(n = num_views_per_type) |>
    pull(view)

stopifnot(all(highlighted_views %in% colnames(factor_weights)))
factor_weights = factor_weights[, highlighted_views]

col_fun_r2 = colorRamp2(
    seq(0, min(100, max(factor_weights, na.rm = TRUE) + 5), length = 50),
    hcl.colors(50, "Oranges", rev = TRUE)
)

################################################################################
#   Donor-level covariates column annotation
################################################################################

#   Test association of various donor-level covariates with each factor
assoc_list = list()
assoc_list[['sex']] = get_associations(
    model = model,
    metadata = samples_metadata(model),
    sample_id_column = "sample",
    test_variable = "sex",
    test_type = "categorical",
    group = FALSE
)
assoc_list[['age']] = get_associations(
    model = model,
    metadata = samples_metadata(model),
    sample_id_column = "sample",
    test_variable = "age",
    test_type = "continuous",
    group = FALSE
)

#   Show unadjusted p-value for consistency with other plots
for (var_name in names(assoc_list)) {
    assoc_list[[var_name]]$adj_pvalue = assoc_list[[var_name]]$p.value
}

assoc_pvals = assoc_list |>
    tibble::enframe(name = "test") |>
    tidyr::unnest(c(value)) |>
    dplyr::mutate(log_adjpval = -log10(.data$adj_pvalue)) |>
    dplyr::select(test, Factor, log_adjpval) |>
    tidyr::pivot_wider(names_from = test, values_from = log_adjpval) |>
    dplyr::select(-Factor) |>
    as.matrix()

col_fun_assoc = colorRamp2(
    seq(0, max(assoc_pvals) + 0.5, length = 20),
    hcl.colors(20, "Purples", rev = TRUE)
)

column_ha = HeatmapAnnotation(
    "R2" = factor_weights,
    "pvalue" = assoc_pvals,
    gap = unit(2.5, "mm"),
    border = TRUE,
    col = list(R2 = col_fun_r2, pvalue = col_fun_assoc)
)
