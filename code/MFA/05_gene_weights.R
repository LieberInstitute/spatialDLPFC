library(SpatialExperiment)
library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(ComplexHeatmap)

dataset = 'DLPFC'
num_factors = 4
specific_factor = 'Factor2'
my_views = c('Excit_L4', 'Excit_L5', 'Excit_L5.6')

model_path = here(
    'processed-data', 'MFA', 'models',
    sprintf('%s_%d.rds', dataset, num_factors)
)
sce_path = here(
    'processed-data', 'MFA', 'combined_rds', sprintf('%s.rds', dataset)
)
plot_path = here(
    'plots', 'MFA',
    sprintf('gene_weights_%s_n%d_%s.pdf', dataset, num_factors, specific_factor)
)

model = load_model(model_path)

sce = readRDS(sce_path)$sce

gene_weights = get_geneweights(model = model, factor = specific_factor) |>
    left_join(
        rowData(sce) |>
            as.data.frame() |>
            select(gene_id, gene_name),
        by = c("feature" = "gene_id")
    ) |>
    as_tibble()

stopifnot(all(my_views %in% unique(gene_weights$ctype)))

top_gene_weights = gene_weights |>
    filter(ctype %in% my_views) |>
    mutate(abs_value = abs(value), weight_pos = value > 0) |>
    group_by(ctype, weight_pos) |>
    arrange(desc(abs_value)) |>
    slice_head(n = 5)

## prep heatmap 
top_gw_value_matrix = gene_weights |>
    filter(feature %in% top_gene_weights$feature) |>
    select(gene_name, ctype, value) |>
    pivot_wider(names_from = ctype, values_from = value) |>
    column_to_rownames("gene_name") |>
    as.matrix()

row_order = top_gene_weights |>
    group_by(gene_name) |> 
    summarise(max_value = max(value)) |>
    arrange(desc(max_value)) |>
    pull(gene_name)

top_gw_value_matrix = top_gw_value_matrix[row_order,]

pdf_width_all = (nrow(top_gw_value_matrix) / 8) + 3
pdf_width_select = (length(my_views) / 8) + 3
pdf_height = length(row_order) / 5

## weights across all clusters
pdf(plot_path, width = pdf_width_all, height = pdf_height)
Heatmap(
    top_gw_value_matrix,
    name = "feature\nweights",
    cluster_rows = FALSE,
    cluster_columns = FALSE
)
dev.off()

session_info()
