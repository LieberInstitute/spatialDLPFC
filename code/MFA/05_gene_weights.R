library(SpatialExperiment)
library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)

dataset = 'combined'
num_factors = 5
specific_factor = 'Factor1'

model_path = here(
    'processed-data', 'MFA', 'models',
    sprintf('%s_%d.rds', dataset, num_factors)
)
sce_path = here(
    'processed-data', 'MFA', 'combined_rds', sprintf('%s.rds', dataset)
)

model = load_model(model_path)

sce = readRDS(sce_path)$sce

top_gene_weights <- gene_weights$Factor3 |>
    mutate(abs_value = abs(value),
           weight_pos = value > 0,
           datatype = ifelse(grepl("_Sp", ctype), "Visium", "snRNA-seq")) |>
    filter(ctype %in% my_views) |>
    group_by(ctype,weight_pos) |>
    arrange(-abs_value) |>
    dplyr::slice(1:5)


## prep heatmap 
top_gw_value_matrix <- gene_weights$Factor3 |>
    filter(feature %in% top_gene_weights$feature) |>
    select(gene_name, ctype, value) |>
    pivot_wider(names_from = ctype, values_from = value) |>
    column_to_rownames("gene_name") |>
    as.matrix()

col_order <- cluster_levels[cluster_levels %in% colnames(top_gw_value_matrix)]

# row_order <- order(top_gw_value_matrix[,"Oligo"])

row_order <- top_gene_weights |>
    group_by(gene_name) |> 
    summarise(max_value = max(value)) |>
    arrange(-max_value) |>
    pull(gene_name)

top_gw_value_matrix <- top_gw_value_matrix[row_order,col_order]

dim(top_gw_value_matrix)
top_gw_value_matrix[1:5, 1:5]


pdf_width_all <- (nrow(top_gw_value_matrix)/8) + 3
pdf_width_select <- (length(my_views)/8) + 3

pdf_height <- length(row_order)/5

## weights across all clusters
pdf(here(plot_dir, sprintf("MOFA_gene_weight_heatmap_%s_all.pdf", opt$datatype)), width = pdf_width_all, height = pdf_height)
Heatmap(top_gw_value_matrix,
        name = "feature\nweights",
        cluster_rows = FALSE,
        cluster_columns = FALSE)
dev.off()
