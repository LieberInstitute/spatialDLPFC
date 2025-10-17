library(SpatialExperiment)
library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

dataset = 'combined'
num_factors = 6
specific_factor = 'Factor3'
fdr_cutoff = 0.05
z_cutoff_views = 1.3 # somewhat arbitrary
z_cutoff_genes = 1.96
max_genes = 2

model_path = here(
    'processed-data', 'MFA', 'models',
    sprintf('%s_%d.rds', dataset, num_factors)
)
sce_path = here(
    'processed-data', 'MFA', 'combined_rds', sprintf('%s.rds', dataset)
)
plot_dir = here('plots', 'MFA')

view_colors = c(
    DLPFC_sn_EndoMural = '#0C7602',
    DLPFC_sn_Micro = '#250A2F',
    DLPFC_sn_Oligo = '#A6E3BC',
    DLPFC_sn_Excit_L3 = '#032C5E',
    DLPFC_sn_Excit_L3.4.5 = '#0077B6',
    HPC_sn_CA1 = '#B8620C',
    HPC_sn_CA2.4 = '#751602',
    HPC_visium_CA2.4 = '#F8644E',
    Multi = 'grey30'
)

dir.create(file.path(plot_dir, 'GO'), showWarnings = FALSE)

model = load_model(model_path)

sce = readRDS(sce_path)$sce

gene_weights = get_geneweights(model = model, factor = specific_factor) |>
    left_join(
        rowData(sce) |>
            as.data.frame() |>
            dplyr::select(gene_id, gene_name),
        by = c("feature" = "gene_id")
    ) |>
    as_tibble() |>
    #   Z-score across genes within each view
    group_by(ctype) |>
    mutate(value = (value - mean(value)) / sd(value)) |>
    ungroup()

stopifnot(!any(is.na(gene_weights$gene_name)))

#   Only a subset of views will be used to find signature genes and
#   corresponding GO results, based on exceeding a z-score cutoff for R^2
#   for the factor of interest
view_r2 = model@cache$variance_explained$r2_per_factor$single_group[
    specific_factor,
]
view_r2_z = (view_r2 - mean(view_r2)) / sd(view_r2)
my_views = names(view_r2_z)[view_r2_z > z_cutoff_views]
stopifnot(length(my_views) > 0)
message(
    sprintf(
        'Only finding signature genes and GO results for views "%s"',
        paste(my_views, collapse = '", "')
    )
)

top_gene_weights = gene_weights |>
    mutate(abs_value = abs(value)) |>
    filter(ctype %in% my_views, abs_value > z_cutoff_genes)

stopifnot(nrow(top_gene_weights) > 0)

#   Only display up to max_genes per view and sign
select_genes = top_gene_weights |>
    group_by(ctype, sign(value)) |>
    arrange(desc(abs_value)) |>
    slice_head(n = max_genes) |>
    pull(feature)

################################################################################
#   Heatmap of top-weighted genes
################################################################################

#   Annotation of rows-- each gene is labeled by the view it is a signature gene
#   for and whether it has a positive weight
view_table = top_gene_weights |>
    filter(feature %in% select_genes) |>
    group_by(gene_name) |> 
    summarise(
        n = n(),
        positive_weight = Reduce('|', value > 0),
        view = paste0(ctype, collapse = ", ")
    ) |>
    mutate(
        view = factor(
            ifelse(n > 1, "Multi", view), levels = c("Multi", my_views)
        )
    ) |>
    dplyr::select(gene_name, positive_weight, view) |>
    column_to_rownames("gene_name") |>
    arrange(-positive_weight, view)

view_table_row = rowAnnotation(
    df = view_table,
    col = list(
        view = view_colors,
        positive_weight = c(`TRUE` = "grey80", `FALSE` = "grey20")
    )
)

## prep heatmap 
top_gw_value_matrix = gene_weights |>
    filter(feature %in% select_genes) |>
    dplyr::select(gene_name, ctype, value) |>
    pivot_wider(names_from = ctype, values_from = value) |>
    column_to_rownames("gene_name") |>
    as.matrix()

## weights across all clusters
pdf(
    file.path(plot_dir, sprintf('gene_weights_all_%s.pdf', specific_factor)),
    width = (ncol(top_gw_value_matrix) / 4) + 3, height = nrow(view_table) / 4
)
Heatmap(
    top_gw_value_matrix[rownames(view_table),],
    name = "feature\nweights\n(Z-score)",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = view_table_row
)
dev.off()

## weights for selected clusters only
pdf(
    file.path(plot_dir, sprintf('gene_weights_select_%s.pdf', specific_factor)),
    width = (length(my_views) / 4) + 3, height = nrow(view_table) / 4
)
Heatmap(
    top_gw_value_matrix[rownames(view_table), my_views],
    name = "feature\nweights\n(Z-score)",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = view_table_row
)
dev.off()

################################################################################
#   GO of top-weighted genes
################################################################################

do_go = function(gene_list, universe, plot_path) {
    for (ont_type in c("BP", "MF", "CC")) {
        go_obj = compareCluster(
            gene_list, fun = "enrichGO", universe = universe,
            OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH",
            pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE,
            keyType = "ENSEMBL"
        )

        if(!is.null(go_obj)) {
            go_obj@compareClusterResult = go_obj@compareClusterResult |>
                filter(p.adjust < fdr_cutoff, ONTOLOGY == ont_type)
            
            if(nrow(go_obj@compareClusterResult) > 0) {
                pdf(
                    sprintf(plot_path, ont_type),
                    height = as.integer(
                        min(11, 2 + nrow(go_obj@compareClusterResult) * 0.7)
                    ),
                )
                print(dotplot(go_obj, showCategory = 15))
                dev.off()
            }
        }
    }
}

#-------------------------------------------------------------------------------
#   One view at a time
#-------------------------------------------------------------------------------

for (this_view in my_views) {
    gene_list = list()
    gene_list[['up']] = top_gene_weights |>
        filter(value > 0, ctype == this_view) |>
        pull(feature)
    gene_list[['down']] = top_gene_weights |>
        filter(value < 0, ctype == this_view) |>
        pull(feature)

    message(
        sprintf(
            'For view "%s", using %d up and %d down signature genes',
            this_view, length(gene_list[['up']]), length(gene_list[['down']])
        )
    )

    #   Only consider genes measured in this view
    universe = gene_weights |>
        filter(ctype == this_view) |>
        pull(feature)

    do_go(
        gene_list,
        universe,
        file.path(
            plot_dir, 'GO', sprintf('%%s_%s_%s.pdf', specific_factor, this_view)
        )
    )
}

#-------------------------------------------------------------------------------
#   Selected views pooled together
#-------------------------------------------------------------------------------

gene_list = list()
gene_list[['up']] = top_gene_weights |>
    filter(value > 0) |>
    pull(feature)
gene_list[['down']] = top_gene_weights |>
    filter(value < 0) |>
    pull(feature)

#   Only consider genes measured in any of the selected views
universe = gene_weights |>
    filter(ctype %in% my_views) |>
    pull(feature) |>
    unique()

do_go(
    gene_list,
    universe,
    file.path(
        plot_dir, 'GO', sprintf('%%s_%s_together.pdf', specific_factor)
    )
)

session_info()
