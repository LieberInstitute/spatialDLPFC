
library("spatialLIBD")
library("tidyverse")
library("xlsx")
library("jaffelab")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

## plot dir
plot_dir <- here("plots", "12_spatial_registration_sn", "05_velm_correlation_annotation")
data_dir <-
    here(
        "processed-data",
        "rdata",
        "spe",
        "12_spatial_registration_sn", 
        "05_velm_correlation_annotation"
    )

if(!dir.exists(data_dir)) dir.create(data_dir)

#### Load Layer and k Data  ####
layer_modeling_results <- fetch_data(type = "modeling_results")

paths <-
    list(k09 = "modeling_results_BayesSpace_k09.Rdata", k16 = "modeling_results_BayesSpace_k16.Rdata")

modeling_results <- lapply(paths, function(x) {
    get(load(
        here(
            "processed-data",
            "rdata",
            "spe",
            "07_layer_differential_expression",
            x
        )
    ))
})

modeling_results <-
    c(list(layer = layer_modeling_results), modeling_results)
names(modeling_results)


correlate_and_annotate <- function(registration_stats, make_cor_plot = FALSE) {
    #### Load Registration Stats ####

    registration_t_stats <-
        registration_stats[, grep("^t_stat", colnames(registration_stats))]
    colnames(registration_t_stats) <-
        gsub("^t_stat_", "", colnames(registration_t_stats))

    # # Fix Cell Types
    # colnames(registration_t_stats) <-
    #     cell_types[colnames(registration_t_stats)]

    #### Correlate with modeling results ####
    cor_top100 <-
        map(
            modeling_results,
            ~ layer_stat_cor(
                registration_t_stats,
                .x,
                model_type = "enrichment",
                reverse = FALSE,
                top_n = 100
            )
        )

    # ## Plot
    if (make_cor_plot) {
        pdf(here(
            plot_dir,
            paste0("spatial_registration_plot_", dataset, ".pdf")
        ))
        map(cor_top100, layer_stat_cor_plot, max = 1)
        dev.off()
    }

    #### Annotate Layers ####
    layer_anno <-
        map2(cor_top100, names(cor_top100), function(cor, name) {
            anno <- annotate_registered_clusters(
                cor_stats_layer = cor,
                confidence_threshold = 0.25,
                cutoff_merge_ratio = 0.10
            )
            colnames(anno) <- gsub("layer", name, colnames(anno))
            return(anno)
        })

    layer_anno <- reduce(layer_anno, left_join, by = "cluster")
    return(list(cor_top100 = cor_top100, layer_anno = layer_anno))
}

#### Load velm data ####
load(here("processed-data", "rdata", "spe", "12_spatial_registration_sn", "sn_velm_registration.Rdata"), verbose = TRUE)

velm_registration <- sn_hc_registration_dx
names(velm_registration)
# [1] "asd"     "control"

## extract enrichment data
velm_registration <- map(velm_registration, ~pluck(.x, "enrichment"))

## Calculate correlations and annotations for each dataset
velm_correlation_annotation <- map(velm_registration , correlate_and_annotate)

save(velm_correlation_annotation,
    file = here(data_dir, "Velmeshev_correlation_annotation.Rdata")
)

#### Save Output to XLSX sheet ####
key <- data.frame(
    data = c("annotation", paste0("cor_", names(modeling_results))),
    description = c(
        "Annotations of spatial registration",
        "Correlation values with manual layer annotations",
        "Correlation values with k09 domains",
        "Correlation values with k16 domains"
    )
)

## Clear file and write key
annotation_xlsx <- here(data_dir, "spatial_annotations_Velmeshev.xlsx")
write.xlsx(
    key,
    file = annotation_xlsx,
    sheetName = "Key",
    append = FALSE,
    row.names = FALSE
)

## write annotations
walk2(
    velm_correlation_annotation,
    names(velm_correlation_annotation),
    ~ write.xlsx(
        .x$layer_anno,
        file = annotation_xlsx,
        sheetName = paste0("annotation_", .y),
        append = TRUE,
        row.names = FALSE
    )
)

## write correlations
walk2(velm_correlation_annotation, names(velm_correlation_annotation), function(data, name) {
    name <- paste0("cor_", name)
    # message(name)
    walk2(
        data$cor_top100,
        names(data$cor_top100),
        ~ write.xlsx(
            t(.x),
            file = annotation_xlsx,
            sheetName = paste0(name, "_", .y),
            append = TRUE,
            row.names = TRUE
        )
    )
})

#### Compare annotations for each cell type ####

## Match with bayesSpace spatial annotations for k09 k16 exploration
bayes_layers <- get(load(here("processed-data", "rdata", "spe", "08_spatial_registration", "bayesSpace_layer_annotations.Rdata"))) |>
  select(Annotation = bayesSpace, layers = cluster, layer_combo) |>
  filter(Annotation %in% c("k09", "k16"))

## Add Layer levels + update factor
spatial_layer_anno <- tibble(
  Annotation = "layer",
  layers = c(paste0("Layer", 1:6), "WM"),
  layer_combo = layers
) |>
  bind_rows(bayes_layers) |>
  mutate(layer_combo = factor(layer_combo, levels = c(paste0("Layer", 1:6), "WM", levels(bayes_layers$layer_combo))),
         Annotation = factor(Annotation, levels = c("k16", "k09", "layer")))

## prep data
source(
    here("code", "analysis", "12_spatial_registration_sn", "utils.R"),
    echo = TRUE,
    max.deparse.length = 500
)

layer_anno <- transpose(velm_correlation_annotation)$layer_anno

layer_anno_all <- do.call("rbind", layer_anno) |>
    rownames_to_column("PrimaryDx") |>
    mutate(
        PrimaryDx = gsub("\\.[0-9]+", "", PrimaryDx),
        layer_label_order = gsub("\\*", "", fix_layer_order2(layer_label))
    ) 
# |> separate(cluster, into = c("cluster", "Dx"), sep = "_") |>
#   mutate(Dataset = ifelse(is.na(Dx), Dataset, paste0(Dataset, "_" ,Dx)))

layer_anno_long <- layer_anno_all |>
    select(PrimaryDx, cluster, ends_with("label")) |>
    pivot_longer(!c(PrimaryDx, cluster),
        names_to = "Annotation",
        values_to = "label"
    ) |>
    mutate(
        confidence = !grepl("\\*", label),
        layers = str_split(gsub("\\*", "", label), "/"),
        Annotation = gsub("_label", "", Annotation),
        cluster = factor(cluster)
    ) |>
    unnest_longer("layers") |>
    # mutate(layers = ifelse(Annotation == "layer" & grepl("^[0-9]",layers), paste0("L",layers), layers))
    mutate(
        layers = ifelse(Annotation == "layer",
            ifelse(
                grepl("^[0-9]", layers),
                paste0("Layer", layers),
                gsub("L","Layer",layers)
            ),
            layers
        ),
        # layers_long = gsub("Layers","L",layers),
        anno_confidence = ifelse(confidence, "high", "low")
    ) |>
  left_join(spatial_layer_anno)

#### Plot layer annotation ####

source(here("code", "analysis" ,"14_spatial_registration_PEC", "registration_dot_plot.R"))

# highlight cell type layer annotations

cell_types = unique(layer_anno_long$cluster)

ct_anno <- tibble(cell_type = cell_types, 
                  ct_cat = rep(c("Non-Neuronal", "Neuron"), c(6,11))
                  )

spatial_layer_anno_long <- spatial_layer_anno |> 
  mutate(ml = str_split(gsub("S.* ~ |ayer","", layer_combo), "/")) |> ## ml = "manual layer"
  unnest_longer("ml") |> 
  mutate(ml = ifelse(grepl("^[0-9]", ml), paste0("L", ml), ml))


cell_type_anno <- tibble(cluster = cell_types) |>
  mutate(layer_label = ifelse(grepl("^L[0-9]", cluster), sub("\\.", "/", cluster), NA)) |>
  filter(!is.na(layer_label)) |>
  mutate(ml = str_split(gsub("\\*", "", layer_label), "/")) |> ## ml = "manual layer"
  unnest_longer("ml") |>
  mutate(
    ml = gsub(".CC", "", ifelse(
      grepl("^[0-9]", ml), paste0("L", ml), ml)
    ),
    match = TRUE) |>
  left_join(spatial_layer_anno_long)

cell_type_anno |> count(ml)

## control dotplot
control_anno_long <- layer_anno_long |> 
  filter(PrimaryDx == "control") |>
  rename(cell_type = cluster) |>
  left_join(ct_anno)

## Not working.....
dotplot_control_datasets <- registration_dot_plot2(control_anno_long, color_by = "Dataset", ct_anno = cell_type_anno) 
ggsave(dotplot_control_datasets, filename = here(plot_dir, "registration_anno_dotplot_control.png"), width = 12)

#### Control only Heatmap ####

## Color setup
## match spatialLIBD color scale
theSeq <- seq(min(cor_all), max(cor_all), by = 0.01)
my.col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq))

# domain colors from polychrome
k_colors <- Polychrome::palette36.colors(28)
names(k_colors) <- c(1:28)

# ## cell color bar ## no cell colors for Velmdata
# cell_color_bar <- columnAnnotation(
#   " " = colnames(cor_all),
#   col = list(" " = cell_type_colors),
#   show_legend = FALSE
# )

## Add intermediate colors to layers
source(here("code", "analysis", "08_spatial_registration", "libd_intermediate_layer_colors.R"))
libd_intermediate_layer_colors <- c(spatialLIBD::libd_layer_colors, libd_intermediate_layer_colors)
names(libd_intermediate_layer_colors) <- gsub("ayer", "", names(libd_intermediate_layer_colors))

layer_anno_colors <- spatial_layer_anno |>
  separate(layer_combo, into = c(NA, "layer_color"), " ~ ", remove = FALSE) |>
  mutate(domain_color = as.integer(gsub("Sp[0-9]+D", "", layers)),
         layer_color = ifelse(is.na(layer_color),gsub("ayer","",layer_combo),layer_color)) |>
  arrange(Annotation) |>
  as.data.frame()

## build row annotation
bayes_color_bar <- rowAnnotation(
  df = layer_anno_colors |>
    select(domain_color, layer_color),
  col = list(
    domain_color = k_colors,
    layer_color = libd_intermediate_layer_colors
  ),
  show_legend = FALSE
)


#### Spatial Registration Heatmap ####
cor_top100 <- velm_correlation_annotation$control$cor_top100

## Build cor all
cor_top100$k09 <- cor_top100$k09[rownames(cor_top100$layer), ]
cor_top100$k16 <- cor_top100$k16[rownames(cor_top100$layer), ]
cor_top100$layer <- cor_top100$layer[, c(paste0("Layer", seq_len(6)), "WM")]

cor_all <- t(do.call("cbind", cor_top100))
corner(cor_all)

## reorder and rename to meet bayes annotations
cor_all <- cor_all[layer_anno_colors$layers,]
rownames(cor_all) <- layer_anno_colors$layer_combo

## Annotation matrix
anno_matrix <- layer_anno_long |>
  filter(PrimaryDx == "control") |>
  mutate(fill = ifelse(confidence, "X", "*")) |>
  select(cluster, layer_combo, fill) |>
  pivot_wider(names_from = "layer_combo", values_from = "fill", values_fill = "") 

## nned a better way to fix this...
setdiff(rownames(cor_all), colnames(anno_matrix))
# [1] "Sp09D09 ~ WM" "Sp16D12 ~ L6" "Sp16D11 ~ WM" "Sp16D13 ~ WM" ## missing need to fix!

anno_matrix <- anno_matrix |>
  mutate("Sp16D12 ~ L6" = "", ## add missing domains
         "Sp16D12 ~ L6" = "",
         "Sp16D11 ~ WM" = "",
         "Sp16D13 ~ WM" = "",
         "Sp09D09 ~ WM" = "") |>
  column_to_rownames("cluster") |>
  t()

fix_velm_ct <- function(ct){
  ct <- gsub("(\\d)\\.(\\d)","\\1/\\2",ct)
  gsub("\\.","-", ct)
}

colnames(cor_all) <- fix_velm_ct(colnames(cor_all))
colnames(anno_matrix) <- fix_velm_ct(colnames(anno_matrix))

setequal(colnames(cor_all), colnames(anno_matrix)) # good

# check rownames match
setequal(rownames(cor_all), rownames(anno_matrix))
setdiff(rownames(anno_matrix), rownames(cor_all))
setdiff(rownames(cor_all), rownames(anno_matrix))

## reorder
anno_matrix <- anno_matrix[rownames(cor_all), colnames(cor_all)]
corner(anno_matrix)
corner(cor_all)


## Ordered heatmap
pdf(here(plot_dir, "spatial_registration_Velemshev_heatmap.pdf"), height = 8, width = 10)
Heatmap(cor_all,
        name = "Cor",
        col = my.col,
        row_split = layer_anno_colors$Annotation,
        rect_gp = gpar(col = "black", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        right_annotation = bayes_color_bar,
        # bottom_annotation = cell_color_bar,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(anno_matrix[i, j], x, y, gp = gpar(fontsize = 10))
        }
)
dev.off()

#### K09 version ###
layer_anno_colors_k9 <- layer_anno_colors |> filter(Annotation == "k09")

bayes_color_bar_k9 <- rowAnnotation(
  df = layer_anno_colors_k9 |>
    select(domain_color, layer_color),
  col = list(
    domain_color = k_colors,
    layer_color = libd_intermediate_layer_colors
  ),
  show_legend = FALSE
)

#### K09 ####

cor_k9 <- cor_all[as.character(layer_anno_colors_k9$layer_combo),]
anno_matrix_k9 <- anno_matrix[rownames(cor_k9),]
## Ordered heatmap

# ## also case?
# case_k9 <- t(velm_correlation_annotation$asd$cor_top100$k09)
# 
# case_k9 <- case_k9[layer_anno_colors_k9$layers,colnames(cor_k9)]
# rownames(case_k9) <- layer_anno_colors_k9$layer_combo
# 
# diff_k09 = (cor_k9 - case_k9)/cor_k9

heatmap_k9 <- Heatmap(cor_k9,
                      name = "Cor",
                      col = tail(my.col,-1), ## if same num weird behavior ?
                      # row_split = layer_anno_colors$Annotation,
                      rect_gp = gpar(col = "black", lwd = 1),
                      cluster_rows = FALSE,
                      cluster_columns = TRUE,
                      right_annotation = bayes_color_bar_k9,
                      # bottom_annotation = cell_color_bar,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(anno_matrix_k9[i, j], x, y, gp = gpar(fontsize = 10))
                      }
)

co <- column_order(heatmap_k9)
# [1] 12 11 17 15 16  8  7 13 14  9 10  1  5  6  2  3  4
velm_col_order_k9 <- colnames(cor_k9)[co]


pdf(here(plot_dir, "spatial_registration_Velemshev_heatmap_k9.pdf"))
draw(heatmap_k9)
dev.off()

## save annotation data
save(anno_matrix_k9, velm_col_order_k9, file = here(data_dir, "Velmeshev_k9_annotation_details.Rdata"))


# Heatmap(case_k9,
#         name = "Cor - case",
#         col = tail(my.col,-1), ## if same num weird behavior ?
#         # row_split = layer_anno_colors$Annotation,
#         rect_gp = gpar(col = "black", lwd = 1),
#         cluster_rows = FALSE,
#         cluster_columns = TRUE,
#         right_annotation = bayes_color_bar_k9
#         # bottom_annotation = cell_color_bar,
#         # cell_fun = function(j, i, x, y, width, height, fill) {
#         #   grid.text(anno_matrix_reorder[i, j], x, y, gp = gpar(fontsize = 10))
#         # }
# )


# ## k09 control - case
# case_k9 <- t(velm_correlation_annotation$asd$cor_top100$k09)
# 
# case_k9 <- case_k9[layer_anno_colors_k9$layers,colnames(cor_k9)]
# rownames(case_k9) <- layer_anno_colors_k9$layer_combo
# 
# diff_k09 = (cor_k9 - case_k9)/cor_k9
# 
# ## Ordered heatmap
# pdf(here(plot_dir, "spatial_registration_Velemshev_heatmap_k9_diff.pdf"))
# Heatmap(diff_k09,
#         name = "Percent diff",
#         col = tail(my.col,-1), ## if same num weird behavior ?
#         # row_split = layer_anno_colors$Annotation,
#         rect_gp = gpar(col = "black", lwd = 1),
#         cluster_rows = FALSE,
#         cluster_columns = TRUE,
#         right_annotation = bayes_color_bar_k9
#         # bottom_annotation = cell_color_bar,
#         # cell_fun = function(j, i, x, y, width, height, fill) {
#         #   grid.text(anno_matrix_reorder[i, j], x, y, gp = gpar(fontsize = 10))
#         # }
# )
# dev.off()


# sgejobs::job_single('03_correlate_spatial', create_shell = TRUE, memory = '25G', command = "Rscript 03_correlate_spatial.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
