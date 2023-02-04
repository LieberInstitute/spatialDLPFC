library("spatialLIBD")
library("tidyverse")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

## plot dir
plot_dir <- here("plots", "14_spatial_registration_PEC", "07_PEC_annotation_Dx_plot")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <-
    here(
        "processed-data",
        "rdata",
        "spe",
        "14_spatial_registration_PEC"
    )

load(here(data_dir, "PEC_correlation_annotation_Dx.Rdata"), verbose = TRUE)
load(here(data_dir, "pec_cell_type_tb.Rdata"), verbose = TRUE)

#### Compare annotations for each cell type ####
## prep data
source(
    here("code", "analysis", "12_spatial_registration_sn", "utils.R"),
    echo = TRUE,
    max.deparse.length = 500
)

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
    mutate(
        layer_combo = factor(layer_combo, levels = c(paste0("Layer", 1:6), "WM", levels(bayes_layers$layer_combo))),
        Annotation = factor(Annotation, levels = c("k16", "k09", "layer"))
    )

## Extract Layer anno, pivot each longer
layer_anno <- transpose(pe_correlation_annotation)$layer_anno

head(layer_anno[[1]])
names(layer_anno)

layer_anno_long <- map(layer_anno, ~ .x |>
    separate(cluster, into = c("cluster", "PrimaryDx"), sep = "-") |>
    left_join(pec_cell_type_tb, by = "cluster") |>
    select(Dataset, PrimaryDx, cell_type, ct_cat, ends_with("label")) |>
    pivot_longer(!c(Dataset, PrimaryDx, cell_type, ct_cat),
        names_to = "Annotation",
        values_to = "label"
    ) |>
    mutate(
        confidence = !grepl("\\*", label),
        layers = str_split(gsub("\\*", "", label), "/"),
        Annotation = gsub("_label", "", Annotation)
    ) |>
    unnest_longer("layers") |>
    # mutate(layers = ifelse(Annotation == "layer" & grepl("^[0-9]",layers), paste0("L",layers), layers))
    mutate(
        layers = ifelse(Annotation == "layer",
            ifelse(
                grepl("^[0-9]", layers),
                paste0("Layer", layers),
                gsub("L", "Layer", layers)
            ),
            layers
        ),
        # layers_long = gsub("Layers","L",layers),
        anno_confidence = ifelse(confidence, "good", "poor"),
        Annotation = factor(Annotation, levels = c("k16", "k09", "layer"))
    ) |>
    left_join(spatial_layer_anno))


sort(unique(unlist(map(layer_anno_long, ~ .x$PrimaryDx))))
# [1] "Autism Spectrum Disorder" "Bipolar"                  "Bipolar Disorder"         "Control"
# [5] "Schizophrenia"            "Williams Syndrome"

dx_colors <- c(
    Control = "#2e282a",
    `Autism Spectrum Disorder` = "#76b041",
    `Bipolar` = "#17bebb",
    `Bipolar Disorder` = "#17bebb",
    Schizophrenia = "#e4572e",
    `Williams Syndrome` = "#CC9C00"
)

## colors
load(here(data_dir, "pec_dataset_colors.Rdata"), verbose = TRUE)


#### Plot layer annotation ####
## prep cell type annotation
spatial_layer_anno_long <- spatial_layer_anno |>
    mutate(ml = str_split(gsub("S.* ~ |ayer", "", layer_combo), "/")) |> ## ml = "manual layer"
    unnest_longer("ml") |>
    mutate(ml = ifelse(grepl("^[0-9]", ml), paste0("L", ml), ml))

cell_type_anno <- pec_cell_type_tb |>
    mutate(layer_label = ifelse(grepl("^L[0-9]", cell_type), sub(" .*", "", cell_type), NA)) |>
    filter(!is.na(layer_label)) |>
    mutate(ml = str_split(gsub("\\*", "", layer_label), "/")) |> ## ml = "manual layer"
    unnest_longer("ml") |>
    mutate(
        ml = gsub("b", "", ifelse(
            grepl("^[0-9]", ml), paste0("L", ml), ml
        )),
        Match = TRUE
    ) |>
    left_join(spatial_layer_anno_long) |>
    select(cell_type, layer_combo, Match)


spatial_layer_anno |>
    mutate(ml = str_split(gsub("S.* ~ |ayer", "", layer_combo), "/")) |> ## ml = "manual layer"
    unnest_longer("ml")

source("registration_dot_plot.R")

## Plot Dx
# walk2(layer_anno_long, names(layer_anno_long), function(data, name){
#
#   temp_dx_colors <- dx_colors[unique(data$PrimaryDx)]
#
#   dotplot <- registration_dot_plot2(data, ct_anno = cell_type_anno) +
#     labs(title = name) +
#     scale_color_manual(values = temp_dx_colors)+
#     labs(x = "PsychENCODE DLPFC Cell Types", y = "Spatial Domains")
#
#     ggsave(dotplot, filename = here(plot_dir, paste0("registration_anno_dotplot_dx_",name,".png")), width = 10)
#
#     # dotplot_conf <- registration_dot_plot2(data, ct_anno = cell_type_anno, conf_only = TRUE) +
#     #   labs(title = name) +
#     #   scale_color_manual(values = temp_dx_colors)
#     # ggsave(dotplot_conf, filename = here(plot_dir, paste0("registration_anno_dotplot_dx_conf_",name,".png")), width = 10)
#     #
# })

#### Combine Control, plot to compare Datasets ####
control_anno_long <- do.call("rbind", map2(
    layer_anno_long, names(layer_anno_long),
    ~ .x |>
        mutate(Dataset = .y) |>
        filter(PrimaryDx == "Control")
))

dotplot_control_datasets <- registration_dot_plot2(control_anno_long, color_by = "Dataset", ct_anno = cell_type_anno) +
    # theme(legend.position = "bottom")+
    scale_color_manual(values = pec_dataset_colors) +
    labs(x = "PsychENCODE DLPFC Cell Types", y = "Spatial Domains")

ggsave(dotplot_control_datasets, filename = here(plot_dir, "registration_anno_dotplot_control.png"), width = 13)
ggsave(dotplot_control_datasets, filename = here(plot_dir, "registration_anno_dotplot_control.pdf"), width = 13)

# #### ASD Astrocyte Subset ####
#
# astro_anno <- do.call("rbind", map2(layer_anno_long[c("DevBrain-snRNAseq", "UCLA-ASD")], c("DevBrain", "UCLA-ASD"),
#      ~.x |>
#        filter(cell_type == "Astro") |>
#        mutate(cell_type = paste0(cell_type, "\n", .y))))
#
# asd_dx_colors <- dx_colors[unique(astro_anno$PrimaryDx)]
#
# astro_asd_datasets <- registration_dot_plot2(astro_anno, color_by = "PrimaryDx") +
#   scale_color_manual(values = asd_dx_colors) +
#   # theme(axis.text.x = element_text(angle = 0))
#   theme_bw()
#
# ggsave(astro_asd_datasets, filename = here(plot_dir, "registration_anno_dotplot_astro_asd.png"), height = 5)
# ggsave(astro_asd_datasets, filename = here(plot_dir, "registration_anno_dotplot_astro_asd.pdf"), height = 5)
#
# #### ASD Astrocyte Heatmap ####
# combine_cor <- function(cor_list){
#   ## combine Annotations by column
#   # cor_list <- map_depth(cor_list, 2, ~do.call("cbind", .x))
#
#   # just want k16
#   cor_list <- map_depth(cor_list, 2, "k16")
#
#   ## combine Dx by row
#
#   cor_list2 <- map(cor_list,  function(data){
#     cor <- do.call("rbind", data)
#     rownames(cor) <- gsub("-", "_", rownames(cor))
#     return(cor)
#   })
#
#   cor <- do.call("rbind", cor_list2)
#   rownames(cor) <- paste0(rep(names(cor_list2), map_int(cor_list2, nrow)), "_", rownames(cor))
#   return(cor)
#   }
#
# cor_asd <- transpose(pe_correlation_annotation)$cor_top100[c("DevBrain-snRNAseq", "UCLA-ASD")]
# cor_asd <- combine_cor(cor_asd)
# dim(cor_asd)
#
# k16_annotation <- bayes_layers |> filter(Annotation == "k16")
# asd_astro <- t(cor_asd[grepl("Astro",rownames(cor_asd)),k16_annotation$layers])
# # asd_astro <- t(cor_asd[,k16_annotation$layers]) ## All
# rownames(asd_astro) <- k16_annotation$layer_combo
#
#
# #### Prep annotation heatmap ####
# asd_dataset_annotation <- data.frame(cn = colnames(asd_astro)) |>
#   separate(cn, into = c("Dataset", "cluster", "PrimaryDx"), sep = "_", extra = "merge") |>
#   left_join(pec_cell_type_tb) |>
#   select(-cluster)
#
# rownames(asd_dataset_annotation) <- colnames(asd_astro)
#
# ## arrange by cell types
# asd_dataset_annotation <- asd_dataset_annotation |>
#   arrange(cell_type, ct_cat)
#
# asd_astro <- asd_astro[,rownames(asd_dataset_annotation)]
#
#
# ## Annotaiton matrix for subset
# anno_matrix <- do.call("cbind", map2(layer_anno_long[c("DevBrain-snRNAseq", "UCLA-ASD")], c("DevBrain-snRNAseq", "UCLA-ASD"),
#                     ~.x |>
#                       mutate(fill = ifelse(confidence, "X", "*"),
#                              name = paste0(.y, ".", PrimaryDx, ".",cell_type)) |>
#                       filter(Annotation == "k16") |>
#                       select(name, layer_combo, fill) |>
#                       pivot_wider(names_from = "layer_combo", values_from = "fill", values_fill = "") |>
#                       filter(!is.na(name)) |>
#                       column_to_rownames("name") |>
#                       t()))
#
# # anno_matrix <- anno_matrix[c("Sp16D02 ~ L1", "Sp16D14 ~ L1"),]
#
# ## colors
# cell_color_bar <- columnAnnotation(
#   df = asd_dataset_annotation,
#   col = list(PrimaryDx = asd_dx_colors,
#              Dataset = dataset_colors,
#              ct_cat = c(Excit = "blue", Inhib = "red", `Non-neuronal` = "orange"))
# )
#
# layer_anno_colors <- k16_annotation |>
#   mutate(domain_color = as.integer(gsub("Sp[0-9]+D", "", layers))) |>
#   separate(layer_combo, into = c(NA, "layer_color"), " ~ ", remove = FALSE) |>
#   select(Annotation, layer_combo, layers, domain_color, layer_color)
#
# # domain colors from polychrome
# k_colors <- Polychrome::palette36.colors(28)
# names(k_colors) <- c(1:28)
#
# source(here("code", "analysis", "08_spatial_registration", "libd_intermediate_layer_colors.R"))
# libd_intermediate_layer_colors <- c(spatialLIBD::libd_layer_colors, libd_intermediate_layer_colors)
# names(libd_intermediate_layer_colors) <- gsub("ayer", "", names(libd_intermediate_layer_colors))
#
#
# bayes_color_bar <- rowAnnotation(
#   df = layer_anno_colors |>
#     select(domain_color, layer_color),
#   col = list(
#     domain_color = k_colors,
#     layer_color = libd_intermediate_layer_colors
#   ),
#   show_legend = FALSE
# )
#
# theSeq <- seq(min(asd_astro), max(asd_astro), by = 0.01)
# my.col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq))
#
#
# pdf(here(plot_dir, "All_asd_heatmap.pdf"), width = 12, height = 11)
#
# Heatmap(asd_astro,
#         name = "Cor",
#         col= my.col,
#         show_column_names = FALSE,
#         bottom_annotation = cell_color_bar,
#         right_annotation = bayes_color_bar,
#         cluster_rows = FALSE
# )
#
# Heatmap(asd_astro,
#         name = "Cor",
#         col= my.col,
#         show_column_names = FALSE,
#         bottom_annotation = cell_color_bar,
#         right_annotation = bayes_color_bar,
#         cluster_rows = FALSE,
#         cluster_columns = FALSE)
# dev.off()
#
#
# pdf(here(plot_dir, "Astro_asd_k16.pdf"), width = 8.5, height = 11)
#
# Heatmap(asd_astro,
#         name = "Cor",
#         col= my.col,
#         show_column_names = FALSE,
#         bottom_annotation = cell_color_bar,
#         right_annotation = bayes_color_bar,
#         cluster_rows = FALSE,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(round(asd_astro, 2)[i, j], x, y, gp = gpar(fontsize = 10))
#         }
# )
#
# dev.off()
#
#
# #### Annotation breakdown ####
# apply(asd_astro, 2, rank)
#
# asd_test <- asd_astro[c("Sp16D02 ~ L1", "Sp16D14 ~ L1"),]
# apply(asd_test, 2, rank)
#
# annotate_registered_clusters(t(asd_astro),
#                              confidence_threshold = 0.25,
#                              cutoff_merge_ratio = 0.10)
#
#
# ## confidence_threshold = 0.25,
# # cutoff_merge_ratio = 0.10
# asd_astro[c("Sp16D02 ~ L1", "Sp16D14 ~ L1"),]
# # DevBrain-snRNAseq.Autism Spectrum Disorder.Astro DevBrain-snRNAseq.Control.Astro
# # Sp16D02 ~ L1                                        0.5373625                       0.5215637
# # Sp16D14 ~ L1                                        0.5981215                       0.5696782
# # DevBrain-snRNAseq.Williams Syndrome.Astro UCLA-ASD.Autism Spectrum Disorder.Astro UCLA-ASD.Control.Astro
# # Sp16D02 ~ L1                                 0.5245672                               0.5535913              0.5535345
# # Sp16D14 ~ L1                                 0.5779144                               0.6091387              0.5963688
#
#
# confidence <- apply(asd_astro, 1, max) > .025
#
# (asd_test["Sp16D14 ~ L1",]- asd_test["Sp16D02 ~ L1",])/asd_test["Sp16D14 ~ L1",]
# # DevBrain-snRNAseq.Autism Spectrum Disorder.Astro                  DevBrain-snRNAseq.Control.Astro
# # 0.10158305                                       0.08445910
# # DevBrain-snRNAseq.Williams Syndrome.Astro          UCLA-ASD.Autism Spectrum Disorder.Astro
# # 0.09230992                                       0.09119004
# # UCLA-ASD.Control.Astro
# # 0.07182516
#
# (asd_test["Sp16D14 ~ L1",]- asd_test["Sp16D02 ~ L1",])/asd_test["Sp16D14 ~ L1",] > .1 ## only DevBrain ASD ??
# # DevBrain-snRNAseq_Astro_Autism Spectrum Disorder                  DevBrain-snRNAseq_Astro_Control
# # TRUE                                            FALSE
# # DevBrain-snRNAseq_Astro_Williams Syndrome          UCLA-ASD_Astro_Autism Spectrum Disorder
# # FALSE                                            FALSE
# # UCLA-ASD_Astro_Control
# # FALSE
#
#


# sgejobs::job_single('07_PEC_annotation_Dx_plots', create_shell = TRUE, memory = '25G', command = "Rscript 07_PEC_annotation_Dx_plots.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
