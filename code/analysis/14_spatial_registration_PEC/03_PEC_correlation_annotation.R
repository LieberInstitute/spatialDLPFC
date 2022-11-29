library("spatialLIBD")
library("tidyverse")
library("xlsx")
library("jaffelab")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

## plot dir
plot_dir <- here("plots", "14_spatial_registration_PEC")
data_dir <- here("processed-data", "rdata", "spe", "14_spatial_registration_PEC")

#### Load Layer and k Data  ####
layer_modeling_results <- fetch_data(type = "modeling_results")

paths <- list(k09 = "modeling_results_BayesSpace_k09.Rdata", k16 = "modeling_results_BayesSpace_k16.Rdata")

modeling_results <- lapply(paths, function(x) {
    get(load(here("processed-data", "rdata", "spe", "07_layer_differential_expression", x)))
})

modeling_results <- c(list(layer = layer_modeling_results), modeling_results)
names(modeling_results)

cell_types <- c(
    "Astro", "Chandelier", "Endo", "L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP",
    "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Lamp5", "Lamp5 Lhx6", "Micro/PVM", "OPC",
    "Oligo", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "VLMC", "Vip"
)

names(cell_types) <- make.names(cell_types)


correlate_and_annotate <- function(dataset) {
    #### Load Registration Stats ####
    message(Sys.time(), " - ", dataset)
    registration_stats <- readRDS(here(
        "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
        paste0("registration_stats_", dataset, ".rds")
    ))

    registration_t_stats <- registration_stats[, grep("^t_stat", colnames(registration_stats))]
    colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

    # Fix Cell Types
    colnames(registration_t_stats) <- cell_types[colnames(registration_t_stats)]

    #### Correlate with modeling results ####
    cor_top100 <- map(modeling_results, ~ layer_stat_cor(registration_t_stats,
        .x,
        model_type = "enrichment",
        reverse = FALSE,
        top_n = 100
    ))

    # ## Plot
    # pdf(here(plot_dir, paste0("spatial_registration_plot_",dataset,".pdf")))
    # map(cor_top100, layer_stat_cor_plot)
    # dev.off()
    #
    #### Annotate Layers ####
    layer_anno <- map2(cor_top100, names(cor_top100), function(cor, name) {
        anno <- annotate_registered_clusters(
            cor_stats_layer = cor,
            confidence_threshold = 0.25,
            cutoff_merge_ratio = 0.25
        )
        colnames(anno) <- gsub("layer", name, colnames(anno))
        return(anno)
    })

    layer_anno <- reduce(layer_anno, left_join, by = "cluster")
    return(list(cor_top100 = cor_top100, layer_anno = layer_anno))
}

# datasets <- c("CMC", "DevBrain-snRNAseq", "IsoHuB", "SZBDMulti", "UCLA-ASD", "Urban-DLPFC")
datasets <- c("CMC", "DevBrain-snRNAseq", "IsoHuB", "UCLA-ASD", "Urban-DLPFC")
names(datasets) <- datasets

## Calculate correlations and annotations for each dataset
pe_correlation_annotation <- map(datasets, correlate_and_annotate)
save(pe_correlation_annotation, file = here(data_dir, "pe_correlation_annotation.Rdata"))

#### Save Output to XLSX sheet ####
key <- data.frame(
    data = c("annotation", paste0("cor_", names(modeling_results))),
    description = c(
        "Annotations of spatial registration",
        "Correlation values with manual layer annotations",
        "Correlation values with k9 domains",
        "Correlation values with k16 domains"
    )
)

## Clear file and write key
annotation_xlsx <- here(data_dir, "PE_spatial_annotations.xlsx")
write.xlsx(key, file = annotation_xlsx, sheetName = "Key", append = FALSE, row.names = FALSE)

## write annotations
walk2(
    pe_correlation_annotation, names(pe_correlation_annotation),
    ~ write.xlsx(.x$layer_anno, file = annotation_xlsx, sheetName = paste0("annotation_", .y), append = TRUE, row.names = FALSE)
)

## write correlations
walk2(pe_correlation_annotation, names(pe_correlation_annotation), function(data, name) {
    name <- paste0("cor_", name)
    # message(name)
    walk2(
        data$cor_top100, names(data$cor_top100),
        ~ write.xlsx(t(.x), file = annotation_xlsx, sheetName = paste0(name, "_", .y), append = TRUE, row.names = TRUE)
    )
})

#### Compare annotations for each cell type ####
## prep data
source(here("code", "analysis", "12_spatial_registration_sn", "utils.R"), echo = TRUE, max.deparse.length = 500)

layer_anno <- transpose(pe_correlation_annotation)$layer_anno

layer_anno_all <- do.call("rbind", layer_anno) |>
    rownames_to_column("Dataset") |>
    mutate(
        Dataset = gsub("\\.[0-9]+", "", Dataset),
        layer_label_order = gsub("\\*", "", fix_layer_order2(layer_label))
    )

layer_anno_long <- layer_anno_all |>
    select(Dataset, cluster, ends_with("label")) |>
    pivot_longer(!c(Dataset, cluster), names_to = "Annotation", values_to = "label") |>
    mutate(
        confidence = !grepl("\\*", label),
        layers = str_split(gsub("\\*", "", label), "/"),
        Annotation = gsub("_label", "", Annotation)
    ) |>
    unnest_longer("layers") |>
    # mutate(layers = ifelse(Annotation == "layer" & grepl("^[0-9]",layers), paste0("L",layers), layers))
    mutate(
        layer_long = ifelse(Annotation == "layer",
            ifelse(grepl("^[0-9]", layers), paste0("Layer", layers), gsub("L", "Layer", layers)),
            layers
        ),
        anno_confidence = ifelse(confidence, "high", "low")
    )

layer_anno_long |> count(layers)
layer_anno_long |> count(layer_long)
layer_anno_long |> count(confidence)
layer_anno_long |> count(Dataset, Annotation)
layer_anno_long |>
    filter(confidence) |>
    count(Dataset, Annotation)
layer_anno_long |>
    filter(Annotation == "layer", confidence) |>
    count(Dataset, Annotation, cluster)

layer_anno_long |>
    filter(Annotation == "layer") |>
    count(layers, layer_long)

#### Plot layer annotation ####
# highlight cell type layer annotaions
cell_type_anno <- tibble(cluster = cell_types) |>
    mutate(layer_label = ifelse(grepl("^L[0-9]", cluster), sub(" .*", "", cluster), NA)) |>
    filter(!is.na(layer_label)) |>
    mutate(layers = str_split(gsub("\\*", "", layer_label), "/")) |>
    unnest_longer("layers") |>
    mutate(
        layers = gsub("b", "", ifelse(grepl("^[0-9]", layers), paste0("L", layers), layers)),
        val = 1
    )

cell_type_anno |> count(layers)

## Match with bayesSpace spatial annotations for k09 k16 exploration
bayes_anno <- read.csv(file = here("processed-data", "rdata", "spe", "08_spatial_registration", "bayesSpace_layer_annotations.csv")) |>
    filter(bayesSpace %in% c("k09", "k16")) |>
    select(Annotation = bayesSpace, layer_long = cluster, layer_annotation, layer_combo)

spatial_layer_anno <- data.frame(
    Annotation = "layer",
    layer_long = c(paste0("Layer", 1:6), "WM"),
    layer_annotation = c(paste0("L", 1:6), "WM")
) |>
    mutate(layer_combo = layer_long) |>
    rbind(bayes_anno) |>
    mutate(layers = str_split(gsub("\\*", "", layer_annotation), "/")) |>
    unnest_longer("layers")

spatial_layer_anno

cell_type_anno_all <- spatial_layer_anno |>
    right_join(cell_type_anno |> select(layers, val, cluster))

layer_combo_factor <- function(x) {
    uniq <- unique(x)
    factor(x, levels = uniq[order(gsub(".* ~ ", "", gsub("ayer", "", uniq)))])
}

cell_type_anno_all$layer_combo <- layer_combo_factor(cell_type_anno_all$layer_combo)

anno_all_test <- cell_type_anno_all |>
    ggplot(aes(x = cluster, y = layer_combo)) +
    geom_tile(fill = "blue", alpha = 0.2) +
    facet_grid(Annotation ~ ., scales = "free_y", space = "free") +
    scale_y_discrete(limits=rev)

ggsave(anno_all_test, filename = here(plot_dir, "anno_all_test.png"))

## Add layer_combo to layer_anno_long to add bayes space details
layer_anno_long <- layer_anno_long |> left_join(spatial_layer_anno |> select(Annotation, layer_long, layer_combo))
layer_anno_long$layer_combo <- layer_combo_factor(layer_anno_long$layer_combo)

#### Dot plots ####
layer_anno_plot <- layer_anno_long |>
    filter(Annotation == "layer") |>
    ggplot(aes(x = cluster, y = layer_combo)) +
    geom_point(aes(color = Dataset, shape = anno_confidence), position = position_dodge(width = .8)) +
    geom_tile(data = cell_type_anno_all |> filter(Annotation == "layer"), fill = "blue", alpha = 0.2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_discrete(limits=rev)

ggsave(layer_anno_plot, filename = here(plot_dir, "PE_datasets_layer_annotation.png"), width = 10, height = 5)

## Filter to high confidence
layer_anno_plot_filter <- layer_anno_long |>
    filter(Annotation == "layer", confidence) |>
    ggplot(aes(x = cluster, y = layer_combo)) +
    geom_point(aes(color = Dataset), position = position_dodge(width = .8)) +
    geom_tile(data = cell_type_anno_all |> filter(Annotation == "layer"), fill = "blue", alpha = 0.2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_discrete(limits=rev)

ggsave(layer_anno_plot_filter, filename = here(plot_dir, "PE_datasets_layer_annotation_confident.png"), width = 10, height = 5)

## All annotations
label_anno_plot <- layer_anno_long |>
    ggplot(aes(x = cluster, y = layer_combo)) +
    geom_point(aes(color = Dataset, shape = anno_confidence), position = position_dodge(width = .8)) +
    geom_tile(data = cell_type_anno_all, fill = "blue", alpha = 0.2) +
    # facet_wrap(~Annotation, ncol = 1, scales = "free_y") +
    facet_grid(Annotation ~ ., scales = "free_y", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(y = "bayesSpace Domain & Annotation", x = "PsychEncode DLPFC Cell Types") +
    scale_y_discrete(limits=rev)

ggsave(label_anno_plot, filename = here(plot_dir, "PE_datasets_all_annotations.png"), height = 10)

label_anno_plot <- layer_anno_long |>
    filter(anno_confidence == "high") |>
    ggplot(aes(x = cluster, y = layer_combo)) +
    geom_point(aes(color = Dataset), position = position_dodge(width = .8)) +
    geom_tile(data = cell_type_anno_all, fill = "blue", alpha = 0.2) +
    # facet_wrap(~Annotation, ncol = 1, scales = "free_y") +
    facet_grid(Annotation ~ ., scales = "free_y", space = "free") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top"
    ) +
    labs(y = "bayesSpace Domain & Annotation", x = "PsychENCODE DLPFC Cell Types") +
    scale_y_discrete(limits=rev)

ggsave(label_anno_plot, filename = here(plot_dir, "PE_datasets_all_annotations_confident.png"), height = 7.25, width = 8)

#### Heatmap ####
cor_top100 <- transpose(pe_correlation_annotation)$cor_top100
corner(cor_top100_2$CMC)

## DEvBrain, IsoHub, and UrbanDLPFC missing
# [1] "Sst Chodl"
cor_top100_2 <- map2(cor_top100, names(cor_top100), function(cor_data, dataset) {
    cor_data <- map(cor_data, function(cd) {
        if (!setequal(rownames(cd), cell_types)) {
            return(setdiff(cell_types, rownames(cd)))
        }
        cd <- cd[cell_types, ]
        rownames(cd) <- paste0(dataset, "_", rownames(cd))
        return(cd)
    })
    return(do.call("cbind", cor_data))
})

map(cor_top100$DevBrain$k09, ~ .x[cell_types, ])

map(cor_top100$DevBrain, ~ setdiff(rownames(.x), cell_types))
map(cor_top100$DevBrain, ~ setdiff(cell_types, rownames(.x)))
map(cor_top100$DevBrain, ~ union(cell_types, rownames(.x)))

# sgejobs::job_single('03_correlate_spatial', create_shell = TRUE, memory = '25G', command = "Rscript 03_correlate_spatial.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
