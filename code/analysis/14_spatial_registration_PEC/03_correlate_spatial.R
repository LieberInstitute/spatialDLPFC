
# library(SingleCellExperiment)
library("SpatialExperiment")
library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

## plot dir
plot_dir <- here("plots","14_spatial_registration_PEC")

#### Load Layer and k Data  ####
layer_modeling_results <- fetch_data(type = "modeling_results")

paths <- list(k9 = "parsed_modeling_results_k9.Rdata", k16 = "parsed_modeling_results_k16.Rdata")

modeling_results <- lapply(paths, function(x) 
  get(load(here("processed-data","rdata","spe","08_layer_differential_expression",x))))

modeling_results <- c(list(layer = layer_modeling_results), modeling_results)
names(modeling_results)

cell_types <- c("Astro","Chandelier","Endo","L2/3 IT","L4 IT","L5 ET","L5 IT","L5/6 NP",   
                "L6 CT","L6 IT","L6 IT Car3","L6b","Lamp5","Lamp5 Lhx6","Micro/PVM","OPC",     
                "Oligo","Pax6","Pvalb","Sncg","Sst","Sst Chodl","VLMC","Vip")

names(cell_types) <- make.names(cell_types)


correlate_and_annotate <- function(dataset){
  #### Load Registration Stats ####
  message(Sys.time(), " - " ,dataset)
  registration_stats <- readRDS(here("processed-data","rdata","spe","14_spatial_registration_PEC",
                                     paste0("registration_stats_",dataset,".rds")))
  
  registration_t_stats <- registration_stats[, grep("^t_stat", colnames(registration_stats))]
  colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))
  
  # Fix Cell Types
  colnames(registration_t_stats) <- cell_types[colnames(registration_t_stats)]
  
  #### Correlate with modeling results ####
  cor_top100 <- map(modeling_results, ~layer_stat_cor(registration_t_stats,
                                                      .x,
                                                      model_type = "enrichment",
                                                      reverse = FALSE,
                                                      top_n = 100))
  
  ## Plot
  pdf(here(plot_dir, paste0("spatial_registration_plot_",dataset,".pdf")))
  map(cor_top100, layer_stat_cor_plot)
  dev.off()
  # 
  #### Annotate Layers ####
  layer_anno <- map2(cor_top100, names(cor_top100), function(cor, name){
    anno <- annotate_registered_clusters(cor_stats_layer = cor,
                                         confidence_threshold = 0.25,
                                         cutoff_merge_ratio = 0.25)
    colnames(anno) <- gsub("layer",name, colnames(anno))
    return(anno)
  })
  
  layer_anno <- reduce(layer_anno, left_join, by = "cluster")
  return(list(cor_top100 = cor_top100, layer_anno = layer_anno))
}

datasets <- c("DevBrain", "IsoHuB","CMC","UCLA-ASD")
# datasets <- c("DevBrain", "IsoHuB","CMC")
names(datasets) <- datasets

## Caluclate correlatiosn and annotaions for each dataset
pe_correlation_annotation <- map(datasets, correlate_and_annotate)

#### Save Output to XLSX sheet ####
data_dir <- here("processed-data","rdata","spe","14_spatial_registration_PEC")
names(pe_correlation_annotation$DevBrain$cor_top100)

library("xlsx")

key <- data.frame(data = c("annotation", paste0("cor_", names(modeling_results))),
                  description = c("Annotations of spatial registration",
                                  "Correlation values with manual layer annotations",
                                  "Correlation values with k9 domains",
                                  "Correlation values with k16 domains"))

## Clear file and write key
annotation_xlsx <- here(data_dir,"PE_spatial_annotations.xlsx")
write.xlsx(key, file=annotation_xlsx, sheetName="Key", append=FALSE, row.names=FALSE)

## write annotations
walk2(pe_correlation_annotation, names(pe_correlation_annotation), 
      ~write.xlsx(.x$layer_anno, file=annotation_xlsx, sheetName= paste0("annotation_", .y), append=TRUE, row.names=FALSE))

## write correlations 
walk2(pe_correlation_annotation, names(pe_correlation_annotation), function(data, name){
  
  name <- paste0("cor_",name)
  # message(name)
  walk2(data$cor_top100, names(data$cor_top100),
        ~write.xlsx(t(.x), file=annotation_xlsx, sheetName= paste0(name, "_", .y), append=TRUE, row.names=TRUE))
  
})

#### Compare annotations for each cell type ####
## prep data
source(here("code","analysis", "12_spatial_registration_sn","utils.R"))

layer_anno <- transpose(pe_correlation_annotation)$layer_anno
layer_anno_all <- do.call("rbind", layer_anno) |>
  rownames_to_column("Dataset") |>
  mutate(Dataset = gsub("\\.[0-9]+","",Dataset),
         layer_label_order = gsub("\\*","", fix_layer_order2(layer_label)))

layer_anno_long <- layer_anno_all |>
  select(Dataset, cluster,layer_confidence, layer_label) |>
  mutate(layers = str_split(gsub("\\*","",layer_label),"/")) |>
  unnest_longer("layers") |>
  mutate(layers = ifelse(grepl("^[0-9]",layers), paste0("L",layers), layers))

layer_anno_long |> count(layers)

cell_type_anno <- tibble(cluster = cell_types) |>
  mutate(layer_label = ifelse(grepl("^L[0-9]", cluster),sub(" .*", "",cluster),NA)) |>
  filter(!is.na(layer_label))|>
  mutate(layers = str_split(gsub("\\*","",layer_label),"/")) |>
  unnest_longer("layers")|>
  mutate(layers = gsub("b","",ifelse(grepl("^[0-9]",layers), paste0("L",layers), layers)),
         val = 1)

cell_type_anno |> count(layers)

## Plot 
layer_anno_plot <- ggplot(layer_anno_long, aes(x = cluster, y = layers)) +
  geom_jitter(aes(color = Dataset, shape = layer_confidence),width = 0.1, height = 0.1) +
  geom_tile(data = cell_type_anno, aes(x = cluster, y = layers), fill = "blue", alpha = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(layer_anno_plot, filename = here(plot_dir, "dataset_layer_annotation.png"))


layer_anno_plot2 <- ggplot(cell_type_anno) +
  # geom_jitter(aes(color = Dataset, shape = layer_confidence),width = 0.1, height = 0.1) +
  geom_tile(aes(x = cluster, y = layers, fill= val)) 

ggsave(layer_anno_plot2, filename = here(plot_dir, "dataset_layer_annotation2.png"))

