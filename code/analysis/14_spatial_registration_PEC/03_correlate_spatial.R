
# library(SingleCellExperiment)
library("SpatialExperiment")
library("spatialLIBD")
library("purrr")
library("here")
library("sessioninfo")

dataset <- "DevBrain"

#### Load Registration Stats ####
registration_stats <- readRDS(here("processed-data","rdata","spe","14_spatial_registration_PEC",
                           paste0("registration_stats_",dataset,".rds")))

registration_t_stats <- registration_stats[, grep("^t_stat", colnames(registration_stats))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

# Fix Cell Types
cell_types <- c("Astro","Chandelier","Endo","L2/3 IT","L4 IT","L5 ET","L5 IT","L5/6 NP",   
                "L6 CT","L6 IT","L6 IT Car3","L6b","Lamp5","Lamp5 Lhx6","Micro/PVM","OPC",     
                "Oligo","Pax6","Pvalb","Sncg","Sst","Sst Chodl","VLMC","Vip")

names(cell_types) <- make.names(cell_types)
colnames(registration_t_stats) <- cell_types[colnames(registration_t_stats)]

#### Load Layer and k Data  ####
layer_modeling_results <- fetch_data(type = "modeling_results")

paths <- list(k9 = "parsed_modeling_results_k9.Rdata", k16 = "parsed_modeling_results_k16.Rdata")

modeling_results <- lapply(paths, function(x) 
  get(load(here("processed-data","rdata","spe","08_layer_differential_expression",x), verbose = TRUE)))

modeling_results <- c(list(layer = layer_modeling_results), modeling_results)
names(modeling_results)

#### Correlate with modeling results ####
cor_top100 <- map(modeling_results, ~layer_stat_cor(registration_t_stats,
                                             .x,
                                             model_type = "enrichment",
                                             reverse = FALSE,
                                             top_n = 100))

#### Annotate Layers ####
layer_anno <- map(cor_top100, ~annotate_registered_clusters(cor_stats_layer = .x,
                                                            confidence_threshold = 0.25,
                                                            cutoff_merge_ratio = 0.25))

## Save data

#### Plot ####
plot_dir <- here("plots","14_spatial_registration_PEC")

pdf(here(plot_dir, paste0("spatial_registration_plot_",dataset,".pdf")))
map(cor, layer_stat_cor_plot)
dev.off()




