
library("SpatialExperiment")
library("spatialLIBD")
library("tidyverse")
library("xlsx")
library("here")
library("sessioninfo")

## Set up plotting
plot_dir <- here("plots", "07_spatial_registration")
data_dir <- here("processed-data", "rdata", "spe", "07_spatial_registration")

## Load data 
# load(here("processed-data", "rdata","spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

## Load Registration Results
k_list <- c(9, 16, 28)
names(k_list) <- paste0("Sp", k_list) ## Use paper naming convention

bayesSpace_registration_fn <- map(k_list, ~here(data_dir, paste0("dlpfc_pseudobulked_bayesSpace_specific_Ts_k",.x,".Rdata")))
bayesSpace_registration <- lapply(bayesSpace_registration_fn, function(x) get(load(x)))

## Select t-stats from the registration enrichment data

registration_t_stats <- map2(bayesSpace_registration, k_list, function(data, k){
 t_stats <- sapply(data, function(x) {
    x$t[, 2, drop = FALSE]
  })
 
 rownames(t_stats) <- rownames(data[[1]]$t)
 colnames(t_stats) <- paste0("Sp", k, "D", colnames(t_stats))
 
 return(t_stats)
})

map(registration_t_stats, jaffelab::corner)

#### Calculate Correlation Matrix ####
## get layer data
layer_modeling_results <- fetch_data(type = "modeling_results")

#### Correlate with modeling results ####
cor_top100 <- map(registration_t_stats, ~layer_stat_cor(.x,
                                                        layer_modeling_results,
                                                        model_type = "enrichment",
                                                        reverse = FALSE,
                                                        top_n = 100))

save(cor_top100, file = here(data_dir, "layer_cor_top100.RDS"))

## Plot all for portability
pdf(here(plot_dir, "cor_top100_spatial_registration.pdf"))
map(cor_top100, layer_stat_cor_plot)
dev.off()

## Plot separately for illustrator
# map2(cor_top100, names(cor_top100), function(data, name){
#   
#   pdf(here(plot_dir, paste0("spatial_registration_plot_sn-",name,".pdf")))
#   layer_stat_cor_plot(data)
#   dev.off()
#   
# })

#### Annotate Layers ####
layer_anno <- map2(cor_top100, names(cor_top100), function(cor, name){
  anno <- annotate_registered_clusters(cor_stats_layer = cor,
                                       confidence_threshold = 0.25,
                                       cutoff_merge_ratio = 0.25)
  return(anno)
})

#### Annotate Cell Types by Layer ####
anno_abby <- data.frame(cluster = paste0("Sp9D",1:9), layer_abby = c("Vas", "L1", "L2/3", "L5", "L3", "WM", "L6A", "L4", "WM"))

layer_anno$Sp9 |> arrange(cluster) |> left_join(anno_abby)
# cluster layer_confidence layer_label layer_abby
# 1   Sp9D1             good          L1        Vas
# 2   Sp9D2             good          L1         L1
# 3   Sp9D3             good          L2       L2/3
# 4   Sp9D4             good          L5         L5
# 5   Sp9D5             good          L3         L3
# 6   Sp9D6             good          WM         WM
# 7   Sp9D7             good          L6        L6A
# 8   Sp9D8             good          L4         L4
# 9   Sp9D9             good       WM/L6         WM

## Add additonal annotaitons 
source(here("code","analysis", "12_spatial_registration_sn","utils.R"))

## layer_annotation is the reordered layer label - removes detail from the ordering process but helps group
layer_anno_all <- do.call("rbind", layer_anno) |>
  mutate(layer_annotation = fix_layer_order2(layer_label))

## Save for reference
write.csv(layer_anno_all, file = here(data_dir, "bayesSpace_layer_annotations.csv"))
# layer_anno_all <- read_csv(here(data_dir, "cellType_layer_annotations.csv")) 


#### Save Output to XLSX sheet ####
key <- data.frame(data = c("annotation", paste0("cor_", names(registration_t_stats))),
                  description = c("Annotations of baySpace Domains",
                                  paste0("Correlation values vs. manual annotation for ", names(registration_t_stats)))
                  )

## Clear file and write key
annotation_xlsx <- here(data_dir,"bayesSpace_layer_cor_annotations.xlsx")
write.xlsx(key, file=annotation_xlsx, sheetName="Key", append=FALSE, row.names=FALSE)

## write annotations
write.xlsx(layer_anno_all, file=annotation_xlsx, sheetName= paste0("annotation"), append=TRUE, row.names=FALSE)

## write correlations 
walk2(cor_top100, names(cor_top100),
        ~write.xlsx(t(.x), file=annotation_xlsx, sheetName= paste0("cor_", .y), append=TRUE, row.names=TRUE))

#### Explore Annotations ####
layer_anno_long <- layer_anno_all |>
  select(cluster, layer_label) |>
  pivot_longer(!cluster, names_to = "Annotation", values_to ="label") |>
  mutate(confidence = !grepl("\\*", label),
         layers = str_split(gsub("\\*","",label),"/"),
         Annotation = gsub("_label","", Annotation)) |>
  unnest_longer("layers") |>
  mutate(layer_short = ifelse(grepl("^[0-9]",layers), paste0("L",layers), layers),
         layer_long = gsub("L","Layer",layer_short))

layer_anno_long |> count(layer_long, layer_short)
layer_anno_long |> count(confidence)

#### Condensed Spatial Registration ####
cor_all <- do.call("rbind", cor_top100)

library("ComplexHeatmap")
library("jaffelab")

## match spatialLIBD color scale
theSeq <- seq(min(cor_all), max(cor_all), by = 0.01)
my.col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq))

## Add split for k
annotation_split <- gsub("9","09",ss(rownames(cor_all),"D"))

## Add annotations
anno_matrix <- layer_anno_long |> 
  mutate(fill = ifelse(confidence, "X","*")) |> 
  select(cluster, layer_long, fill) |> 
  pivot_wider(names_from = "layer_long", values_from = "fill", values_fill = "") |>
  column_to_rownames("cluster")

anno_matrix <- anno_matrix[rownames(cor_all), colnames(cor_all)]

corner(anno_matrix)

## all match
setdiff(colnames(anno_matrix),colnames(cor_all))
setdiff(rownames(anno_matrix),rownames(cor_all))

pdf(here(plot_dir,"spatial_registration_heatmap.pdf"))
Heatmap(cor_all,
        name = "Cor",
        col = my.col,
        row_split = annotation_split,
        rect_gp = gpar(col = "black", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(anno_matrix[i,j], x, y, gp = gpar(fontsize = 10))
        }
        )
dev.off()

pdf(here(plot_dir,"spatial_registration_heatmap-cluster.pdf"))
Heatmap(cor_all,
        name = "Cor",
        col = my.col,
        # row_split = annotation_split,
        rect_gp = gpar(col = "black", lwd = 1),
        cluster_rows = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(anno_matrix[i,j], x, y, gp = gpar(fontsize = 10))
        })
dev.off()

## Add annotation colors 


# sgejobs::job_single('02_cellType_correlation_annotation', create_shell = TRUE, memory = '5G', command = "Rscript 02_cellType_correlation_annotation.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
