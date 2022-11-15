
library("SpatialExperiment")
library("spatialLIBD")
library("tidyverse")
library("xlsx")
library("jaffelab")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

## Set up plotting
plot_dir <- here("plots", "07_spatial_registration")
data_dir <- here("processed-data", "rdata", "spe", "07_spatial_registration")

## Load data 
# load(here("processed-data", "rdata","spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

## Load Registration Results
k_list <- c(7, 9, 16, 28)
names(k_list) <- paste0("k", k_list) ## Use paper naming convention

bayesSpace_registration_fn <- map(k_list, ~here(data_dir, paste0("dlpfc_pseudobulked_bayesSpace_specific_Ts_k",.x,".Rdata")))
bayesSpace_registration <- lapply(bayesSpace_registration_fn, function(x) get(load(x)))

## Select t-stats from the registration enrichment data

registration_t_stats <- map2(bayesSpace_registration, k_list, function(data, k){
 t_stats <- sapply(data, function(x) {
    x$t[, 2, drop = FALSE]
  })
 
 rownames(t_stats) <- rownames(data[[1]]$t)
 colnames(t_stats) <- paste0("Sp", k, "D", str_pad(colnames(t_stats),nchar(k),pad = "0"))
 
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

save(cor_top100, file = here(data_dir, "bayesSpacce_layer_cor_top100.Rdata"))

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
layer_anno_easy <- map2(cor_top100, names(cor_top100), function(cor, name){
  anno <- annotate_registered_clusters(cor_stats_layer = cor,
                                       confidence_threshold = 0.25,
                                       cutoff_merge_ratio = 0.25)
  return(anno)
})

layer_anno_strict <- map2(cor_top100, names(cor_top100), function(cor, name){
  anno <- annotate_registered_clusters(cor_stats_layer = cor,
                                       confidence_threshold = 0.25,
                                       cutoff_merge_ratio = 0.1)
  return(anno)
})

layer_anno <- c(layer_anno_easy[c("k7", "k9")], layer_anno_strict[c("k16", "k28")])

#### Annotate Cell Types by Layer ####
anno_abby <- data.frame(cluster = paste0("Sp9D",1:9), layer_abby = c("Vas", "L1", "L2/3", "L5", "L3", "WM", "L6A", "L4", "WM"))

layer_anno_strict$k9 |> arrange(cluster) |> left_join(anno_abby)
# cluster layer_confidence layer_label layer_abby
# 1   Sp9D1             good          L1        Vas
# 2   Sp9D2             good          L1         L1
# 3   Sp9D3             good          L2       L2/3
# 4   Sp9D4             good          L5         L5
# 5   Sp9D5             good          L3         L3
# 6   Sp9D6             good          WM         WM
# 7   Sp9D7             good          L6        L6A
# 8   Sp9D8             good          L4         L4
# 9   Sp9D9             good          WM         WM

layer_anno$k9 |> arrange(cluster) |> left_join(anno_abby)
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
  mutate(layer_annotation = fix_layer_order2(layer_label),
         layer_combo = factor(paste(cluster, "~", layer_annotation)),
         layer_combo2 = paste(layer_annotation, cluster),
         bayesSpace = factor(gsub("Sp","k",ss(cluster, "D")), levels = c("k7", "k9", "k16", "k28")), .before = cluster)  |>
  mutate(layer_combo = fct_reorder(layer_combo, layer_combo2, .desc = FALSE)) |>
  arrange(layer_combo) |>
  arrange(bayesSpace)

levels(layer_anno_all$layer_combo)
rownames(layer_anno_all) <- layer_anno_all$cluster

layer_anno_all |> count(layer_annotation)

## Save for reference
write.csv(layer_anno_all, file = here(data_dir, "bayesSpace_layer_annotations.csv"), row.names = FALSE)
save(layer_anno_all, file = here(data_dir, "bayesSpace_layer_annotations.Rdata")) ## save to preserve factors
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
  select(bayesSpace, layer_combo, cluster, layer_label) |>
  pivot_longer(!c(bayesSpace, layer_combo, cluster), names_to = "Annotation", values_to ="label") |>
  mutate(confidence = !grepl("\\*", label),
         layers = str_split(gsub("\\*","",label),"/"),
         Annotation = gsub("_label","", Annotation)) |>
  unnest_longer("layers") |>
  mutate(layer_short = ifelse(grepl("^[0-9]",layers), paste0("L",layers), layers),
         layer_long = gsub("L","Layer",layer_short)) |>
  select(-layers)

layer_anno_long |> count(layer_long, layer_short)
layer_anno_long |> count(confidence)

## Spot_plots 
bayes_layer_anno_plot <- layer_anno_long |>
  ggplot(aes(x = layer_short, y = layer_combo, color = layer_long)) +
  geom_point() +
  facet_grid(bayesSpace~., scales = "free_y", space = "free") +
  # facet_wrap(bayesSpace, scales = "free_y", ncol = 1) +
  scale_color_manual(values = libd_layer_colors) + 
  scale_y_discrete(limits=rev) ## WM on bottom

ggsave(bayes_layer_anno_plot, filename = here(plot_dir, "bayesSpace_layer_anno.png"))

#### bayesSpace Spatial Registration heatmaps ####

## match spatialLIBD color scale
theSeq <- seq(min(cor_kplus), max(cor_kplus), by = 0.01)
my.col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq))

# domain colors from plychrome
k_colors <- Polychrome::palette36.colors(28)
names(k_colors) <- c(1:28)

## Add intermediate colors to layers
source("libd_intermediate_layer_colors.R")
libd_intermediate_layer_colors <- c(spatialLIBD::libd_layer_colors, libd_intermediate_layer_colors)
names(libd_intermediate_layer_colors) <- gsub("ayer","",names(libd_intermediate_layer_colors))
libd_intermediate_layer_colors
# L1            L2            L3            L4            L5            L6            WM            NA 
# "#F0027F"     "#377EB8"     "#4DAF4A"     "#984EA3"     "#FFD700"     "#FF7F00"     "#1A1A1A" "transparent" 
# WM2          L1/2          L2/3          L3/4          L4/5          L5/6         L6/WM 
# "#666666"     "#BF3889"     "#50DDAC"     "#8278B0"     "#BD8339"     "#FFB300"     "#7A3D00" 

## build annotation matrix
anno_matrix <- layer_anno_long |> 
  mutate(fill = ifelse(confidence, "X","*")) |> 
  select(cluster, layer_long, fill) |> 
  pivot_wider(names_from = "layer_long", values_from = "fill", values_fill = "") |>
  column_to_rownames("cluster")

layer_anno_colors <- layer_anno_all |>
  mutate(domain_color = as.integer(gsub("Sp[0-9]+D","",cluster))) |>
  select(bayesSpace, cluster, layer_combo, domain_color, layer_annotation)

layer_color_bar <- columnAnnotation(" " = colnames(cor_top100$k7), 
                                    col = list(" " = spatialLIBD::libd_layer_colors),
                                    show_legend = FALSE)

## Just k7 plot 
layer_anno_k7 <- layer_anno_colors |> filter(bayesSpace == "k7")

cor_k7 <- cor_top100$k7[layer_anno_k7$cluster,]
rownames(cor_k7) <- layer_anno_k7$layer_combo

anno_matrix_k7 <- anno_matrix[grepl("Sp7", rownames(anno_matrix)), colnames(cor_k7)]
anno_matrix_k7 <- anno_matrix_k7[layer_anno_k7$cluster,]
rownames(anno_matrix_k7) <- layer_anno_k7$layer_combo

k7_color_bar <- rowAnnotation(df = layer_anno_k7 |>
                               select(domain_color, layer_anno = layer_annotation),
                             col = list(domain_color = k_colors,
                                        layer_anno = libd_intermediate_layer_colors),
                             show_legend = FALSE)

pdf(here(plot_dir,"bayesSpace_k7_spatial_registration_heatmap_color.pdf"), height = 4, width = 5.5)
Heatmap(cor_k7,
        name = "Cor",
        col = my.col,
        # row_split = layer_anno_all$bayesSpace,
        rect_gp = gpar(col = "black", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        right_annotation = k7_color_bar,
        bottom_annotation = layer_color_bar,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(anno_matrix_k7[i,j], x, y, gp = gpar(fontsize = 10))
        }
)
dev.off()

## 'kplus' = k9, 16, 28 in main plot
layer_anno_colors <- layer_anno_colors |> filter(bayesSpace != "k7")

cor_kplus <- do.call("rbind", cor_top100[c("k9", "k16", "k28")])
## order by bayesSpace annos
cor_kplus <- cor_kplus[layer_anno_colors$cluster,]
rownames(cor_kplus) <- layer_anno_colors$layer_combo

anno_matrix_kplus <- anno_matrix[!grepl("Sp7", rownames(anno_matrix)), colnames(cor_kplus)]

anno_matrix_kplus <- anno_matrix_kplus[layer_anno_colors$cluster,]
rownames(anno_matrix_kplus) <- layer_anno_colors$layer_combo


kplus_color_bar <- rowAnnotation(df = layer_anno_colors |>
                                   select(domain_color, layer_anno = layer_annotation),
                                 col = list(domain_color = k_colors,
                                            layer_anno = libd_intermediate_layer_colors),
                                 show_legend = FALSE)


pdf(here(plot_dir,"bayesSpace_kplus_spatial_registration_heatmap_color.pdf"), height = 10)
Heatmap(cor_kplus,
        name = "Cor",
        col = my.col,
        row_split = layer_anno_colors$bayesSpace,
        rect_gp = gpar(col = "black", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        right_annotation = kplus_color_bar,
        bottom_annotation = layer_color_bar,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(anno_matrix_kplus[i,j], x, y, gp = gpar(fontsize = 10))
        }
)
dev.off()

# sgejobs::job_single('02_cellType_correlation_annotation', create_shell = TRUE, memory = '5G', command = "Rscript 02_cellType_correlation_annotation.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
