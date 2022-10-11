
library("SpatialExperiment")
library("spatialLIBD")
library("tidyverse")
library("xlsx")
library("here")
library("sessioninfo")

## Set up plotting
plot_dir <- here("plots", "12_spatial_registration_sn")

## Load Registration Results
load(here("processed-data", "rdata", "spe", "12_spatial_registration_sn", "sn_hc_registration.RDS"), verbose = TRUE)

## Select t-stats from the registration enrichment data
registration_t_stats <- sn_hc_registration$enrichment[, grep("^t_stat", colnames(sn_hc_registration$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

#### Calculate Correlation Matrix ####
## get layer data
#### Load Layer and k Data  ####
layer_modeling_results <- fetch_data(type = "modeling_results")

paths <- list(k9 = "parsed_modeling_results_k9.Rdata", k16 = "parsed_modeling_results_k16.Rdata")

modeling_results <- lapply(paths, function(x) 
  get(load(here("processed-data","rdata","spe","08_layer_differential_expression",x))))

modeling_results <- c(list(layer = layer_modeling_results), modeling_results)
names(modeling_results)

#### Correlate with modeling results ####
cor_top100 <- map(modeling_results, ~layer_stat_cor(registration_t_stats,
                                                    .x,
                                                    model_type = "enrichment",
                                                    reverse = FALSE,
                                                    top_n = 100))
  
## Plot
pdf(here(plot_dir, paste0("spatial_registration_plot_sn.pdf")))
map(cor_top100, layer_stat_cor_plot)
dev.off()

#### Annotate Layers ####
layer_anno <- map2(cor_top100, names(cor_top100), function(cor, name){
  anno <- annotate_registered_clusters(cor_stats_layer = cor,
                                       confidence_threshold = 0.25,
                                       cutoff_merge_ratio = 0.25)
  colnames(anno) <- gsub("layer",name, colnames(anno))
  return(anno)
})

layer_anno_all <- reduce(layer_anno, left_join, by = "cluster")


modeling_results_pilot <- fetch_data(type = "modeling_results_pilot")

## Save correlation matrix
save(cor_top100, file = here("processed-data", "rdata", "spe", "12_spatial_registration_sn", "sn_hc_cor_top100.RDS"))

##  Plot layer correlation


pdf(here(plot_dir, "spatial_registration_plot_sn_v_manual_top100.pdf"))
layer_stat_cor_plot(cor_top100, max = max(cor_top100))
dev.off()


#### Annotate Cell Types by Layer ####

layer_anno <- annotate_registered_clusters(
    cor_stats_layer = cor_top100,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
)

layer_anno$layer |> arrange(cluster)
#         cluster layer_confidence layer_label
# 1         Astro             good          L1
# 2  EndoMural_01             good       L1/WM
# 3  EndoMural_02             good          L1
# 4      Excit_01             good          L3
# 5      Excit_02             good        L5/6
# 6      Excit_03             good          L4
# 7      Excit_04             good          L5
# 8      Excit_05             good          L3
# 9      Excit_06             good          L6
# 10     Excit_07             good          L5
# 11     Excit_08             good          L6
# 12     Excit_09             good      L4/3/5
# 13     Excit_10             good          L4
# 14     Excit_11             good      L4/5/3
# 15     Excit_12             poor       L4/5*
# 16     Excit_13             poor         L4*
# 17     Excit_14             good        L3/2
# 18     Excit_15             poor         L1*
# 19     Inhib_01             good          L2
# 20     Inhib_02             good          L4
# 21     Inhib_03             poor       L4/3*
# 22     Inhib_04             poor         L2*
# 23     Inhib_05             good          L2
# 24     Inhib_06             poor         L2*
# 25        Micro             good       WM/L1
# 26     Oligo_01             poor       L3/4*
# 27     Oligo_02             good          WM
# 28     Oligo_03             good          WM
# 29          OPC             good          WM

layer_anno |>
    filter(grepl("Excit", cluster), layer_confidence == "good") |>
    count(layer_label)
#   layer_label n
# 1          L3 2
# 2        L3/2 1
# 3          L4 2
# 4      L4/3/5 1
# 5      L4/5/3 1
# 6          L5 2
# 7        L5/6 1
# 8          L6 2

##  Mess with annotation plot
layer_order <- layer_anno |>
    filter(grepl("Excit", cluster), layer_confidence == "good") |>
    arrange(layer_label)

cor_temp <- cor_top100[layer_order$cluster, ]
cor_temp[cor_temp < 0.0] <- 0

pdf(here(plot_dir, "spatial_annotation_plot_sn_v_manual_top100.pdf"))
layer_stat_cor_plot(cor_temp, max = max(cor_temp))
dev.off()


#### Add Layer Annotations to colData ####
source(here("code","analysis", "12_spatial_registration_sn","utils.R"))

## layer_annotation is the reordered layer label - removes detail from the ordering process but helps group
layer_anno_all <- layer_anno_all |>
    arrange(cluster) |>
    mutate(
        layer_annotation = fix_layer_order2(layer_label),
        cellType_broad = gsub("_.*", "", cluster),
        cellType_layer = case_when(
            layer_confidence == "good" & grepl("Excit", cluster) ~ paste0(cellType_broad, "_", layer_annotation),
            grepl("Excit", cluster) ~ as.character(NA),
            TRUE ~ cellType_broad
        )
    ) |>
    select(-cellType_broad)

## Save for refrence
write.csv(layer_anno_all, file = here("processed-data", "rdata", "spe", "12_spatial_registration_sn", "cellType_layer_annotations.csv"))

## Add to sce object for future use
load(file = "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata", verbose = TRUE)

sce$cellType_layer <- factor(layer_anno2$cellType_layer[match(sce$cellType_hc, layer_anno2$cluster)],
    levels = c(
        "Astro", "EndoMural", "Micro", "Oligo", "OPC",
        "Excit_L2/3", "Excit_L3", "Excit_L3/4/5", "Excit_L4", "Excit_L5",
        "Excit_L5/6", "Excit_L6", "Inhib"
    )
)
sce$layer_annotation <- factor(layer_anno2$layer_annotation[match(sce$cellType_hc, layer_anno2$cluster)])

table(sce$cellType_layer)
# Astro    EndoMural        Micro        Oligo          OPC   Excit_L2/3     Excit_L3 Excit_L3/4/5     Excit_L4
# 3979         2157         1601        32051         1940           82        10459         3043         2388
# Excit_L5   Excit_L5/6     Excit_L6        Inhib
# 2505         2487         1792        11067

table(sce$layer_annotation)
# L1    L1*     L2    L2*   L2/3     L3  L3/4* L3/4/5     L4    L4*  L4/5*     L5   L5/6     L6     WM  WM/L1
# 5690     66   6558   1932     82  10459  24335   3043   3655   1567    420   2505   2487   1792  10966   2047

# save(sce, file = "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata")

#### Save Output to XLSX sheet ####
data_dir <- here("processed-data","rdata","spe","12_spatial_registration_sn")

key <- data.frame(data = c("annotation", paste0("cor_", names(modeling_results))),
                  description = c("Annotations of spatial registration, with coresponding layer cell type lables used in spatial deconvolution",
                                  "Correlation values with manual layer annotations",
                                  "Correlation values with k9 domains",
                                  "Correlation values with k16 domains"))

## Clear file and write key
annotation_xlsx <- here(data_dir,"sn_spatial_annotations.xlsx")
write.xlsx(key, file=annotation_xlsx, sheetName="Key", append=FALSE, row.names=FALSE)

## write annotations
write.xlsx(layer_anno_all, file=annotation_xlsx, sheetName= paste0("annotation"), append=TRUE, row.names=FALSE)

## write correlations 
walk2(cor_top100, names(cor_top100),
        ~write.xlsx(t(.x), file=annotation_xlsx, sheetName= paste0("cor_", .y), append=TRUE, row.names=TRUE))


#### Explore Annotations ####

layer_anno_long <- layer_anno_all |>
  select(cluster,ends_with("label")) |>
  pivot_longer(!cluster, names_to = "Annotation", values_to ="label") |>
  mutate(confidence = !grepl("\\*", label),
         layers = str_split(gsub("\\*","",label),"/"),
         Annotation = gsub("_label","", Annotation)) |>
  unnest_longer("layers") |>
  # mutate(layers = ifelse(Annotation == "layer" & grepl("^[0-9]",layers), paste0("L",layers), layers))
  mutate(layers = ifelse(Annotation == "layer",
                         ifelse(grepl("^[0-9]",layers), paste0("L",layers), layers),
                         paste0("k",str_pad(layers, 2,pad= "0"))
  ))

layer_anno_long |> count(confidence)


label_anno_plot <- layer_anno_long |>
  ggplot(aes(x = cluster, y = layers)) +
  geom_point(aes(color = confidence), size = 3) +
  # geom_tile(aes( fill = confidence), color = "black") +
  facet_wrap(~Annotation, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# ggsave(label_anno_plot, filename = here(plot_dir, "spatial_annotations_sn_all_tile.png"), height = 10)
ggsave(label_anno_plot, filename = here(plot_dir, "spatial_annotations_sn_all.png"), height = 10)

## which are specific?
n_anno <- layer_anno_long|> filter(confidence) |> group_by(Annotation, cluster) |> summarize(n_anno = n())

label_anno_plot_specific <- layer_anno_long |>
  left_join(n_anno) |>
  filter(confidence) |>
  ggplot(aes(x = cluster, y = layers)) +
  geom_point(aes(color = n_anno == 1), size = 3) +
  # geom_tile(aes( fill = confidence), color = "black") +
  facet_wrap(~Annotation, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(label_anno_plot_specific, filename = here(plot_dir, "spatial_annotations_sn_all_specific.png"), height = 10)
