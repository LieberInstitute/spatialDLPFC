
library("SpatialExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")

## Load Registration Results
load(here("processed-data","rdata","spe","12_spatial_registration_sn","sn_hc_registration.RDS"), verbose = TRUE)

## Select t-stats from the registration enrichment data
registration_t_stats <- sn_hc_registration$enrichment[, grep("^t_stat", colnames(sn_hc_registration$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))

#### Calculate Correlation Matrix ####
## get layer data
modeling_results <- fetch_data(type = "modeling_results")

## Correlate t-stats vs. manual layer annotations for top 100 genes
cor_top100 <- layer_stat_cor(
  registration_t_stats,
  modeling_results,
  model_type = "enrichment",
  reverse = FALSE,
  top_n = 100
)

## Save correlation matrix
save(cor_top100, file = here("processed-data","rdata","spe","12_spatial_registration_sn","sn_hc_cor_top100.RDS"))

##  Plot layer correlation
plot_dir <- here("plots","12_spatial_registration_sn")

pdf(here(plot_dir, "spatial_registration_plot_sn_v_manual_top100.pdf"))
layer_stat_cor_plot(cor_top100, max = max(cor_top100))
dev.off()


#### Annotate Cell Types by Layer ####

layer_anno <- annotate_registered_clusters(cor_stats_layer = cor_top100,
                                           confidence_threshold = 0.25,
                                           cutoff_merge_ratio = 0.25)

layer_anno[order(layer_anno$cluster),]
#         cluster layer_confidence layer_label
# 3         Astro             good          L1
# 4  EndoMural_01             good          L1
# 5  EndoMural_02             good          L1
# 19     Excit_01             good          L3
# 17     Excit_02             good        L5/6
# 13     Excit_03             good          L4
# 16     Excit_04             good          L5
# 20     Excit_05             good          L3
# 8      Excit_06             good          L6
# 18     Excit_07             good          L5
# 9      Excit_08             good          L6
# 14     Excit_09             good      L4/3/5
# 15     Excit_10             good          L4
# 11     Excit_11             good        L4/5
# 12     Excit_12             poor       L4/5*
# 22     Excit_13             poor       L3/1*
# 26     Excit_14             poor       L3/2*
# 21     Excit_15             poor         L1*
# 27     Inhib_01             good        L2/3
# 10     Inhib_02             good          L4
# 28     Inhib_03             poor         L4*
# 24     Inhib_04             poor         L2*
# 25     Inhib_05             good          L2
# 29     Inhib_06             poor   L2/4/3/5*
# 6         Micro             good       WM/L1
# 23     Oligo_01             poor         L3*
# 1      Oligo_02             good          WM
# 2      Oligo_03             good          WM
# 7           OPC             good          WM


## compare to old results
library(dplyr)
old_anno <- read.csv("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/05_explore_sce/DLPFC_HC_layer_annotation.csv")

old_anno <- old_anno |> select(cluster = cellType_hc, layer_label_old = layer_label)

layer_comapre <- layer_anno |> left_join(old_anno)

layer_comapre |> filter(grepl("Excit", cluster), layer_label == layer_label_old)
# cluster layer_confidence layer_label layer_label_old
# 1 Excit_06             good          L6              L6
# 2 Excit_08             good          L6              L6
# 3 Excit_04             good          L5              L5
# 4 Excit_07             good          L5              L5
# 5 Excit_05             good          L3              L3
# 6 Excit_14             poor       L3/2*           L3/2*

layer_comapre |> filter(grepl("Excit", cluster), layer_label != layer_label_old, layer_confidence == "good")
