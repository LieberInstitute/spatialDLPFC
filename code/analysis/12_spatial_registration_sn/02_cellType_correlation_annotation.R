
library("SpatialExperiment")
library("spatialLIBD")
library("dplyr")
library("here")
library("sessioninfo")

## Load Registration Results
load(here("processed-data", "rdata", "spe", "12_spatial_registration_sn", "sn_hc_registration.RDS"), verbose = TRUE)

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
save(cor_top100, file = here("processed-data", "rdata", "spe", "12_spatial_registration_sn", "sn_hc_cor_top100.RDS"))

##  Plot layer correlation
plot_dir <- here("plots", "12_spatial_registration_sn")

pdf(here(plot_dir, "spatial_registration_plot_sn_v_manual_top100.pdf"))
layer_stat_cor_plot(cor_top100, max = max(cor_top100))
dev.off()


#### Annotate Cell Types by Layer ####

layer_anno <- annotate_registered_clusters(
    cor_stats_layer = cor_top100,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
)

layer_anno[order(layer_anno$cluster), ]
#         cluster layer_confidence layer_label
# 3         Astro             good          L1
# 4  EndoMural_01             good       L1/WM
# 5  EndoMural_02             good          L1
# 8      Excit_01             good          L3
# 17     Excit_02             good        L5/6
# 13     Excit_03             good          L4
# 16     Excit_04             good          L5
# 9      Excit_05             good          L3
# 19     Excit_06             good          L6
# 18     Excit_07             good          L5
# 20     Excit_08             good          L6
# 14     Excit_09             good      L4/3/5
# 15     Excit_10             good          L4
# 11     Excit_11             good      L4/5/3
# 12     Excit_12             poor       L4/5*
# 22     Excit_13             poor         L4*
# 27     Excit_14             good        L3/2
# 21     Excit_15             poor         L1*
# 24     Inhib_01             good          L2
# 10     Inhib_02             good          L4
# 28     Inhib_03             poor       L4/3*
# 25     Inhib_04             poor         L2*
# 26     Inhib_05             good          L2
# 29     Inhib_06             poor         L2*
# 6         Micro             good       WM/L1
# 23     Oligo_01             poor       L3/4*
# 1      Oligo_02             good          WM
# 2      Oligo_03             good          WM
# 7           OPC             good          WM

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

## compare to old results
old_anno <- read.csv("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/05_explore_sce/DLPFC_HC_layer_annotation.csv")

old_anno <- old_anno |> select(cluster = cellType_hc, layer_label_old = layer_label)

layer_compare <- layer_anno |> left_join(old_anno)

layer_compare |> filter(grepl("Excit", cluster), layer_label == layer_label_old)
#    cluster layer_confidence layer_label layer_label_old
# 1 Excit_05             good          L3              L3
# 2 Excit_04             good          L5              L5
# 3 Excit_02             good        L5/6            L5/6
# 4 Excit_07             good          L5              L5
# 5 Excit_06             good          L6              L6
# 6 Excit_08             good          L6              L6

layer_compare |> filter(grepl("Excit", cluster), layer_label != layer_label_old)
#    cluster layer_confidence layer_label layer_label_old
# 1 Excit_01             good          L3            L3/2
# 2 Excit_11             good      L4/5/3              L4
# 3 Excit_12             poor       L4/5*             L4*
# 4 Excit_03             good          L4            L4/5
# 5 Excit_09             good      L4/3/5            L5/4
# 6 Excit_10             good          L4            L4/5
# 7 Excit_15             poor         L1*             L4*
# 8 Excit_13             poor         L4*         L4/3/5*
# 9 Excit_14             good        L3/2           L3/2*

layer_compare |> filter(grepl("Excit", cluster), layer_label != layer_label_old, layer_confidence == "good")
#    cluster layer_confidence layer_label layer_label_old
# 1 Excit_01             good          L3            L3/2
# 2 Excit_11             good      L4/5/3              L4
# 3 Excit_03             good          L4            L4/5
# 4 Excit_09             good      L4/3/5            L5/4
# 5 Excit_10             good          L4            L4/5
# 6 Excit_14             good        L3/2           L3/2*

layer_compare |>
    filter(grepl("Excit", cluster), layer_confidence == "good") |>
    arrange(layer_label) |>
    mutate(change = layer_label != layer_label_old)


##  Mess with annotation plot
layer_order <- layer_anno |>
    filter(grepl("Excit", cluster), layer_confidence == "good") |>
    arrange(layer_label)

cor_temp <- cor_top100[layer_order$cluster, ]
cor_temp[cor_temp < 0.0] <- 0

pdf(here(plot_dir, "spatial_annotation_plot_sn_v_manual_top100.pdf"))
layer_stat_cor_plot(cor_temp, max = max(cor_temp))
dev.off()

## make annotation grid
anno_matrix <- layer_anno |>
    mutate(layers = strsplit(layer_label, "/")) |>
    unnest_longer(layers) |>
    mutate(
        layers = gsub("\\*", "", gsub("^([1-9])", "L\\1", layers)),
        # layer_confidence = layer_confidence == "good",
        layer_confidence = ifelse(layer_confidence == "good", 1, .5)
    ) |>
    select(-layer_label) |>
    pivot_wider(names_from = "layers", values_from = "layer_confidence") |>
    column_to_rownames("cluster")

anno_matrix <- anno_matrix[rownames(cor_top100), gsub("ayer", "", colnames(cor_top100))]
colnames(anno_matrix) <- colnames(cor_top100)


anno_matrix2 <- cor_top100
anno_matrix2[is.na(anno_matrix)] <- NA

pdf(here(plot_dir, "spatial_annotation_plot_sn_v_manual_top100.pdf"))
# layer_matrix_plot(anno_matrix)
layer_stat_cor_plot(anno_matrix2)
dev.off()


#### Add Layer Annotations to colData ####

## Polish annotations
fix_layer_order <- function(l) {
    star <- ifelse(grepl("\\*", l), "*", "")
    l <- gsub("L|\\*", "", l)
    l <- sort(unlist(strsplit(l, "/")))

    if (all(l == "WM")) {
        return(paste0("WM", star))
    }

    l[[1]] <- paste0("L", l[[1]])
    if ("WM" %in% l) l <- c("WM", l[l != "WM"])
    fix <- paste0(paste0(l, collapse = "/"), star)

    return(fix)
}

fix_layer_order2 <- Vectorize(fix_layer_order)

## layer_annotation is the reordered layer label - removes detail from the ordering process but helps group
layer_anno2 <- layer_anno |>
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
write.csv(layer_anno2, file = here("processed-data", "rdata", "spe", "12_spatial_registration_sn", "cellType_layer_annotations.csv"))

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
