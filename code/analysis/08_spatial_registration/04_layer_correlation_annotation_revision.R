library("spatialLIBD")
library("tidyverse")
library("xlsx")
library("jaffelab")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

## Input dir
dir_input <- here::here(
    "processed-data",
    "rdata",
    "spe",
    "07_layer_differential_expression"
)

## Set up plotting
plot_dir <- here("plots", "08_spatial_registration","cor_all_gene")
data_dir <-
    here("processed-data", "rdata", "spe", "08_spatial_registration")

## Load data
# load(here("processed-data", "rdata","spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

## Load Registration Results
k_list <- c(2, 7, 9, 16, 28)
names(k_list) <-
    paste0("k", sprintf("%02d", k_list)) ## Use paper naming convention

bayesSpace_registration_fn <-
    map(k_list, ~ here(
        dir_input,
        paste0(
            "modeling_results_BayesSpace_k",
            sprintf("%02d", .x),
            ".Rdata"
        )
    ))
bayesSpace_registration <-
    lapply(bayesSpace_registration_fn, function(x) {
        get(load(x))
    })

## Select t-stats from the registration enrichment data

registration_t_stats <-
    map(bayesSpace_registration, function(data) {
        x <- data$enrichment
        t_stats <- x[, grep("^t_stat_", colnames(x))]
        colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))
        return(t_stats)
    })

map(registration_t_stats, jaffelab::corner)

#### Calculate Correlation Matrix ####
## get layer data
layer_modeling_results <- fetch_data(type = "modeling_results")

#### Correlate with modeling results ####
## layer stat cor with ALL genes
cor_ALL <- map(
    registration_t_stats,
    ~ layer_stat_cor(
        .x,
        layer_modeling_results,
        model_type = "enrichment",
        reverse = FALSE
    )
)

save(cor_ALL,
    file = here(data_dir, "bayesSpacce_layer_cor_ALL.Rdata")
)

## Plot all for portability
# pdf(here(plot_dir, "cor_ALL_spatial_registration.pdf"))
# map(cor_ALL, layer_stat_cor_plot, max = 1)
# dev.off()

## Plot separately for illustrator
# map2(cor_ALL, names(cor_ALL), function(data, name){
#
#   pdf(here(plot_dir, paste0("spatial_registration_plot_sn-",name,".pdf")))
#   layer_stat_cor_plot(data)
#   dev.off()
#
# })

#### Annotate Layers ####
layer_anno_easy <-
    map2(cor_ALL, names(cor_ALL), function(cor, name) {
        anno <- annotate_registered_clusters(
            cor_stats_layer = cor,
            confidence_threshold = 0.25,
            cutoff_merge_ratio = 0.25
        )
        return(anno)
    })

layer_anno_strict <-
    map2(cor_ALL, names(cor_ALL), function(cor, name) {
        anno <- annotate_registered_clusters(
            cor_stats_layer = cor,
            confidence_threshold = 0.25,
            cutoff_merge_ratio = 0.1
        )
        return(anno)
    })

## We'll use the strict annotations
layer_anno <- layer_anno_strict

#### Annotate Cell Types by Layer ####
anno_abby <-
    data.frame(
        cluster = paste0("Sp09D0", 1:9),
        layer_abby = c("Vas", "L1", "L2/3", "L5", "L3", "WM", "L6A", "L4", "WM")
    )

layer_anno_strict$k09 |>
    arrange(cluster) |>
    left_join(anno_abby)
# cluster layer_confidence layer_label layer_abby
# 1 Sp09D01             good          L1        Vas
# 2 Sp09D02             good          L1         L1
# 3 Sp09D03             good          L2       L2/3
# 4 Sp09D04             good          L5         L5
# 5 Sp09D05             good          L3         L3
# 6 Sp09D06             good          WM         WM
# 7 Sp09D07             good          L6        L6A
# 8 Sp09D08             good          L4         L4
# 9 Sp09D09             good          WM         WM

layer_anno_easy$k09 |>
    arrange(cluster) |>
    left_join(anno_abby)
#   cluster layer_confidence layer_label layer_abby
# 1 Sp09D01             good          L1        Vas
# 2 Sp09D02             good          L1         L1
# 3 Sp09D03             good        L2/3       L2/3
# 4 Sp09D04             good          L5         L5
# 5 Sp09D05             good          L3         L3
# 6 Sp09D06             good          WM         WM
# 7 Sp09D07             good          L6        L6A
# 8 Sp09D08             good          L4         L4
# 9 Sp09D09             good          WM         WM

## Add additonal annotaitons
source(
    here("code", "analysis", "12_spatial_registration_sn", "utils.R"),
    echo = TRUE,
    max.deparse.length = 500
)

## layer_annotation is the reordered layer label - removes detail from the ordering process but helps group
layer_anno_all <- do.call("rbind", layer_anno) |>
    mutate(
        layer_annotation = fix_layer_order2(layer_label),
        layer_combo = factor(paste(cluster, "~", layer_annotation)),
        layer_combo2 = paste(layer_annotation, cluster),
        bayesSpace = factor(gsub("Sp", "k", ss(cluster, "D")), levels = names(layer_anno)),
        .before = cluster
    ) |>
    mutate(layer_combo = fct_reorder(layer_combo, layer_combo2, .desc = FALSE)) |>
    arrange(layer_combo) |>
    arrange(bayesSpace)

levels(layer_anno_all$layer_combo)
rownames(layer_anno_all) <- layer_anno_all$cluster

layer_anno_all |> count(layer_annotation)
# layer_annotation  n
# 1                L1 12
# 2                L2  3
# 3              L2/3  1
# 4                L3  8
# 5              L3/4  4
# 6                L4  4
# 7              L4/5  1
# 8                L5  7
# 9              L5/6  2
# 10               L6  4
# 11               WM 14

## Save for reference
write.csv(
    layer_anno_all,
    file = here(data_dir, "bayesSpace_layer_annotations_ALLgene.csv"),
    row.names = FALSE
)
save(layer_anno_all,
    file = here(data_dir, "bayesSpace_layer_annotations_ALLgene.Rdata")
) ## save to preserve factors
# layer_anno_all <- read_csv(here(data_dir, "cellType_layer_annotations.csv"))


#### Save Output to XLSX sheet ####
# key <- data.frame(
#     data = c("annotation", paste0("cor_", names(
#         registration_t_stats
#     ))),
#     description = c(
#         "Annotations of baySpace Domains",
#         paste0(
#             "Correlation values vs. manual annotation for ",
#             names(registration_t_stats)
#         )
#     )
# )
# 
# ## Clear file and write key
# annotation_xlsx <-
#     here(data_dir, "bayesSpace_layer_cor_annotations.xlsx")
# write.xlsx(
#     key,
#     file = annotation_xlsx,
#     sheetName = "Key",
#     append = FALSE,
#     row.names = FALSE
# )
# 
# ## write annotations
# write.xlsx(
#     layer_anno_all,
#     file = annotation_xlsx,
#     sheetName = paste0("annotation"),
#     append = TRUE,
#     row.names = FALSE
# )
# 
# ## write correlations
# walk2(
#     cor_ALL,
#     names(cor_ALL),
#     ~ write.xlsx(
#         t(.x),
#         file = annotation_xlsx,
#         sheetName = paste0("cor_", .y),
#         append = TRUE,
#         row.names = TRUE
#     )
# )

#### Explore Annotations ####
layer_anno_long <- layer_anno_all |>
    select(bayesSpace, layer_combo, cluster, layer_label) |>
    pivot_longer(!c(bayesSpace, layer_combo, cluster),
        names_to = "Annotation",
        values_to = "label"
    ) |>
    mutate(
        confidence = !grepl("\\*", label),
        layers = str_split(gsub("\\*", "", label), "/"),
        Annotation = gsub("_label", "", Annotation)
    ) |>
    unnest_longer("layers") |>
    mutate(
        layer_short = ifelse(grepl("^[0-9]", layers), paste0("L", layers), layers),
        layer_long = gsub("L", "Layer", layer_short)
    ) |>
    select(-layers)

layer_anno_long |> count(layer_long, layer_short)
# # A tibble: 7 × 3
#   layer_long layer_short     n
#   <chr>      <chr>       <int>
# 1 Layer1     L1             12
# 2 Layer2     L2              5
# 3 Layer3     L3             13
# 4 Layer4     L4              9
# 5 Layer5     L5              7
# 6 Layer6     L6              8
# 7 WM         WM             14
layer_anno_long |> count(confidence)
# # A tibble: 1 × 2
#   confidence     n
#   <lgl>      <int>
# 1 TRUE          68

# ## Spot_plots
# bayes_layer_anno_plot <- layer_anno_long |>
#     ggplot(aes(x = layer_short, y = layer_combo, color = layer_long)) +
#     geom_point() +
#     facet_grid(bayesSpace ~ ., scales = "free_y", space = "free") +
#     # facet_wrap(bayesSpace, scales = "free_y", ncol = 1) +
#     scale_color_manual(values = libd_layer_colors) +
#     scale_y_discrete(limits = rev) ## WM on bottom
# 
# ggsave(bayes_layer_anno_plot,
#     filename = here(plot_dir, "bayesSpace_layer_anno.png")
# )

## To switch the order in order to have L1 to L6, then WM on the x-axis
layer_order <- c(paste0("Layer", 1:6), "WM")
cor_ALL <- lapply(cor_ALL, function(x) {
    x[, layer_order]
})

#### bayesSpace Spatial Registration heatmaps ####
## color set up
## match spatialLIBD color scale
cor_kplus <- do.call("rbind", cor_ALL[c("k09", "k16", "k28")])
max(cor_kplus)
# [1] 0.9452202
theSeq <- seq(min(cor_kplus), max(cor_kplus), by = 0.01)
my.col <-
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq))

# domain colors from polychrome
k_colors <- Polychrome::palette36.colors(28)
names(k_colors) <- c(1:28)

## Add intermediate colors to layers
source(
    here(
        "code",
        "analysis",
        "08_spatial_registration",
        "libd_intermediate_layer_colors.R"
    ),
    echo = TRUE,
    max.deparse.length = 500
)
libd_intermediate_layer_colors <-
    c(
        spatialLIBD::libd_layer_colors,
        libd_intermediate_layer_colors
    )
names(libd_intermediate_layer_colors) <-
    gsub("ayer", "", names(libd_intermediate_layer_colors))
libd_intermediate_layer_colors
# L1            L2            L3            L4            L5            L6            WM            NA
# "#F0027F"     "#377EB8"     "#4DAF4A"     "#984EA3"     "#FFD700"     "#FF7F00"     "#1A1A1A" "transparent"
# WM2          L1/2          L2/3          L3/4          L4/5          L5/6         L6/WM
# "#666666"     "#BF3889"     "#50DDAC"     "#8278B0"     "#BD8339"     "#FFB300"     "#7A3D00"

## build annotation matrix
anno_matrix <- layer_anno_long |>
    mutate(fill = ifelse(confidence, "X", "*")) |>
    select(cluster, layer_long, fill) |>
    pivot_wider(
        names_from = "layer_long",
        values_from = "fill",
        values_fill = ""
    ) |>
    column_to_rownames("cluster")

layer_anno_colors <- layer_anno_all |>
    mutate(domain_color = as.integer(gsub("Sp[0-9]+D", "", cluster))) |>
    select(
        bayesSpace,
        cluster,
        layer_combo,
        domain_color,
        layer_annotation
    )

layer_color_bar <- columnAnnotation(
    " " = colnames(cor_ALL$k07),
    col = list(" " = spatialLIBD::libd_layer_colors),
    show_legend = FALSE
)

## Just one k at a time
registration_one_k <- function(k) {
    k_long <- paste0("k", sprintf("%02d", k))
    layer_anno_subset <-
        layer_anno_colors |> filter(bayesSpace == k_long)

    cor_subset <- cor_ALL[[k_long]][layer_anno_subset$cluster, ]
    rownames(cor_subset) <- layer_anno_subset$layer_combo

    anno_matrix_subset <-
        anno_matrix[grepl(paste0("Sp", sprintf("%02d", k)), rownames(anno_matrix)), colnames(cor_subset)]
    anno_matrix_subset <-
        anno_matrix_subset[layer_anno_subset$cluster, ]
    rownames(anno_matrix_subset) <- layer_anno_subset$layer_combo


    k_colors_current <- k_colors
    if (k == 2) {
        names(k_colors_current)[1:2] <- 2:1
    }

    subset_color_bar <- rowAnnotation(
        df = layer_anno_subset |>
            select(domain = domain_color, layer = layer_annotation),
        col = list(
            domain = k_colors_current,
            layer = libd_intermediate_layer_colors
        ),
        show_legend = FALSE
    )

    pdf(
        here(
            plot_dir,
            paste0(
                "bayesSpace_",
                k_long,
                "_spatial_registration_heatmap_color_ALLgene.pdf"
            )
        ),
        height = 4,
        width = 5.5
    )
    p <- Heatmap(
        cor_subset,
        name = "Cor",
        col = my.col,
        # row_split = layer_anno_all$bayesSpace,
        rect_gp = gpar(col = "black", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        right_annotation = subset_color_bar,
        bottom_annotation = layer_color_bar,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(anno_matrix_subset[i, j], x, y, gp = gpar(fontsize = 10))
        }
    )
    print(p)
    dev.off()
}

walk(k_list, registration_one_k)

## 'kplus' = k09, 16, 28 in main plot
layer_anno_colors <-
    layer_anno_colors |> filter(bayesSpace %in% c("k09", "k16", "k28"))

## order by bayesSpace annos
cor_kplus <- cor_kplus[layer_anno_colors$cluster, ]
rownames(cor_kplus) <- layer_anno_colors$layer_combo

anno_matrix_kplus <-
    anno_matrix[!grepl("Sp02|Sp07", rownames(anno_matrix)), colnames(cor_kplus)]

anno_matrix_kplus <- anno_matrix_kplus[layer_anno_colors$cluster, ]
rownames(anno_matrix_kplus) <- layer_anno_colors$layer_combo


kplus_color_bar <- rowAnnotation(
    df = layer_anno_colors |>
        select(domain = domain_color, layer = layer_annotation),
    col = list(
        domain = k_colors,
        layer = libd_intermediate_layer_colors
    ),
    show_legend = FALSE
)


pdf(
    here(
        plot_dir,
        "bayesSpace_kplus_spatial_registration_heatmap_color_ALLgene.pdf"
    ),
    height = 10
)
Heatmap(
    cor_kplus,
    name = "Cor",
    col = my.col,
    row_split = layer_anno_colors$bayesSpace,
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = kplus_color_bar,
    bottom_annotation = layer_color_bar,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(anno_matrix_kplus[i, j], x, y, gp = gpar(fontsize = 10))
    }
)
dev.off()


hm_color_long <- Heatmap(
    cor_kplus,
    name = "Cor",
    col = my.col,
    row_split = layer_anno_colors$bayesSpace,
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = kplus_color_bar,
    bottom_annotation = layer_color_bar,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(anno_matrix_kplus[i, j], x, y, gp = gpar(fontsize = 10))
    }
)

pdf(
    here(
        plot_dir,
        "bayesSpace_kplus_spatial_registration_heatmap_color_long_ALLgene.pdf"
    ),
    height = 10, width = 5
    # height = 9.5, width = 5
)

draw(hm_color_long, show_heatmap_legend = FALSE)

dev.off()

## Will create seperate legend and place in Ai
# lgd2 = Legend(col_fun = circlize::colorRamp2(seq(-1, 1, by = 0.5),
#                                   RColorBrewer::brewer.pal(5, "PRGn")),
#               title = "cor",
#               direction = "horizontal")
#
# pdf(here(plot_dir, "cor_legend.pdf"), height = 1, width = 2)
# draw(lgd2)
# dev.off()

#### Compare correlations ####
# load(here(data_dir, "bayesSpacce_layer_cor_ALL.Rdata"), verbose = TRUE)
# load(here(data_dir, "bayesSpace_layer_annotations_ALLgene.Rdata"))

## load top100 data
load(here(data_dir, "bayesSpacce_layer_cor_top100.Rdata"), verbose = TRUE)
# load(here(data_dir, "bayesSpace_layer_annotations.Rdata"), verbose = TRUE)
layer_anno_top100 <- read.csv(here(data_dir, "bayesSpace_layer_annotations.csv"))

cor_ALL_long <- do.call("rbind", map2(cor_ALL, names(cor_ALL), ~.x |> 
                               as.data.frame() |>
                               rownames_to_column("cluster") |>
                               pivot_longer(!cluster, names_to = "layer_long", values_to = "cor_ALL") |>
                               mutate(k = .y)
)) |> left_join(layer_anno_long |> select(cluster, layer_long, anno_ALL = confidence)) |>
  replace_na(list(anno_ALL = FALSE))

## top 100
layer_anno_top100_long <- layer_anno_top100 |>
  select(bayesSpace, layer_combo, cluster, layer_label) |>
  pivot_longer(!c(bayesSpace, layer_combo, cluster),
               names_to = "Annotation",
               values_to = "label"
  ) |>
  mutate(
    confidence = !grepl("\\*", label),
    layers = str_split(gsub("\\*", "", label), "/"),
    Annotation = gsub("_label", "", Annotation)
  ) |>
  unnest_longer("layers") |>
  mutate(
    layer_short = ifelse(grepl("^[0-9]", layers), paste0("L", layers), layers),
    layer_long = gsub("L", "Layer", layer_short)
  ) |>
  select(cluster, layer_long, anno_top100 = confidence)

cor_top100_long <- do.call("rbind", map2(cor_top100, names(cor_top100), ~.x |> 
                                        as.data.frame() |>
                                        rownames_to_column("cluster") |>
                                        pivot_longer(!cluster, names_to = "layer_long", values_to = "cor_top100") |>
                                        mutate(k = .y)
)) |> left_join(layer_anno_top100_long)|>
  replace_na(list(anno_top100 = FALSE))


cor_compare <- cor_top100_long |> 
  left_join(cor_ALL_long) |>
  mutate(anno_match = case_when(anno_ALL & anno_top100 ~ "match",
                                anno_ALL ~ "only ALL",
                                anno_top100 ~ "only top100",
                                TRUE ~ "none"))

cor_compare |> count(anno_match)
cor_compare |> filter(anno_match %in% c("ALL","top100"))

cor_compare |> filter(anno_match != "none")|> count(k, anno_match)

cor_compare |> 
  mutate(diff = cor_ALL - cor_top100)  |>
  group_by(k) |>
  summarize(mean_d = mean(diff),
            max_d = max(abs(diff)))

## plot scatter plot

cor_compare_scatter <- cor_compare |>
  ggplot(aes(x = cor_top100, y = cor_ALL, color = anno_match)) +
  geom_point() +
  facet_wrap(~k) +
  geom_abline() +
  coord_equal() +
  theme_bw()

ggsave(cor_compare_scatter, filename = here(plot_dir, "cor_compare_scatter.png"), width = 10)

cor_diff_scatter <- cor_compare |> 
  mutate(diff = cor_ALL - cor_top100) |>
  ggplot(aes(x = cor_top100, y = diff, color = anno_match)) +
  geom_point() +
  facet_wrap(~k) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  labs(y = "cor_ALL - cor_top100")

ggsave(cor_diff_scatter , filename = here(plot_dir, "cor_diff_scatter.png"))


# sgejobs::job_single('04_layer_correlation_annotation_revision', create_shell = TRUE, memory = '5G', command = "Rscript 04_layer_correlation_annotation_revision.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
