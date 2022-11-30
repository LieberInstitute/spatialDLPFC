library("spatialLIBD")
library("ggplot2")
library("ggpubr")
library("here")
library("purrr")
library("tidyr")
library("sessioninfo")

load(
    here(
        "processed-data",
        "rdata",
        "spe",
        "07_layer_differential_expression",
        "modeling_results_BayesSpace_k07.Rdata"
    ),
    verbose = TRUE
)

sce_pseudo <- readRDS(
    here(
        "processed-data",
        "rdata",
        "spe",
        "07_layer_differential_expression",
        "sce_pseudo_BayesSpace_k07.rds"
    )
)


allstats <-
    list("k07" = modeling_results, pilot = fetch_data("modeling_results"))
pseudo <-
    list("k07" = sce_pseudo, pilot = fetch_data(type = "sce_layer"))

sig_genes <-
    map2(
        allstats,
        pseudo,
        ~ sig_genes_extract(
            n = nrow(.y),
            modeling_results = .x,
            model_type = "enrichment",
            sce_layer = .y
        )
    )
top_pilot_genes <-
    unique(subset(sig_genes$pilot, top <= 100)$ensembl)
length(top_pilot_genes)
# [1] 692

ensembl_both <- intersect(sig_genes$k07$ensembl, top_pilot_genes)
length(ensembl_both)
# [1] 584

stats_long_WM_Sp07D07 <-
    map2_dfr(sig_genes, names(sig_genes), function(x, y) {
        res <- subset(x, ensembl %in% ensembl_both)
        if (y == "k07") {
            res <- subset(res, test == "Sp07D07")
        } else {
            res <- subset(res, test == "WM")
        }
        res$set <- y
        return(res)
    })
dim(stats_long_WM_Sp07D07)
# [1] 1168   10
stopifnot(nrow(stats_long_WM_Sp07D07) / 2 == length(ensembl_both))

stats_wide_WM_Sp07D07 <-
    tidyr::pivot_wider(
        stats_long_WM_Sp07D07,
        id_cols = ensembl,
        names_from = set,
        values_from = stat
    ) |>
    left_join(subset(sig_genes$pilot, top <= 100)[, c("ensembl", "test")]) |>
    mutate(
        layer = test,
        layer_color = spatialLIBD::libd_layer_colors[layer]
    ) |>
    select(-c(test))

head(stats_wide_WM_Sp07D07)

pdf(
    file = here::here(
        "plots",
        "08_spatial_registration",
        "ggplot_t_cor_k7_wm_colored.pdf"
    )
)
ggplot(
    stats_wide_WM_Sp07D07,
    aes(x = k07, y = pilot)
) +
    geom_point(aes(color = layer), alpha = 0.8, size = 3) +
    scale_color_manual(
        values = spatialLIBD::libd_layer_colors,
        name = "Marker for:",
        guide = guide_legend(override.aes = list(size = 6))
    ) +
    labs(
        x = "t-statistic Sp07D07 > rest\nunsurpervised clustered new data\nobserved 584/692 layer marker genes",
        y = "t-statistic WM > rest\nmanually annotated prior data\ntop 100 layer marker genes; 692 unique",
        title = "spatial registration: enrichment stats"
    ) +
    theme_bw(base_size = 20) +
    geom_smooth(
        method = lm,
        se = FALSE,
        colour = "grey30"
    ) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.24),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
    ) +
    # guides(colour = guide_legend(override.aes = list(size = 4))) +
    stat_cor(
        mapping = aes(x = k07, y = pilot),
        method = "pearson",
        label.x = -10,
        label.y = 16,
        size = 7
    )
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
