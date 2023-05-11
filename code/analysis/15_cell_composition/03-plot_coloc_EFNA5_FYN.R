library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(here)
library(escheR)

spe_dat <- readRDS(here(
    "processed-data/rdata/spe/01_build_spe",
    "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
))

# Config ------------------------------------------------------------------

gene_1 <- "FYN"
gene_2 <- "EFNA5"

source(here("code", "analysis", "15_cell_composition", "fun_coloc.R"))
coloc_spe <- calc_coloc(spe_dat, gene_1, gene_2, sample_id = "Br8667_mid")


# Spatial Colocalization Plot ---------------------------------------------

vis_coloc(coloc_spe, gene_1, gene_2,
    sample_id = "Br8667_mid",
    save.path = here("plots/15_cell_composition", "coloc")
)

vars_spd <- c(
    "BayesSpace_harmony_09",
    "BayesSpace_harmony_16"
) |>
    set_names()

# Factor the Spd Labels ----
bayes_layers <- here(
    "processed-data", "rdata", "spe",
    "08_spatial_registration",
    "bayesSpace_layer_annotations.Rdata"
) |>
    load() |>
    get() |>
    select(Annotation = bayesSpace, layer_long = cluster, layer_combo) |>
    filter(Annotation %in% c("k09", "k16"))
#
# vis_coloc_spd(coloc_spe, gene_1, gene_2, sample_id = "Br8667_mid",
#               save.path = here("plots/15_cell_composition", "coloc"))



# Boxplot SpD -------------------------------------------------------------

full_coloc_spe <- calc_coloc(spe_dat, gene_1, gene_2, sample_id = NULL)

fnl_dat <- colData(full_coloc_spe) |> data.frame()

spd_dat <- vars_spd |>
    map_dfr(.f = function(var) {
        fnl_dat |>
            select(position, subject, sample_id,
                spd = var, coloc
            )
    }, .id = "spd_method") |>
    # group_split(spd_method, position, .keep = TRUE) |>
    group_split(spd_method, spd, sample_id, .keep = TRUE) |>
    map_dfr(.f = function(dat) {
        # browser()

        meta_df <- dat |>
            select(
                spd_method, spd, sample_id
            ) |>
            unique() |>
            mutate(
                k = str_sub(spd_method, start = -2, end = -1) |>
                    as.integer(),
                layer_long = sprintf("Sp%02dD%02d", k, spd)
            ) |>
            left_join(
                bayes_layers,
                by = "layer_long"
            )

        stopifnot(nrow(meta_df) == 1)

        data.frame(
            meta_df,
            n_spots = nrow(dat),
            prop.table(table(dat$coloc))
        )
    })



# * Plot K=9 without Neither ------------------------------------------------
layer_df <- bayes_layers |>
    filter(Annotation == "k09")

k <- 9

plot_dat <- spd_dat |>
    filter(
        Annotation == "k09",
        Var1 != "Neither"
    ) |>
    mutate(
        layer_combo = fct_drop(layer_combo),
        Var1 = factor(Var1,
            levels = c("co-localize", gene_1, gene_2),
            labels = c(
                paste0(gene_1, " & ", gene_2),
                paste0(gene_1, " only"),
                paste0(gene_2, " only")
            )
        )
    )

# Tests -------------------------------------------------------------------

test_dat <- plot_dat |>
    filter(Var1 == paste0(gene_1, " & ", gene_2)) |>
    transmute(
        L6_group = factor(layer_long == "Sp09D07"),
        Freq
    )

wilcox.test(Freq ~ L6_group, test_dat) # p=7.199e-05
t.test(Freq ~ L6_group, test_dat) # p=0.0045


# Box Plot ----------------------------------------------------------------

limit_fun <- function(y) {
    stopifnot(max(y) < 1)
    stopifnot(max(y) > 0)

    if ((max(y) * 100) %% 10 >= 5) {
        ul <- ceiling(max(y) * 10) / 10
    } else {
        ul <- floor(max(y) * 10) / 10 + 0.05
    }

    return(c(0, ul))
}

break_fun <- function(y) {
    ul <- limit_fun(y)[2]

    if (ul > 0.1) {
        return(seq(0, ul, 0.1))
    } else {
        return(c(0.00, 0.05, 0.10))
    }
}



ret_plot <- plot_dat |>
    mutate(
        Var1 = factor(Var1,
            labels = c(
                paste0(gene_1, " & ", gene_2, " (p=0.0046)"),
                paste0(gene_1, " only"),
                paste0(gene_2, " only")
            )
        )
    ) |>
    ggplot(
        aes(x = layer_combo, y = Freq, fill = Var1)
    ) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~Var1, scale = "free_y", ncol = 1) +
    scale_fill_manual(
        values = c("#d7191c", "#abd9e9", "#2c7bb6")
        # set_names(
        # Polychrome::palette36.colors(k)[str_sub(layer_df$layer_long, -2, -1) |>
        #                                     as.integer()],
        # layer_df$layer_combo
        # )
    ) +
    guides(x = guide_axis(angle = 90)) +
    scale_y_continuous(
        breaks = break_fun,
        limits = limit_fun,
        minor_breaks = NULL
    ) +
    ylab("Proportion of Spots") +
    xlab("") +
    theme_set(theme_bw(base_size = 20)) +
    theme(panel.grid.minor.y = element_blank())


ggsave(
    filename = here(
        "plots/15_cell_composition",
        # "layer_comp",
        paste0("coloc_comp_", gene_1, "-", gene_2, "_spd_9.pdf")
    ),
    plot = ret_plot,
    height = 7,
    width = 8
)
