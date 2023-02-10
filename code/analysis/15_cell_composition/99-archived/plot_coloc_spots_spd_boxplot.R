# Load Packages
library(tidyverse)
library(here)
library(SpatialExperiment)
library(ggpubr)

# Load Spe Object
spe_dat <- readRDS(
    here(
        "processed-data/rdata/spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
    )
)


coloc_spe <- calc_coloc(spe_dat, "EFNA5", "EPHA5", sample_id = NULL)

fnl_dat <- colData(coloc_spe) |> data.frame()


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


# Curate Data -------------------------------------------------------------
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


# Plots -------------------------------------------------------------------

ret_plot_list <- spd_dat |>
    # mutate(
    #     layer_combo = as.character(layer_combo)
    # ) |>
    group_split( Annotation, .keep = TRUE) |>
    map(.f = function(dat){
        # browser()
        dat |>
            mutate(
                layer_combo = fct_drop(layer_combo)
            ) |>
            ggplot(
                aes(x = layer_combo, y = Freq, fill = Var1)
            ) +
            geom_boxplot() +
            # scale_fill_manual(
            #     name = "BayesSpace Domain",
            #     values = set_names(
            #         Polychrome::palette36.colors(k)[str_sub(layer_df$layer_long, -2, -1) |>
            #                                             as.integer()],
            #         layer_df$layer_combo
            #     )
            # ) +
            guides(x = guide_axis(angle = 90)) +
            facet_wrap(~Var1,  scale = "free_y", ncol = 1)
        # ylab("Proportion of Spots") +
        # xlab("") +
        # theme_set(theme_bw(base_size = 20))
    })


ggsave(
    filename = here(
        "plots/15_cell_composition",
        # "layer_comp",
        paste0("coloc_comp_spd.pdf")
    ),
    plot = ggpubr::ggarrange(
        # nrow = 2,
        ncol = 2,
        plotlist = ret_plot_list,
        common.legend = TRUE),
    width = 10
)


# Plot K=9 without Neither ------------------------------------------------
ret_plot <- spd_dat |>
    filter(
        Annotation == "k09",
        Var1 != "Neither"
    ) |>
    mutate(
        layer_combo = fct_drop(layer_combo),
        Var1 = factor(
            Var1,
            levels = c("co-localize", "EFNA5", "EPHA5"),
            labels = c("Both", "EFNA5 only", "EPHA5 only"))
    ) |>
    ggplot(
        aes(x = layer_combo, y = Freq, fill = Var1)
    ) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(~Var1,  scale = "free_y", ncol = 1) +
    guides(x = guide_axis(angle = 90)) +
    ylab("Proportion of Spots") +
    xlab("") +
    theme_set(theme_bw(base_size = 20))


ggsave(
    filename = here(
        "plots/15_cell_composition",
        # "layer_comp",
        paste0("coloc_comp_spd_9.pdf")
    ),
    plot = ret_plot,
    height = 7,
    width = 8
)


