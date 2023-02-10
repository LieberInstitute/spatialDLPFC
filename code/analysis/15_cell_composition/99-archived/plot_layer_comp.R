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

fnl_dat <- colData(spe_dat) |> data.frame()


vars_spd <- c(
    "BayesSpace_harmony_09",
    "BayesSpace_harmony_16"
) |>
    set_names()


# Factor the Spd Labels
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
layer_dat <- vars_spd |>
    map_dfr(.f = function(var) {
        fnl_dat |>
            select(position, subject,
                spd = var
            )
    }, .id = "spd_method") |>
    group_split(spd_method, position, .keep = TRUE) |>
    map_dfr(.f = function(dat) {
        # browser()
        dat |>
            group_by(subject, spd) |>
            summarize(n_spots = n()) |>
            mutate(
                tot_spots = sum(n_spots),
                layer_perc = n_spots / tot_spots,
                spd_method = unique(dat$spd_method),
                position = unique(dat$position)
            ) |>
            ungroup()
    })


# Generate a List of ggplots --------------------------------------------------------
layer_plots <- layer_dat |>
    group_split(spd_method, .keep = TRUE) |>
    map(.f = function(sum_dat) {
        k <- sum_dat$spd_method |>
            unique() |>
            str_sub(start = -2, end = -1) |>
            as.integer()

        layer_df <- bayes_layers |>
            filter(Annotation == sprintf("k%02d", k))

        ret_plot <- sum_dat |>
            mutate(
                spd = factor(spd,
                    levels = str_sub(layer_df$layer_long, -2, -1) |>
                        as.integer(),
                    labels = layer_df$layer_combo
                ),
                position = str_to_sentence(position)
            ) |>
            ggplot(aes(x = subject, y = layer_perc, fill = spd)) +
            geom_bar(position = "fill", stat = "identity") +
            scale_fill_manual(
                name = "BayesSpace Domain",
                values = set_names(
                    Polychrome::palette36.colors(k)[str_sub(layer_df$layer_long, -2, -1) |>
                        as.integer()],
                    layer_df$layer_combo
                )
            ) +
            guides(x = guide_axis(angle = -45)) +
            facet_wrap(~position) +
            ylab("Proportion of Spots") +
            xlab("") +
            theme_set(theme_bw(base_size = 20))
    })

# Create Plot Grid --------------------------------------------------------
layer_grid_plot <- ggarrange(
    plotlist = layer_plots,
    ncol = 2, nrow = 1
)

ggsave(
    filename = here(
        "plots/15_cell_composition",
        "layer_comp",
        paste0("layer_comp_grid.pdf")
    ),
    plot = layer_grid_plot,
    height = 7,
    width = 12
)
