library(tidyverse)
library(here)
library(SpatialExperiment)

spe_dat <- readRDS(here("processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters_and_deconvolution_results.rds"))

fnl_dat <- colData(spe_dat) |> data.frame()

# sp9_colors <- Polychrome::palette36.colors(9)
# names(sp9_colors) <- paste0("Sp9D", 1:9)
#
# sp16_colors <- Polychrome::palette36.colors(16)
# names(sp16_colors) <- paste0("Sp16D", 1:16)

fnl_dat |>
    group_by(position) |>
    summarize(
        n_sample = n_distinct(subject),
        n_spots = n(),
        n_sp9_miss = sum(is.na(BayesSpace_harmony_09)),
        n_sp16_miss = sum(is.na(BayesSpace_harmony_16))
    )

## Hypothesis testing: Are the number of spots different across pos
## NOTE: no great power
# spot_lm <- fnl_dat |>
#     group_by(subject, position) |>
#     summarize(n_spots = n()) |>
#     with(lm(n_spots~position))
# summary(spot_lm)
# anova(spot_lm)


# Visualize the number of in-tissue spots per sample-------------------------------------------

stopifnot(all(fnl_dat$in_tissue))

ggsave(
    filename = here("plots/15_cell_composition/", "n_spots.pdf"),
    fnl_dat |>
        group_by(subject, position) |>
        summarize(n_spots = n()) |>
        ggplot(aes(y = n_spots, x = position)) +
        # geom_violin() +
        geom_boxplot() +
        geom_point() +
        geom_line(aes(group = subject, color = subject)) +
        theme_bw()
)


# Layer proportion - Sp9 --------------------------------------------------

vars_spd <- c(
    "BayesSpace_harmony_09",
    "BayesSpace_harmony_16"
) |>
    set_names()

# TODO: create pooled figure
vars_spd |>
    map_dfr(
        .f = function(var) {
            fnl_dat |>
                select(position, spd = var, subject) |>
                group_by(position, spd) |>
                summarize(
                    n_sample = n_distinct(subject),
                    n_spots = n()
                )
        },
        .id = "spd_method"
    ) |>
    filter(n_sample != 10)

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

vars_spd |>
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
    }) |>
    group_split(spd_method, .keep = TRUE) |>
    walk(.f = function(sum_dat) {
        # browser()
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
                    # TODO: decide which order for color palette we use
                    Polychrome::palette36.colors(k)[str_sub(layer_df$layer_long, -2, -1) |>
                        as.integer()],
                    layer_df$layer_combo
                )
            ) +
            guides(x = guide_axis(angle = -45)) +
            facet_wrap(~position) +
            ylab("Proportion of Spots") +
            xlab("") +
            # theme(axis.title.x = element_blank()) +
            theme_set(theme_bw(base_size = 20))

        ggsave(
            filename = here(
                "plots/15_cell_composition",
                "layer_comp",
                paste0("layer_comp_sp", k, ".pdf")
            ),
            plot = ret_plot,
            height = 7,
            width = 12
        )
    })
