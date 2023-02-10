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
            select(position, subject,
                   spd = var, coloc
            )
    }, .id = "spd_method") |>
    # group_split(spd_method, position, .keep = TRUE) |>
    group_split(spd_method, spd, .keep = TRUE) |>
    map_dfr(.f = function(dat) {
        spd_method <- dat$spd_method |>
            unique()

        spd <- dat$spd |> unique()

        k <- dat$spd_method |>
            unique() |>
            str_sub(start = -2, end = -1) |>
            as.integer()

        # layer_df <- bayes_layers |>
        #     filter(Annotation == sprintf("k%02d", k))
        # browser()
        # Set up meta for this spd
        spd_df <- bayes_layers |>
            inner_join(
                data.frame(
                    # Annotation = sprintf("k%02d", k),
                    layer_long = sprintf("Sp%02dD%02d", k, spd)
                )
            )


        data.frame(
            spd_df,
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
            geom_bar(position = "fill", stat = "identity") +
            # scale_fill_manual(
            #     name = "BayesSpace Domain",
            #     values = set_names(
            #         Polychrome::palette36.colors(k)[str_sub(layer_df$layer_long, -2, -1) |>
            #                                             as.integer()],
            #         layer_df$layer_combo
            #     )
            # ) +
            guides(x = guide_axis(angle = 90)) +
            facet_wrap(~Annotation) +
            ylab("Proportion of Spots") +
            xlab("") +
            theme_set(theme_bw(base_size = 20))
    })


ggsave(
    filename = here(
        "plots/15_cell_composition",
        # "layer_comp",
        paste0("coloc_comp_spd.pdf")
    ),
    plot = ggpubr::ggarrange(
        nrow = 1,
        plotlist = ret_plot_list,
        common.legend = TRUE),
    width = 10
)

