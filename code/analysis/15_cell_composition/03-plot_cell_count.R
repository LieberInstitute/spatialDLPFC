library(tidyverse)
library(SpatialExperiment)
library(here)
library("Polychrome")
library(ggpubr)
data(palette36)

method_colors = palette36[6:8]
names(method_colors) = c('Tangram', 'Cell2location', 'SPOTlight')

spe_dat <- readRDS(
    here(
        "processed-data/rdata/spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
    )
)

fnl_dat <- colData(spe_dat) |> data.frame()

# Create Panel A ----------------------------------------------------------
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

# * Curate Data -------------------------------------------------------------
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


# * Generate a List of ggplots --------------------------------------------------------
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

# * Create Plot Grid --------------------------------------------------------
layer_grid_plot <- ggarrange(
    plotlist = layer_plots,
    ncol = 2, nrow = 1
)





# Create Panel B ----------------------------------------------------------

deconv_comb <- expand_grid(
    res = c("broad", "layer"),
    deconv = c("tangram", "cell2location", "spotlight")
)

deconv_df <- fnl_dat |>
    select(starts_with(c("broad", "layer")))

deconv_com_indx_mat <- deconv_comb |>
    pmap_dfc(.f = function(res, deconv) {
        str_starts(names(deconv_df), paste(res, deconv, sep = "_")) |>
            as.integer() |>
            data.frame() |>
            set_names(paste(res, deconv, sep = "_"))
    }) |>
    as.matrix()

# Check if the correct number of colums are detected
stopifnot(
    colSums(deconv_com_indx_mat) == ifelse(deconv_comb$res == "broad", 7, 13)
)

deconv_cell_counts <- (deconv_df |> as.matrix()) %*% deconv_com_indx_mat


#* Curate Data ------------------------------------------

cell_count_dat <- deconv_comb |>
    pmap_dfr(.f = function(res, deconv) {
        # factor_cell_type <- factor_cell_type_layer
        # cell_type_colors <- cell_type_colors_layer
        # if(res == "broad"){
        #     factor_cell_type <- factor_cell_type_broad
        #     cell_type_colors <- cell_type_colors_broad
        # }


        # Fetch the total cell counts per spot
        deconv_count <- deconv_cell_counts |>
            data.frame() |>
            pull(paste(res, deconv, sep = "_"))


        # browser()
        ret_plot <- fnl_dat |>
            dplyr::select(
                position, subject,
                starts_with(paste(res, deconv, sep = "_"))
            ) |>
            cbind(deconv_count) |>
            group_by(position, subject) |>
            summarise(
                n_cell_deconv = sum(deconv_count) # ,
                # across(starts_with(paste(res, deconv, sep = "_")), .fns = sum)
            ) |>
            ungroup() |>
            mutate(
                res = res,
                method = deconv
            )
    })

cell_count_plots <- cell_count_dat |>
    group_split(res) |>
    map(.f = function(dat) {
        res <- dat |>
            pull(res) |>
            head(1)
        ret_plot <- dat |>
            mutate(
                position = str_to_sentence(position),
                method = factor(method,
                    levels = c("cell2location", "spotlight", "tangram"),
                    labels = c("Cell2location", "SPOTlight", "Tangram")
                )
            ) |>
            ggplot(aes(x = subject, y = n_cell_deconv, fill = method)) +
            geom_col(position = "dodge") +
            facet_wrap(~position) +
            guides(
                x = guide_axis(angle = -45)
            ) +
            labs(
                title = paste(res, "level") |> str_to_title(),
                x = "",
                y = "Estimated Cell Counts",
                fill = "Method"
            ) +
            coord_cartesian(ylim = c(0, 50100)) +
            scale_y_continuous(
                breaks = (0:5) * 10000,
                labels = c("0", paste0((1:5) * 10, "K"))
            ) +
            scale_fill_manual(values = method_colors) +
            theme_set(theme_bw(base_size = 20)) +
            theme(
                plot.title = element_text(hjust = 0.5),
                axis.title.x = element_blank()
            )
    })


# * Create Plot Grid --------------------------------------------------------
cell_count_grid_plot <- ggarrange(
    plotlist = cell_count_plots,
    ncol = 2, nrow = 1 # ,
    # legend = "right",
    # common.legend = TRUE
)


ggsave(
    here(
        "plots/15_cell_composition/",
        "cell_count/",
        "cell_count_grid.pdf"
    ),
    plot = cell_count_grid_plot,
    height = 15,
    width = 12
)



# Supplementary Figure ----------------------------------------------------


sfigu_deconvo_visium_composition <- ggarrange(
    layer_grid_plot, cell_count_grid_plot,
    ncol = 1, nrow = 2,
    labels = "AUTO",
    font.label = list(
        size = 30, color = "black",
        face = "bold", family = NULL
    )
)

ggsave(
    here(
        "plots/15_cell_composition/",
        paste0("sfigu_deconvo_visium_composition.pdf")
    ),
    plot = sfigu_deconvo_visium_composition,
    height = 15,
    width = 26
)
