library(tidyverse)
library(SpatialExperiment)
library(here)
library("Polychrome")
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

fnl_dat |>
    group_by(region) |>
    summarize(
        n_sample = n_distinct(subject),
        n_spots = n(),
        n_sp9_miss = sum(is.na(sp9)),
        n_sp16_miss = sum(is.na(sp16))
    )

# (Vis) - Per-spot cell count estimation ------------------------------------------

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


#* (Vis) Per-sample # of Cell ------------------------------------------

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


# Create Individual Plots -------------------------------------------------




# Create Plot Grid --------------------------------------------------------
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
