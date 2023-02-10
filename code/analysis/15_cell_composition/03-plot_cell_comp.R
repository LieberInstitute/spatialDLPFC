# Config Color ------------------------------------------------------------
sce <- readRDS("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC_annotated/se.rds")
cell_type_colors_layer <- metadata(sce)$cell_type_colors_layer[levels(sce$cellType_layer)]
cell_type_colors_broad <- metadata(sce)$cell_type_colors_broad[levels(sce$cellType_broad_hc)]


factor_cell_type_broad <- function(vec) {
    factor(vec,
        levels = c("astro", "endomural", "micro", "oligo", "opc", "excit", "inhib"),
        labels = c("Astro", "EndoMural", "Micro", "Oligo", "OPC", "Excit", "Inhib")
    )
}

factor_cell_type_layer <- function(vec) {
    factor(vec,
        # levels = c("Astro",  "EndoMural", "Micro", "Oligo", "OPC",
        #            "Excit_L2_3", "Excit_L3", "Excit_L3_4_5", "Excit_L4","Excit_L5",
        #            "Excit_L5_6","Excit_L6","Inhib"),
        levels = c(
            "astro", "endomural", "micro", "oligo", "opc",
            "excit_l2_3", "excit_l3", "excit_l3_4_5", "excit_l4", "excit_l5",
            "excit_l5_6", "excit_l6", "inhib"
        ),
        labels = c(
            "Astro", "EndoMural", "Micro", "Oligo", "OPC",
            "Excit_L2/3", "Excit_L3", "Excit_L3/4/5", "Excit_L4", "Excit_L5",
            "Excit_L5/6", "Excit_L6", "Inhib"
        )
    )
}


# Curate Data -------------------------------------------------------------


cell_comp_dat <- deconv_comb |>
    pmap_dfr(.f = function(res, deconv) {
        factor_cell_type <- factor_cell_type_layer
        cell_type_colors <- cell_type_colors_layer
        if (res == "broad") {
            factor_cell_type <- factor_cell_type_broad
            cell_type_colors <- cell_type_colors_broad
        }


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
                n_cell_deconv = sum(deconv_count),
                across(starts_with(paste(res, deconv, sep = "_")), .fns = sum)
            ) |>
            ungroup() |>
            mutate(across(starts_with(paste(res, deconv, sep = "_")),
                .fns = ~ .x / n_cell_deconv
            )) |>
            pivot_longer(starts_with(paste(res, deconv, sep = "_")),
                names_to = "cell_type",
                values_to = "cell_perc"
            ) |>
            mutate(
                cell_type = str_remove(
                    cell_type,
                    paste0(res, "_", deconv, "_")
                ) # |>
                # factor_cell_type()
                ,
                position = str_to_sentence(position),
                res = res,
                method = deconv
            )
    })

group_cell_comp_dat <- cell_comp_dat |>
    filter(method %in% c("tangram", "cell2location")) |>
    group_by(res, method)


cell_comp_plots <- group_cell_comp_dat |>
    group_split() |>
    map(.f = function(dat) {
        res <- unique(dat$res)
        method <- unique(dat$method)

        factor_cell_type <- factor_cell_type_layer
        cell_type_colors <- cell_type_colors_layer
        if (res == "broad") {
            factor_cell_type <- factor_cell_type_broad
            cell_type_colors <- cell_type_colors_broad
        }

        ret_plot <-
            dat |>
            mutate(cell_type = factor_cell_type(cell_type)) |>
            ggplot(aes(x = subject, y = cell_perc, fill = cell_type)) +
            geom_bar(position = "stack", stat = "identity") +
            facet_wrap(~position) +
            scale_fill_manual(
                name = "Cell Type",
                values = cell_type_colors
            ) +
            guides(x = guide_axis(angle = -45)) +
            labs(
                title = paste(
                    method |> str_to_title(),
                    "at",
                    paste(res, "Resolution") |> str_to_title()
                ),
                y = "Proportion of Counts",
                x = ""
            ) +
            theme_set(theme_bw(base_size = 20)) +
            theme(
                plot.title = element_text(hjust = 0.5),
                axis.title.x = element_blank(),
                plot.margin = margin(0.1, 0.9, 0.1, 0.3, "cm")
            )
    })

group_cell_comp_dat |>
    group_keys()



# ggsave(here("plots/15_cell_composition/",
#             paste0("cell_comp_",
#                    res, "_",deconv,".pdf")),
#        plot = ret_plot,
#        height = 7,
#        width = 12
# )
# })

broad_plots <- ggarrange(
    plotlist = cell_comp_plots[2:1],
    ncol = 1, nrow = 2,
    legend = "right",
    common.legend = TRUE,
    labels = "AUTO",
    font.label = list(
        size = 30, color = "black",
        face = "bold", family = NULL
    )
)

layer_plots <- ggarrange(
    plotlist = cell_comp_plots[4:3],
    ncol = 1, nrow = 2,
    legend = "right",
    common.legend = TRUE # ,
    # labels = "AUTO",
    # font.label = list(size = 30, color = "black",
    #                   face = "bold", family = NULL)
)

# Create Plot Grid --------------------------------------------------------


sfigu_deconvo_visium_composition_results <- ggarrange(
    broad_plots, layer_plots,
    # plotlist = cell_comp_plots,
    ncol = 2, nrow = 1 # ,
    # labels = "AUTO",
    # font.label = list(size = 30, color = "black",
    #                   face = "bold", family = NULL)
)

ggsave(
    here(
        "plots/15_cell_composition/",
        paste0("sfigu_deconvo_visium_composition_results.pdf")
    ),
    plot = sfigu_deconvo_visium_composition_results,
    height = 15,
    width = 26
)
