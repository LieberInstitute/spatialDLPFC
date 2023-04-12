# Load Packages
library(tidyverse)
library(here)
library(SpatialExperiment)
library(ggpubr)
library(SingleCellExperiment)

# Load Spe Object
spe_dat <- readRDS(
    here(
        "processed-data/rdata/spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
    )
)

# sampleid <- "Br8667_mid"
# spe <- spe_dat[,spe_dat$sample_id == "Br8667_mid"]

# Config Color ------------------------------------------------------------
sce <- readRDS("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC_annotated/se.rds")
# cell_type_colors_layers were removed from metadata(sce) [as of 2023/04/11]
# cell_type_colors_layer <- metadata(sce)$cell_type_colors_layer[levels(sce$cellType_layer)]
cell_type_colors_broad <- metadata(sce)$cell_type_colors_broad[levels(sce$cellType_broad_hc)]

# Hard Code of cell type color
cell_type_colors_layer <- c(
    "Astro" = "#3BB273",
    "EndoMural" = "#FF56AF",
    "Micro" = "#663894",
    "Oligo" = "#E07000",
    "OPC" = "#D2B037",
    "Excit_L2/3" = "#D0D1E6",
    "Excit_L3" = "#A6BDDB",
    "Excit_L3/4/5" ="#74A9CF",
    "Excit_L4" = "#3690C0",
    "Excit_L5" = "#0570B0",
    "Excit_L5/6" = "#045A8D",
    "Excit_L6" = "#023858",
    "Excit_ambig" = "#A0A7A7",
    "Inhib" = "#E94F37"
)


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

# Calc Colocalization -----------------------------------------------------
source(here("code/analysis/15_cell_composition/fun_coloc.R"))
# coloc_spe <- calc_coloc(spe_dat, "EFNA5", "EPHA5", sample_id = "Br8667_mid")
coloc_spe <- calc_coloc(spe_dat, "EFNA5", "EPHA5", sample_id = NULL)

fnl_dat <- colData(coloc_spe) |> data.frame()

# Calculate Total Number of Cells per spot --------------------------------
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
        # NOTE: Br 2720_ant Doesn't have colocalized spots
        # Run this after the ret_plot
        # setdiff( unique(fnl_dat$sample_id),
        #          ret_plot |>
        #              filter( coloc == "co-localize") |>
        #              pull(sample_id) |> unique()
        # )

        ret_plot <- fnl_dat |>
            dplyr::select(
                sample_id, coloc,
                starts_with(paste(res, deconv, sep = "_"))
            ) |>
            cbind(deconv_count) |>
            group_by(sample_id, coloc) |>
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
                ) |>
                    factor_cell_type(),
                coloc = coloc,
                res = res,
                method = deconv
            )
    })


group_cell_comp_dat <- cell_comp_dat |>
    filter(method == "cell2location", res == "layer")

# res <- dat |> pull(res) |> head(1)
dat <- group_cell_comp_dat
viol_plot <- dat |>
    mutate(
        cell_type = fct_relevel(cell_type,
                                "Inhib",
                                after = Inf
        ),
        coloc = factor(coloc,
                       levels = c("co-localize", "EFNA5", "EPHA5", "Neither"),
                       labels = c("EFNA5 & EPHA5", "EFNA5 only", "EPHA5 only", "Neither")
        )
    ) |>
    # Remove this to make full
    # dplyr::filter(cell_type %in% c("Excit_L5/6", "Excit_L6") )) |>
    mutate(cell_type = fct_drop(cell_type)) |>
    # mutate(position = str_to_sentence(position),
    #        method = factor(method,
    #                        levels = c("cell2location", "spotlight", "tangram"),
    #                        labels = c("Cell2location", "SPOTlight", "Tangram"))
    # ) |>
    ggplot(
        aes(x = coloc, y = cell_perc, fill = cell_type)
    ) +
    # geom_boxplot(show.legend = FALSE) +
    geom_violin(show.legend = FALSE) +
    # geom_bar(position = "stack", stat = "identity") +
    facet_wrap(~cell_type, scale = "free_y", ncol = 3) +
    scale_fill_manual(
        name = "Cell Type",
        # TODO: edit this
        limits = names(cell_type_colors_layer),
        values = cell_type_colors_layer
    ) +
    guides(x = guide_axis(angle = 90)) +
    labs(
        #     title = paste(method |> str_to_title(),
        #                   "at",
        #                   paste(res, "Resolution") |> str_to_title()),
        y = "Predicted Proportion in Spots",
        x = "Spot Co-expression"
    ) +
    scale_y_continuous(n.breaks = 3) +
    theme_set(theme_bw(base_size = 20))

ggsave(
    filename = here(
        "plots/15_cell_composition",
        # "layer_comp",
        paste0("coloc_comp_all_cells_c2l.pdf")
    ),
    plot = ret_plot,
    height = 10,
    width = 8
)



# * Compositional Plot ----------------------------------------------------

comp_plot <- dat |>
    mutate(
        cell_type = fct_relevel(cell_type,
                                "Inhib",
                                after = Inf
        ),
        coloc = factor(coloc,
                       levels = c("co-localize", "EFNA5", "EPHA5", "Neither"),
                       labels = c("EFNA5 & EPHA5", "EFNA5 only", "EPHA5 only", "Neither")
        )
    ) |>
    # Remove this to make full
    # dplyr::filter(cell_type %in% c("Excit_L5/6", "Excit_L6") )) |>
    mutate(cell_type = fct_drop(cell_type)) |>
    group_by(coloc, cell_type) |>
    summarize(cell_perc = median(cell_perc)) |>
    ggplot(
        aes(x = coloc, y = cell_perc, fill = cell_type)
    ) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(
        name = "Cell Type",
        # TODO: edit this
        limits = names(cell_type_colors_layer),
        values = cell_type_colors_layer,
        guide = "none"
    ) +
    guides(x = guide_axis(angle = 90)) +
    labs(
        #     title = paste(method |> str_to_title(),
        #                   "at",
        #                   paste(res, "Resolution") |> str_to_title()),
        y = "Proportion Medians",
        x = ""
        # x = "Spot Co-expression"
    ) +
    # scale_y_continuous(n.breaks = 3) +
    theme_set(theme_bw(base_size = 20))

# ggsave(
#     filename = here(
#         "plots/15_cell_composition",
#         # "layer_comp",
#         paste0("coloc_comp_all_cells_c2l_median.pdf")
#     ),
#     plot = ret_plot,
#     height = 10,
#     width = 8
# )


ggsave(
    filename = here(
        "plots/15_cell_composition",
        paste0("coloc_comp_combined.pdf")
    ),
    plot = ggpubr::ggarrange(
        viol_plot, comp_plot,
        nrow = 1,
        labels = "AUTO",
        font.label = list(size = 20,
                          color = "black",
                          face = "bold",
                          family = "Helvetica"),
        widths = c(2, 1)
    ),
    height = 8,
    width = 10
)



# Tests -------------------------------------------------------------------
dat |>
    group_split(cell_type, .keep = TRUE) |>
    map_dfr(.f = function(cell_dat) {
        stopifnot(length(unique(cell_dat$cell_type)) == 1)

        test_dat <- cell_dat |>
            transmute(
                coloc_group = factor(coloc == "co-localize",
                                     levels = c("TRUE", "FALSE")
                ),
                cell_perc
            )
        # browser()
        data.frame(
            cell_type = unique(cell_dat$cell_type),
            wilcox = wilcox.test(cell_perc ~ coloc_group, test_dat,
                                 alternative = "greater"
            )$p.value,
            t_test = t.test(cell_perc ~ coloc_group, test_dat,
                            alternative = "greater"
            )$p.value
        )
    })
