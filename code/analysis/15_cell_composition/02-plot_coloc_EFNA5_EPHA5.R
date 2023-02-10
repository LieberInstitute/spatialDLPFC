library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(here)

spe_dat <- readRDS(here(
    "processed-data/rdata/spe/01_build_spe",
    "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
))

# Config ------------------------------------------------------------------

gene_1 <- "EFNA5"
gene_2 <- "EPHA5"

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

wilcox.test(Freq ~ L6_group, test_dat) # p=3.993e-09
t.test(Freq ~ L6_group, test_dat) # p=1.699e-05


# Plot --------------------------------------------------------------------

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
        return(c(0.0, 0.1, 0.2))
    } else {
        return(c(0.00, 0.05, 0.10))
    }
}



ret_plot <- plot_dat |>
    mutate(
        Var1 = factor(Var1,
            labels = c(
                paste0(gene_1, " & ", gene_2, " (p=1.7e-05)"),
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


# ret_plot_list <- plot_dat |>
#     group_split(Var1, .keep = TRUE) |>
#     map(.f=function(dat){
#         # browser()
#         group_var <- unique(dat$Var1)
#
#         ret <- dat |>
#             ggplot(
#                 aes(x = layer_combo, y = Freq, fill = Var1)
#             ) +
#             geom_boxplot(show.legend = FALSE) +
#             facet_wrap(~Var1,  scale = "free_y", ncol = 1) +
#             scale_fill_manual(
#                 values = c("#d7191c", "#abd9e9", "#2c7bb6") |>
#                     set_names(c(paste0(gene_1, " & ", gene_2),
#                                 paste0(gene_1, " only"),
#                                 paste0(gene_2, " only")))
#                 #set_names(
#                 # Polychrome::palette36.colors(k)[str_sub(layer_df$layer_long, -2, -1) |>
#                 #                                     as.integer()],
#                 # layer_df$layer_combo
#                 # )
#             ) +
#             guides(x = guide_axis(angle = 90)) +
#             ylab("") +
#             xlab("") +
#             theme_set(theme_bw(base_size = 20)) +
#             theme(panel.grid.minor.y = element_blank(),
#                    axis.title.x = element_blank(),
#                    axis.title.y = element_blank())
#
#         if(group_var %in% levels(group_var)[1:(nlevels(group_var)-1)])
#             ret <- ret +
#             theme(axis.text.x = element_blank())
#
#         return(ret)
#     })


# ret_plot <- ggpubr::ggarrange(
#     plotlist = ret_plot_list,
#     ncol = 1,
#     label.y = "Proportion of Spots"
# )


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


# Plot Cell Composition ---------------------------------------------------
# * Config Color ------------------------------------------------------------
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

fnl_dat <- colData(full_coloc_spe) |> data.frame()

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
ret_plot <- dat |>
    mutate(
        cell_type = fct_relevel(cell_type,
            "Inhib",
            after = Inf
        ),
        coloc = factor(coloc,
            levels = c("co-localize", gene_1, gene_2, "Neither"),
            labels = c(
                paste0(gene_1, " & ", gene_2),
                paste0(gene_1, " only"),
                paste0(gene_2, " only"),
                "Neither"
            )
        )
    ) |>
    dplyr::filter(cell_type %in% c("Excit_L5/6", "Excit_L6")) |>
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
    facet_wrap(~cell_type, # scale = "free_y",
        nrow = 1, ncol = 2
    ) +
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
        y = "Proportion of Cells",
        x = ""
    ) +
    theme_set(theme_bw(base_size = 20)) +
    theme(
        # plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank() # ,
        #       plot.margin = margin(0.1,0.9,0.1,0.3, "cm")
    )

ggsave(
    filename = here(
        "plots/15_cell_composition",
        # "layer_comp",
        paste0("coloc_comp_c2l.pdf")
    ),
    plot = ret_plot,
    height = 7, # 3.75*5,
    width = 5 # 7.5*5
)
