library(tidyverse)
library(here)
library(SpatialExperiment)

spe_dat <- readRDS(here("processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters_and_deconvolution_results.rds"))

fnl_dat <- colData(spe_dat) |> data.frame()

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

n_pred_Neuropil <- c(
    map_dbl(deconv_cell_counts |> data.frame(), ~ sum(.x == 0)),
    VistoSeg = sum(fnl_dat$VistoSeg_count == 0)
)
# NOTE: cell2location produces expectation. So, wouldn't be exactly 0.
# see https://jhu-genomics.slack.com/archives/C01EA7VDJNT/p1670017513046739?thread_ts=1670009703.418879&cid=C01EA7VDJNT


deconv_comb |>
    pwalk(.f = function(res, deconv) {
        deconv_count <- deconv_cell_counts |>
            data.frame() |>
            pull(paste(res, deconv, sep = "_"))
        diff <- fnl_dat$VistoSeg_count - deconv_count

        ret_plot <- ggplot(fnl_dat) +
            geom_histogram(aes(x = diff)) +
            facet_grid(
                rows = vars(position),
                cols = vars(subject)
            ) +
            labs(
                title = paste0(
                    "Histogram of difference estimated cell count per spot",
                    "(VistoSeg (new) - ", paste(res, deconv, sep = "_"), ")."
                )
            ) +
            theme_bw()

        ggsave(
            here(
                "plots/15_cell_composition/",
                paste0(
                    "cell_count_diff_",
                    res, "_", deconv, ".pdf"
                )
            ),
            plot = ret_plot,
            height = 5,
            width = 10
        )
    })


# (Test) - Cell Proportion---------------------------------------------------------


#* (Test) - SpD Strat Cell Proportion to spd compo---------------------------------------------------------
# Hypothesis testing
# If spD stratified cell proportion is related to the spatial domain composition
# Find Colors

for (spd in vars_spd) {
    deconv_comb |>
        pwalk(.f = function(res, deconv) {
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

            spd_fct <- factor(fnl_dat[, spd])
            # Calculate cell composition per spatial domain

            # browser()
            cont_dat <- fnl_dat |>
                dplyr::select(
                    position, subject,
                    starts_with(paste(res, deconv, sep = "_"))
                ) |>
                cbind(deconv_count, spd_fct) |>
                group_by(position, subject, spd_fct) |>
                summarise(
                    n_cell_deconv = sum(deconv_count),
                    across(starts_with(paste(res, deconv, sep = "_")), .fns = sum)
                ) |>
                ungroup() |>
                mutate(across(starts_with(paste(res, deconv, sep = "_")),
                    .fns = ~ .x / n_cell_deconv
                )) |>
                # Calc centered logratio transformation
                # NOTE: doesn't work because some geometric means are 0
                # rowwise() |>
                # transmute(geo_mean = prod(c_across(starts_with(paste(res, deconv, sep = "_")))))
                mutate(across(starts_with(paste(res, deconv, sep = "_")),
                    .fns = ~ .x / broad_tangram_excit
                )) |>
                select(-n_cell_deconv, -broad_tangram_excit) |>
                dplyr::rename_with(
                    .cols = starts_with(paste(res, deconv, sep = "_")),
                    .fn = ~ str_remove(.x,
                        pattern = paste0(res, "_", deconv, "_")
                    )
                ) |>
                group_split(spd_fct) |>
                map(.f = function(spd_dat) {
                    browser()
                    mdl <- lm(cbind(astro, endomural, inhib, micro, oligo, opc) ~ position, data = spd_dat)
                    anova(mdl)
                })


            # dplyr::rename_with(.cols = starts_with(paste(res, deconv, sep = "_")),
            #                    .fn = ~str_remove(.x,
            #                                      pattern = paste0(res,"_", deconv, "_")
            #                    )
            # )

            # TODO:
        })
}


#* (Test) - SpD Strat Cell Proportion---------------------------------------------------------
## Hypothesis testing
# Find Colors

for (spd in vars_spd) {
    deconv_comb |>
        pwalk(.f = function(res, deconv) {
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

            spd_fct <- factor(fnl_dat[, spd])
            # Calculate cell composition per spatial domain

            # browser()
            cont_dat <- fnl_dat |>
                dplyr::select(
                    position, subject,
                    starts_with(paste(res, deconv, sep = "_"))
                ) |>
                cbind(deconv_count, spd_fct) |>
                group_by(position, subject, spd_fct) |>
                summarise(
                    n_cell_deconv = sum(deconv_count),
                    across(starts_with(paste(res, deconv, sep = "_")), .fns = sum)
                ) |>
                ungroup() |>
                mutate(across(starts_with(paste(res, deconv, sep = "_")),
                    .fns = ~ .x / n_cell_deconv
                )) |>
                # Calc centered logratio transformation
                # NOTE: doesn't work because some geometric means are 0
                # rowwise() |>
                # transmute(geo_mean = prod(c_across(starts_with(paste(res, deconv, sep = "_")))))
                mutate(across(starts_with(paste(res, deconv, sep = "_")),
                    .fns = ~ .x / broad_tangram_excit
                )) |>
                select(-n_cell_deconv, -broad_tangram_excit) |>
                dplyr::rename_with(
                    .cols = starts_with(paste(res, deconv, sep = "_")),
                    .fn = ~ str_remove(.x,
                        pattern = paste0(res, "_", deconv, "_")
                    )
                ) |>
                group_split(spd_fct) |>
                map(.f = function(spd_dat) {
                    browser()
                    mdl <- lm(cbind(astro, endomural, inhib, micro, oligo, opc) ~ position, data = spd_dat)
                    anova(mdl)
                })


            # dplyr::rename_with(.cols = starts_with(paste(res, deconv, sep = "_")),
            #                    .fn = ~str_remove(.x,
            #                                      pattern = paste0(res,"_", deconv, "_")
            #                    )
            # )

            # TODO:
        })
}

#* (Test) - Independency between cell comp and layer comp---------------------------------------------------------
for (spd in vars_spd) {
    deconv_comb |>
        # filter(deconv=="tangram") |>
        pwalk(.f = function(res, deconv) {
            # Fetch the total cell counts per spot
            deconv_count <- deconv_cell_counts |>
                data.frame() |>
                pull(paste(res, deconv, sep = "_"))

            spd_fct <- factor(fnl_dat[, spd])
            # Calculate cell composition per spatial domain

            # browser()

            layer_prob <- fnl_dat |>
                cbind(spd_fct) |>
                group_by(position, subject, spd_fct) |>
                summarize(n = n()) |>
                mutate(total = sum(n)) |>
                ungroup() |>
                transmute(position, subject, spd_fct, prop = n / total)


            cell_prob <- fnl_dat |>
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
                select(-n_cell_deconv)


            expect_dat <- layer_prob |>
                left_join(cell_prob,
                    by = c("position", "subject")
                ) |>
                mutate(across(starts_with(paste(res, deconv, sep = "_")),
                    .fns = ~ .x * prop
                )) |>
                select(-prop) |>
                pivot_longer(
                    cols = starts_with(paste(res, deconv, sep = "_")),
                    names_to = "cell_type",
                    values_to = "expect_prop"
                )

            obs_dat <- fnl_dat |>
                dplyr::select(
                    position, subject,
                    starts_with(paste(res, deconv, sep = "_"))
                ) |>
                cbind(deconv_count, spd_fct) |>
                group_by(position, subject, spd_fct) |>
                summarise(
                    n_cell_deconv = sum(deconv_count),
                    across(starts_with(paste(res, deconv, sep = "_")), .fns = sum)
                ) |>
                mutate(sample_n_cell = sum(n_cell_deconv)) |>
                ungroup() |>
                mutate(across(starts_with(paste(res, deconv, sep = "_")),
                    .fns = ~ .x / sample_n_cell
                )) |>
                select(-n_cell_deconv, -sample_n_cell) |>
                pivot_longer(
                    cols = starts_with(paste(res, deconv, sep = "_")),
                    names_to = "cell_type",
                    values_to = "obs_prop"
                )

            stopifnot(identical(dim(expect_dat), dim(obs_dat)))

            diff_dat <- inner_join(
                expect_dat, obs_dat
            ) |>
                transmute(
                    position,
                    id = glue::glue("{subject}_{spd_fct}"),
                    cell_type,
                    diff = obs_prop - expect_prop
                )

            if (!dir.exists(here("plots/15_cell_composition/tmp_cell_prop_diff"))) {
                dir.create(here("plots/15_cell_composition/tmp_cell_prop_diff"),
                    recursive = TRUE
                )
            }

            ggsave(
                here(
                    "plots/15_cell_composition/tmp_cell_prop_diff",
                    paste0(
                        "cell_diff_heatmap_",
                        res, "_", deconv, ".pdf"
                    )
                ),
                ggplot(diff_dat) +
                    geom_tile(aes(x = cell_type, y = id, fill = diff)) +
                    facet_wrap(~position) +
                    scale_fill_viridis(discrete = FALSE),
                height = 10,
                width = 5
            )
        })
}


# (Vis) Cell Proportion ---------------------------------------------------------

# Find Colors
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
#* (Vis) Overall Cell Proportion ------------------------------------------
deconv_comb |>
    pmap_dfr(.f = function(res, deconv) {
        # browser()
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
            # group_by(position, subject) |>
            summarise(
                n_cell_deconv = sum(deconv_count),
                across(starts_with(paste(res, deconv, sep = "_")), .fns = sum)
            ) |>
            # ungroup() |>
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
                res = res,
                method = deconv
            )
    }) |>
    group_split(res) |>
    walk(.f = function(dat) {
        # browser()
        res <- dat |>
            pull(res) |>
            head(1)
        cell_type_colors <- cell_type_colors_layer
        sce_cell_type <- sce$cellType_layer
        if (res == "broad") {
            cell_type_colors <- cell_type_colors_broad
            sce_cell_type <- sce$cellType_broad_hc
        }


        if ("drop" %in% levels(sce_cell_type)) {
            sce_cell_type <- droplevels(sce_cell_type[sce_cell_type != "drop"])
        }


        sce_count_df <- sce_cell_type |>
            table() |>
            proportions() |>
            data.frame() |>
            rename(
                sce_cell_type = "cell_type",
                Freq = "cell_perc"
            ) |>
            mutate(
                n_cell_deconv = length(sce_cell_type),
                res = res,
                method = "snRNA"
            )
        # browser()
        ret_plot <- dat |>
            rbind(sce_count_df) |>
            ggplot(aes(
                x = method,
                y = cell_perc, fill = cell_type
            )) +
            geom_col() +
            # geom_bar(position = "stack", stat = "identity") +
            # facet_wrap(~position) +
            scale_fill_manual(values = cell_type_colors) +
            coord_polar("y") +
            scale_x_discrete(
                # NOTE: the order is from the inside to outside
                limits = c("", "cell2location", "tangram", "snRNA")
            ) +
            theme_void() +
            theme(legend.position = "none")



        ggsave(
            here(
                "plots/15_cell_composition/",
                paste0(
                    "cell_comp_overall_",
                    res, "_no_legend.pdf"
                )
            ),
            plot = ret_plot # ,
            # height = 10,
            # width = 5
        )
    })


#* (Vis) Per-sample # of Cell ------------------------------------------

deconv_comb |>
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
                n_cell_deconv = sum(deconv_count) # ,
                # across(starts_with(paste(res, deconv, sep = "_")), .fns = sum)
            ) |>
            ungroup() |>
            # mutate(across(starts_with(paste(res, deconv, sep = "_")),
            #               .fns = ~.x/n_cell_deconv)) |>
            # pivot_longer(starts_with(paste(res, deconv, sep = "_")),
            #              names_to = "cell_type",
            #              values_to = "cell_perc") |>
            mutate(
                res = res,
                method = deconv
            )
    }) |>
    group_split(res) |>
    walk(.f = function(dat) {
        # browser()

        res <- dat |>
            pull(res) |>
            head(1)

        ret_plot <- dat |>
            ggplot(aes(x = subject, y = n_cell_deconv, fill = method)) +
            geom_col(
                position = "dodge" # , stat = "identity"
            ) +
            facet_wrap(~position) +
            # scale_fill_manual(values = cell_type_colors) +
            guides(x = guide_axis(angle = 90)) +
            # scale_y_continuous(trans='log2') +
            # scale_y_log10() +
            # labs(
            #     title = paste0(
            #         "Cell composition at ", res, " level using ", deconv
            #     )
            # ) +
            theme_bw()

        ggsave(
            here(
                "plots/15_cell_composition/",
                paste0(
                    "cell_num_",
                    res, ".pdf"
                )
            ),
            plot = ret_plot # ,
            # height = 10,
            # width = 5
        )
    })

#* (Vis) Per-sample Cell Proportion ------------------------------------------

deconv_comb |>
    pwalk(.f = function(res, deconv) {
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
            mutate(cell_type = str_remove(
                cell_type,
                paste0(res, "_", deconv, "_")
            ) |>
                factor_cell_type()) |>
            ggplot(aes(x = subject, y = cell_perc, fill = cell_type)) +
            geom_bar(position = "stack", stat = "identity") +
            facet_wrap(~position) +
            scale_fill_manual(values = cell_type_colors) +
            guides(x = guide_axis(angle = 90)) +
            labs(
                title = paste0(
                    "Cell composition at ", res, " level using ", deconv
                )
            ) +
            theme_bw()

        ggsave(
            here(
                "plots/15_cell_composition/",
                paste0(
                    "cell_comp_",
                    res, "_", deconv, ".pdf"
                )
            ),
            plot = ret_plot # ,
            # height = 10,
            # width = 5
        )
    })

#* (Vis) - SpD Stratified Cell Proportion ------------------------------------------
for (spd in vars_spd) {
    deconv_comb |>
        pwalk(.f = function(res, deconv) {
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

            spd_fct <- factor(fnl_dat[, spd])
            # Calculate cell composition per spatial domain

            # browser()
            ret_plot <- fnl_dat |>
                dplyr::select(
                    position, subject,
                    starts_with(paste(res, deconv, sep = "_"))
                ) |>
                cbind(deconv_count, spd_fct) |>
                group_by(position, subject, spd_fct) |>
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
                mutate(cell_type = str_remove(
                    cell_type,
                    paste0(res, "_", deconv, "_")
                ) |>
                    factor_cell_type()) |>
                ggplot(aes(x = subject, y = cell_perc, fill = cell_type)) +
                geom_bar(position = "stack", stat = "identity") +
                facet_grid(rows = vars(spd_fct), cols = vars(position)) +
                scale_fill_manual(values = cell_type_colors) +
                guides(x = guide_axis(angle = 90)) +
                labs(
                    title = paste0(
                        "Cell composition at ", res, " level using ", deconv
                    )
                ) +
                theme_bw()

            ggsave(
                here(
                    "plots/15_cell_composition/",
                    paste0(
                        "cell_comp_", spd, "_",
                        res, "_", deconv, ".pdf"
                    )
                ),
                plot = ret_plot,
                height = 10,
                width = 5
            )
        })
}
