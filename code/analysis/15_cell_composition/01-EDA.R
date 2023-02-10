library(tidyverse)
library(SingleCellExperiment)
library(SpatialExperiment)
library(here)
# NOTE: this line only relevents to Boyi
# .libPaths(c(.libPaths(), "/users/bguo/R/4.2"))

# Question about dissection variation
## Are the number of spots different across pos
## Are the percentage of each domain change across pos
## Per pos, is the cell composition of layers the same
## Cell composition change across samples

# TODO: edit this path
# if(!"fnl_dat" %in% ls()) fnl_dat <- readRDS(here("processed-data/15_cell_composition/cell_comp_full_dat.rds"))

spe_dat <- readRDS(here("processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters_and_deconvolution_results.rds"))


# Find Colors
sce <- readRDS("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC_annotated/se.rds")
cell_type_colors_layer <- metadata(sce)$cell_type_colors_layer[levels(sce$cellType_layer)]

sp9_colors <- Polychrome::palette36.colors(9)
names(sp9_colors) <- paste0("Sp9D", 1:9)

sp16_colors <- Polychrome::palette36.colors(16)
names(sp16_colors) <- paste0("Sp16D", 1:16)

factor_cell_type_layer <- function(vec) {
    factor(vec,
        levels = c(
            "Astro", "EndoMural", "Micro", "Oligo", "OPC",
            "Excit_L2_3", "Excit_L3", "Excit_L3_4_5", "Excit_L4", "Excit_L5",
            "Excit_L5_6", "Excit_L6", "Inhib"
        ),
        labels = c(
            "Astro", "EndoMural", "Micro", "Oligo", "OPC",
            "Excit_L2/3", "Excit_L3", "Excit_L3/4/5", "Excit_L4", "Excit_L5",
            "Excit_L5/6", "Excit_L6", "Inhib"
        )
    )
}

# # of Samples per pos ------------------------------------------------------------

fnl_dat |>
    group_by(region) |>
    summarize(
        n_sample = n_distinct(subject),
        n_spots = n(),
        n_sp9_miss = sum(is.na(sp9)),
        n_sp16_miss = sum(is.na(sp16))
    )
## Hypothesis testing: Are the number of spots different across pos
## NOTE: no great power
spot_lm <- fnl_dat |>
    group_by(subject, region) |>
    summarize(n_spots = n()) |>
    with(lm(n_spots ~ region))
summary(spot_lm)
anova(spot_lm)

ggsave(
    filename = here("plots/15_cell_composition/", "n_spots.pdf"),
    fnl_dat |>
        group_by(subject, region) |>
        summarize(n_spots = n()) |>
        ggplot(aes(y = n_spots, x = region)) +
        # geom_violin() +
        geom_boxplot() +
        geom_point() +
        geom_line(aes(group = subject, color = subject)) +
        theme_bw()
)



## Are the cell counts matches
# cell_count_fnl_dat <- fnl_dat |>
#     group_split()
#     rowwise() |>
#     mutate(n_cell_deconv = sum(c_across(Astro:OPC))) |>
#     ungroup() |>
#     dplyr::rename(n_cell_VS = count)

# TODO: move this up
deconv_method <- c(
    "tangram", "cell2location" # ,
    # "SPOTlight"
)

deconv_method |>
    walk(.f = function(dc_mthd) {
        # browser()
        ggsave(
            here(
                "plots/15_cell_composition/",
                paste0("cell_count_diff_", dc_mthd, ".pdf")
            ),
            plot = fnl_dat |>
                dplyr::rename(n_cell_deconv = paste0("n_cell_", dc_mthd)) |>
                mutate(diff = n_cell_vs - n_cell_deconv) |>
                ggplot() +
                geom_histogram(aes(x = diff)) +
                facet_grid(
                    rows = vars(region),
                    cols = vars(subject)
                ) +
                labs(title = paste0("n_cell_VS-n_cell_", dc_mthd)) +
                theme_bw(),
            height = 5,
            width = 10
        )
    })
# The discrepency of cell counts estimated

deconv_method <- c(
    "tangram", "cell2location" # ,
    # "SPOTlight"
)

c( # "tangram", "cell2location"#,
    "SPOTlight"
) |>
    walk(.f = function(dc_mthd) {
        # browser()
        ggsave(
            here(
                "plots/15_cell_composition/",
                paste0("cell_count_diff_zoomin_", dc_mthd, ".pdf")
            ),
            plot = fnl_dat |>
                dplyr::rename(n_cell_deconv = paste0("n_cell_", dc_mthd)) |>
                mutate(diff = n_cell_vs - n_cell_deconv) |>
                ggplot() +
                geom_histogram(aes(x = diff)) +
                facet_grid(
                    rows = vars(region),
                    cols = vars(subject)
                ) +
                xlim(c(-10, 50)) +
                labs(title = paste0("n_cell_VS-n_cell_", dc_mthd)) +
                theme_bw(),
            height = 5,
            width = 10
        )
    })



c("tangram", "cell2location", "SPOTlight") |>
    map_dfc(.f = function(dc_mthd) {
        # browser()
        fnl_dat |>
            dplyr::rename(n_cell_deconv = paste0("n_cell_", dc_mthd)) |>
            transmute(n_cell_deconv,
                diff = n_cell_vs - n_cell_deconv
            ) |>
            rename_all(~ paste0(.x, "_", dc_mthd)) |>
            sapply(summary) |>
            cbind()
    })

# NOTE: new hyposthesis: using paired snRNA can imporve the deconvolution
# NOTE: Q: is it possible to see which are badly aligned spots?


## Are there any spots missing spd labels?
fnl_dat |>
    group_by(region, subject) |>
    summarize(
        n_spots = n(),
        n_sp9_miss = sum(is.na(sp9)),
        n_sp16_miss = sum(is.na(sp16)),
        avg_sp9_miss = n_sp9_miss / n_spots,
        avg_sp16_miss = n_sp16_miss / n_spots
    ) |>
    filter(avg_sp9_miss != 0)
## NOTE: all spots have spd labels



## Are the percentage of each domain change across pos
# fnl_dat |>
#     group_by(region) |>
#     summarize(n_spots = n(),
#               n_sp9_miss = sum(is.na(sp9)),
#               n_sp16_miss = sum(is.na(sp16)),
#               avg_sp9_miss = n_sp9_miss/n_spots,
#               avg_sp16_miss = n_sp16_miss/n_spots) |>
#     filter(avg_sp9_miss!=0)



## Per pos, is the cell composition of layers the same

## Cell composition change across samples

c("tangram", "cell2location", "SPOTlight") |>
    walk(.f = function(dc_mthd) {
        # browser()
        # str_offset <- str_length(paste0("_", dc_mthd))
        ggsave(
            here(
                "plots/15_cell_composition/",
                paste0("cell_comp_overall_", dc_mthd, ".pdf")
            ),
            plot = fnl_dat |>
                select(
                    region, subject,
                    ends_with(dc_mthd)
                ) |>
                rename_with(
                    .cols = ends_with(dc_mthd),
                    .f = ~ str_remove(.x, paste0("_", dc_mthd))
                ) |>
                group_by(region, subject) |>
                summarize(
                    across(Astro:OPC, sum),
                    n_cell_deconv = sum(n_cell)
                ) |>
                ungroup() |>
                mutate(across(Astro:OPC, ~ .x / n_cell_deconv)) |>
                pivot_longer(Astro:OPC, names_to = "cell_type", values_to = "cell_perc") |>
                mutate(cell_type = factor_cell_type_layer(cell_type)) |>
                ggplot(aes(x = subject, y = cell_perc, fill = cell_type)) +
                geom_bar(position = "stack", stat = "identity") +
                facet_wrap(~region) +
                # TODO: check if the color mapping is correct
                scale_fill_manual(values = cell_type_colors_layer) +
                guides(x = guide_axis(angle = 90)) +
                labs(title = dc_mthd) +
                theme_bw()
        )
    })

# ggsave(here("plots/15_cell_composition/","cell_comp_overall.pdf"),
#        plot =
#            cell_count_fnl_dat |>
#            group_by(region, subject) |>
#            summarize(
#                across(Astro:OPC, sum),
#                n_cell_deconv = sum(n_cell_deconv)) |>
#            ungroup() |>
#            mutate(across(Astro:OPC,  ~.x/n_cell_deconv)) |>
#            pivot_longer(Astro:OPC, names_to = "cell_type", values_to = "cell_perc") |>
#            mutate(cell_type = factor_cell_type_layer(cell_type)) |>
#            ggplot(aes( x = subject, y = cell_perc, fill = cell_type)) +
#            geom_bar(position = "stack", stat = "identity") +
#            facet_wrap(~region) +
#            #TODO: check if the color mapping is correct
#            scale_fill_manual(values = cell_type_colors_layer) +
#            guides(x =  guide_axis(angle = 90)) +
#            theme_bw()
# )
# NOTE: Overall cell_composition doesn't seems to be very different


# Sp9 Analysis ------------------------------------------------------------

# If all sample have 9 layers
fnl_dat |>
    group_by(region, sp9) |>
    summarize(
        n_sample = n_distinct(subject),
        n_spots = n()
    ) |>
    filter(n_sample != 10)
# NOTE: ALl sample have 9 layers


## TODO: Are the number of spots different across pos

# Within each pieace of the brain, is layer composition the same across indviduals
sp9_comp_pc_smp <- fnl_dat |>
    group_split(region, .keep = TRUE) |>
    map(.f = function(dat) {
        # browser()
        dat |>
            group_by(subject, sp9) |>
            summarize(n_spots = n()) |>
            mutate(
                tot_spots = sum(n_spots),
                layer_perc = n_spots / tot_spots
            ) |>
            # pivot_wider(id_cols = subject,
            #             names_from = sp9,
            #             values_from = layer_perc) |>
            mutate(region = head(dat$region, 1))
    })

# NOTE: Visually, the composition of Sp9D are pretty different




# TODO: bueatify this
ggsave(
    filename = here("plots/15_cell_composition/", "layer_comp_sp9.pdf"),
    sp9_comp_pc_smp |>
        bind_rows() |>
        ggplot(aes(x = subject, y = layer_perc, fill = sp9)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = sp9_colors) +
        guides(x = guide_axis(angle = 90)) +
        facet_wrap(~region) +
        theme_bw()
)


## Data table
sp9_comp_pc_smp |>
    map_dfr(.f = function(dat) {
        dat |>
            pivot_wider(
                id_cols = subject,
                names_from = sp9,
                values_from = layer_perc
            ) |>
            mutate(region = head(dat$region, 1))
    })



# Cell composition of Sp9D

ggsave(here("plots/15_cell_composition/", "cell_comp_sp9.pdf"),
    plot =
        cell_count_fnl_dat |>
            group_by(region, sp9, subject) |>
            summarise(
                n_cell_deconv = sum(n_cell_deconv),
                across(Astro:OPC, sum)
            ) |>
            ungroup() |>
            mutate(across(Astro:OPC, ~ .x / n_cell_deconv)) |>
            pivot_longer(Astro:OPC, names_to = "cell_type", values_to = "cell_perc") |>
            mutate(cell_type = factor_cell_type_layer(cell_type)) |>
            ggplot(aes(x = subject, y = cell_perc, fill = cell_type)) +
            geom_bar(position = "stack", stat = "identity") +
            facet_grid(rows = vars(sp9), cols = vars(region)) +
            scale_fill_manual(values = cell_type_colors_layer) +
            guides(x = guide_axis(angle = 90)) +
            theme_bw(),
    height = 10, width = 5
)

# Sp16 Analysis ------------------------------------------------------------
# If all sample have 16 layers
fnl_dat |>
    group_by(region, sp16) |>
    summarize(
        n_sample = n_distinct(subject),
        n_spots = n()
    ) |>
    filter(n_sample != 10)

setdiff(
    fnl_dat$subject |> unique(),
    fnl_dat |> filter(region == "anterior", sp16 == "Sp16D9") |>
        pull(subject) |> unique()
)
# NOTE: Br6471_ant doesn't have Layer 9

# Within each pieace of the brain, is layer composition the same across indviduals
sp16_comp_pc_smp <- fnl_dat |>
    group_split(region, .keep = TRUE) |>
    map(.f = function(dat) {
        # browser()
        dat |>
            group_by(subject, sp16) |>
            summarize(n_spots = n()) |>
            mutate(
                tot_spots = sum(n_spots),
                layer_perc = n_spots / tot_spots
            ) |>
            # pivot_wider(id_cols = subject,
            #             names_from = sp9,
            #             values_from = layer_perc) |>
            mutate(region = head(dat$region, 1))
    })

# NOTE: Visually, the composition of Sp16D are pretty different

ggsave(
    filename = here("plots/15_cell_composition/", "layer_comp_sp16.pdf"),
    sp16_comp_pc_smp |>
        bind_rows() |>
        ggplot(aes(x = subject, y = layer_perc, fill = sp16)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = sp16_colors) +
        guides(x = guide_axis(angle = 90)) +
        # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        facet_wrap(~region) +
        theme_bw()
)


## Data table
sp16_comp_pc_smp |>
    map_dfr(.f = function(dat) {
        dat |>
            pivot_wider(
                id_cols = subject,
                names_from = sp16,
                values_from = layer_perc
            ) |>
            mutate(region = head(dat$region, 1))
    })


ggsave(here("plots/15_cell_composition/", "cell_comp_sp16.pdf"),
    plot =
        cell_count_fnl_dat |>
            group_by(region, sp16, subject) |>
            summarise(
                n_cell_deconv = sum(n_cell_deconv),
                across(Astro:OPC, sum)
            ) |>
            ungroup() |>
            mutate(across(Astro:OPC, ~ .x / n_cell_deconv)) |>
            pivot_longer(Astro:OPC, names_to = "cell_type", values_to = "cell_perc") |>
            mutate(cell_type = factor_cell_type_layer(cell_type)) |>
            ggplot(aes(x = subject, y = cell_perc, fill = cell_type)) +
            geom_bar(position = "stack", stat = "identity") +
            facet_grid(rows = vars(sp16), cols = vars(region)) +
            scale_fill_manual(values = cell_type_colors_layer) +
            guides(x = guide_axis(angle = 90)) +
            theme_bw(),
    height = 25, width = 10
)
