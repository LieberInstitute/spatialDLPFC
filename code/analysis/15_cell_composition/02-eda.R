library(tidyverse)
library(SingleCellExperiment)

# NOTE: this line only relevents to Boyi
# .libPaths(c(.libPaths(), "/users/bguo/R/4.2"))

# Question about dissection variation
## Are the number of spots different across pos
## Are the percentage of each domain change across pos
## Per pos, is the cell composition of layers the same
## Cell composition change across samples

# TODO: edit this path
if(!"fnl_dat" %in% ls()) fnl_dat <- readRDS("~/cell_comp_full_dat.rds")

# Find Colors
sce <- readRDS("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC_annotated/se.rds")
cell_type_colors_layer <- metadata(sce)$cell_type_colors_layer[levels(sce$cellType_layer)]

factor_cell_type_layer <- function(vec){
    factor(vec,
           levels = c("Astro",  "EndoMural", "Micro", "Oligo", "OPC",
                      "Excit_L2_3", "Excit_L3", "Excit_L3_4_5", "Excit_L4","Excit_L5",
                      "Excit_L5_6","Excit_L6","Inhib"),
           labels = c("Astro",  "EndoMural", "Micro", "Oligo", "OPC",
                      "Excit_L2/3", "Excit_L3", "Excit_L3/4/5", "Excit_L4","Excit_L5",
                      "Excit_L5/6","Excit_L6","Inhib")
    )
}

# # of Samples per pos ------------------------------------------------------------

fnl_dat |>
    group_by(region) |>
    summarize(n_sample = n_distinct(subject),
              n_spots = n(),
              n_sp9_miss = sum(is.na(sp9)),
              n_sp16_miss = sum(is.na(sp16)))
## Hypothesis testing: Are the number of spots different across pos
## NOTE: no great power
spot_lm <- fnl_dat |>
    group_by(subject, region) |>
    summarize(n_spots = n()) |>
    with(lm(n_spots~region))
summary(spot_lm)
anova(spot_lm)


## Are the cell counts matches
# TODO: move this cell count to data prep step
cell_count_fnl_dat <- fnl_dat |>
    rowwise() |>
    mutate(n_cell_deconv = sum(c_across(Astro:OPC))) |>
    ungroup() |>
    dplyr::rename(n_cell_VS = count)
#NOTE: count should be the cell count from vistoseg
# see https://jhu-genomics.slack.com/archives/D04BRDW7Q8Y/p1669149653581219



# The discrepency of cell counts estimated
ggsave("~/cell_comp/cell_count_diff.pdf",
       plot = cell_count_fnl_dat |>
           mutate(diff = n_cell_VS - n_cell_deconv) |>
           ggplot() +
           geom_histogram(aes(x = diff)) +
           ylim(-10, 10) +
           facet_grid(rows = vars(region),
                      cols = vars(subject)) +
           labs(title = "count-n_cell_deconv"),
       height = 5,
       width = 10
)
# TODO: Possibly highlight which have matching snRNA


cell_count_fnl_dat |>
    transmute(n_cell_VS, n_cell_deconv,
              diff = n_cell_VS - n_cell_deconv) |>
    summary()
# NOTE: new hyposthesis: using paired snRNA can imporve the deconvolution
# NOTE: Q: is it possible to see which are badly aligned spots?


## Are there any spots missing spd labels?
fnl_dat |>
    group_by(region, subject) |>
    summarize(n_spots = n(),
              n_sp9_miss = sum(is.na(sp9)),
              n_sp16_miss = sum(is.na(sp16)),
              avg_sp9_miss = n_sp9_miss/n_spots,
              avg_sp16_miss = n_sp16_miss/n_spots) |>
    filter(avg_sp9_miss!=0)
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

ggsave("~/cell_comp/cell_comp_overall.pdf",
       plot =
           cell_count_fnl_dat |>
           group_by(region, subject) |>
           summarize(
               across(Astro:OPC, sum),
               n_cell_deconv = sum(n_cell_deconv)) |>
           ungroup() |>
           mutate(across(Astro:OPC,  ~.x/n_cell_deconv)) |>
           pivot_longer(Astro:OPC, names_to = "cell_type", values_to = "cell_perc") |>
           mutate(cell_type = factor_cell_type_layer(cell_type)) |>
           ggplot(aes( x = subject, y = cell_perc, fill = cell_type)) +
           geom_bar(position = "stack", stat = "identity") +
           facet_wrap(~region) +
           #TODO: check if the color mapping is correct
           scale_fill_manual(values = cell_type_colors_layer) +
           theme_bw()
)
# NOTE: Overall cell_composition doesn't seems to be very different
# TODO: combine total number of cells and total number of spots into this plot


# Sp9 Analysis ------------------------------------------------------------

# If all sample have 9 layers
fnl_dat |>
    group_by(region, sp9) |>
    summarize(n_sample = n_distinct(subject),
              n_spots = n()) |>
    filter(n_sample!=10)
# NOTE: ALl sample have 9 layers


## TODO: Are the number of spots different across pos

# Within each pieace of the brain, is layer composition the same across indviduals
sp9_comp_pc_smp <- fnl_dat |>
    group_split(region, .keep = TRUE) |>
    map(.f = function(dat){
        # browser()
        dat |>
            group_by(subject, sp9) |>
            summarize(n_spots = n()) |>
            mutate(tot_spots = sum(n_spots),
                   layer_perc = n_spots/tot_spots) |>
            # pivot_wider(id_cols = subject,
            #             names_from = sp9,
            #             values_from = layer_perc) |>
            mutate(region = head(dat$region,1))
    })

# NOTE: Visually, the composition of Sp9D are pretty different


# TODO: bueatify this
ggsave(filename = "~/cell_comp/layer_comp_sp19.pdf",
       sp9_comp_pc_smp |>
           bind_rows() |>
           ggplot(aes(x = subject, y = layer_perc, fill = sp9)) +
           geom_bar(position = "fill", stat = "identity") +
           facet_wrap(~region))


## Data table
sp9_comp_pc_smp |>
    map_dfr(.f = function(dat){
        dat |>
            pivot_wider(id_cols = subject,
                        names_from = sp9,
                        values_from = layer_perc) |>
            mutate(region = head(dat$region,1))
    })



# Cell composition of Sp9D

ggsave("~/cell_comp/cell_comp_sp9.pdf",
       plot =
           cell_count_fnl_dat |>
           group_by(region, sp9, subject) |>
           summarise(
               n_cell_deconv = sum(n_cell_deconv),
               across(Astro:OPC, sum)) |>
           ungroup() |>
           mutate(across(Astro:OPC, ~.x/n_cell_deconv)) |>
           pivot_longer(Astro:OPC, names_to = "cell_type", values_to = "cell_perc") |>
           mutate(cell_type = factor_cell_type_layer(cell_type)) |>
           ggplot(aes( x = subject, y = cell_perc, fill = cell_type)) +
           geom_bar(position = "stack", stat = "identity") +
           facet_grid(rows = vars(sp9), cols = vars(region)) +
           scale_fill_manual(values = cell_type_colors_layer) +
           theme_bw(),
       height = 10, width = 5
)

# TODO: possibly horizontal plot
# TODO: deconv bash effect


# Density of cells


# Sp16 Analysis ------------------------------------------------------------
# If all sample have 16 layers
fnl_dat |>
    group_by(region, sp16) |>
    summarize(n_sample = n_distinct(subject),
              n_spots = n()) |>
    filter(n_sample!=10)

setdiff(fnl_dat$subject |> unique(),
        fnl_dat |> filter(region == "anterior", sp16 == 9) |>
            pull(subject) |> unique())
# NOTE: Br6471_ant doesn't have Layer 9

# Within each pieace of the brain, is layer composition the same across indviduals
sp16_comp_pc_smp <- fnl_dat |>
    group_split(region, .keep = TRUE) |>
    map(.f = function(dat){
        # browser()
        dat |>
            group_by(subject, sp16) |>
            summarize(n_spots = n()) |>
            mutate(tot_spots = sum(n_spots),
                   layer_perc = n_spots/tot_spots) |>
            # pivot_wider(id_cols = subject,
            #             names_from = sp9,
            #             values_from = layer_perc) |>
            mutate(region = head(dat$region,1))
    })

# NOTE: Visually, the composition of Sp16D are pretty different
# TODO:make heatmap


# TODO: bueatify this
ggsave(filename = "~/cell_comp/layer_comp_sp16.pdf",
       sp16_comp_pc_smp |>
           bind_rows() |>
           ggplot(aes(x = subject, y = layer_perc, fill = sp16)) +
           geom_bar(position = "fill", stat = "identity") +
           facet_wrap(~region))


## Data table
sp16_comp_pc_smp |>
    map_dfr(.f = function(dat){
        dat |>
            pivot_wider(id_cols = subject,
                        names_from = sp16,
                        values_from = layer_perc) |>
            mutate(region = head(dat$region,1))
    })


ggsave("~/cell_comp/cell_comp_sp16.pdf",
       plot =
           cell_count_fnl_dat |>
           group_by(region, sp16, subject) |>
           summarise(
               n_cell_deconv = sum(n_cell_deconv),
               across(Astro:OPC, sum)) |>
           ungroup() |>
           mutate(across(Astro:OPC, ~.x/n_cell_deconv)) |>
           pivot_longer(Astro:OPC, names_to = "cell_type", values_to = "cell_perc") |>
           mutate(cell_type = factor_cell_type_layer(cell_type)) |>
           ggplot(aes( x = subject, y = cell_perc, fill = cell_type)) +
           geom_bar(position = "stack", stat = "identity") +
           facet_grid(rows = vars(sp16), cols = vars(region)) +
           scale_fill_manual(values = cell_type_colors_layer) +
           theme_bw(),
       height = 25, width = 10
)


