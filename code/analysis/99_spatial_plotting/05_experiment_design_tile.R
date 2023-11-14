library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")
# library("patchwork")

plot_dir <-
    here(
        "plots",
        "99_spatial_plotting",
        "05_experiment_design_tile"
    )
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
}

## position standards
pos_df <- tibble(Position = c("anterior", "middle", "posterior"),
                 pos = c("Ant", "Mid", "Post"))

#### single nuc data ####
## no local access
load(here("processed-data", "rdata", "spe", "12_spatial_registration_sn" ,"sce_DLPFC.Rdata"), verbose = TRUE)

# Load full data
load(
  here(
    "processed-data",
    "rdata",
    "spe",
    "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
  )
)

sce_pd <- as.data.frame(colData(sce))

sn_samples <- sce_pd |>
    dplyr::count(Sample, SAMPLE_ID, Position, pos, round, BrNum, age, sex)

sn_n_samp <- sn_samples |>
    dplyr::count(Sample, Position, BrNum) |>
    mutate(data_type = "snRNA-seq",
           Position = tolower(Position))


#### spatial data ####
spe_sample_info <- read.csv(here("processed-data", "rdata","spe", "01_build_spe", "spe_sample_info.csv"))

spe_n_samp <- spe_sample_info |>
    mutate(Sample = gsub("_2","",sample_id)) |>
    dplyr::count(Sample, BrNum = subjects, Position = regions, data_type = "Visium")


#### Visium IF data ####
spe_IF_sample_info <- read.csv(here("processed-data", "rdata","spe_IF", "01_build_spe_IF", "spe_IF_sample_info.csv"))

spe_IF_n_samp <- spe_IF_sample_info |>
    separate(subject, into = c("BrNum", "pos", NA), sep = "_") |>
    mutate(Sample = paste0(BrNum,"_", tolower(pos)),
           data_type = "Visium-SPG") |>
    left_join(pos_df) |>
    dplyr::count(Sample, BrNum , Position, data_type)


#### ALL data ####

all_dlpfc <- spe_n_samp |>
    bind_rows(sn_n_samp) |>
    bind_rows(spe_IF_n_samp)

## Do we have matched assays for each sample?
matched <- all_dlpfc |> dplyr::count(Sample)|> mutate(Matched = n >= 2)

sn_pair <- sn_n_samp |>
    left_join(pos_df) |>
    group_by(BrNum) |>
    summarize(pos_pair = factor(paste(pos, collapse = ","),
                                levels = c("Ant", "Ant,Post","Ant,Mid", "Mid,Post"))) |>
    mutate(BrNum_order = fct_reorder(BrNum, -as.integer(pos_pair))) |>
    arrange(BrNum_order)


all_dlpfc_detail <- all_dlpfc |>
    left_join(matched |> select(-n)) |>
    left_join(pos_df) |>
    left_join(sn_pair) |>
    mutate(data_type = factor(data_type, levels = c("Visium", "snRNA-seq", "Visium-SPG")))

experiment_tile <- all_dlpfc_detail |>
    mutate(n = as.factor(n)) |>
    ggplot(aes(pos, BrNum, fill = n))+
    geom_tile()+
    geom_text(aes(label = n, color = Matched))+
    scale_color_manual(values =c(`FALSE` = "red", `TRUE` = "black")) +
    facet_wrap(~data_type, nrow = 1) +
    theme_bw()

ggsave(experiment_tile, filename = here(plot_dir, "spatialDLPFC_experiment_tile.png"))


experiment_tile2 <- all_dlpfc_detail |>
    mutate(n = as.factor(n)) |>
    ggplot(aes(pos, BrNum_order, fill = data_type))+
    geom_tile(color = "grey25")+
    geom_text(aes(label = n, color = Matched))+
    # geom_text(aes(label = n), color = "black")+
    scale_color_manual(values =c(`TRUE` = "black", `FALSE` = "grey75"), name ="Additional Assay") +
    scale_fill_manual(values =c( Visium = "#E16365", `snRNA-seq` = "#A967A9", `Visium-SPG` = "#66A7DB"), name ="Additional Assay") +
    facet_wrap(~data_type, nrow = 1) +
    theme_bw() +
    labs(x = "Position", y = "BrNum") +
    guides(fill = FALSE)

ggsave(experiment_tile2, filename = here(plot_dir, "spatialDLPFC_experiment_tile2.png"))
ggsave(experiment_tile2, filename = here(plot_dir, "spatialDLPFC_experiment_tile2.pdf"))


## add more detail with color
other_data <- all_dlpfc_detail |>
    group_by(pos, BrNum) |>
    summarize(n = n(),
              matched_data = case_when(any(data_type =="snRNA-seq") ~"snRNA-seq",
                                    any(data_type == "IF") ~ "IF",
                                    all(data_type == "Visium") ~ "Visium only",
                                    TRUE ~ "ERROR"))

experiment_tile_other <- other_data |>
    mutate(n = as.factor(n)) |>
    ggplot(aes(pos, BrNum, fill = matched_data))+
    geom_tile(color = "grey50")+
    geom_text(aes(label = n), color = "black")+
    theme_bw() +
    labs(x = "Position")

ggsave(experiment_tile_other, filename = here(plot_dir, "spatialDLPFC_experiment_tile_other.png"))

