# Load Packages
library(tidyverse)
library(here)
library(SpatialExperiment)
library(ggpubr)
library(spatialLIBD)

# Load Spe Object
spe_dat <- readRDS(
    here(
        "processed-data/rdata/spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
    )
)

# Subset
sampleid <- "Br8667_mid"
spe <- spe_dat[, spe_dat$sample_id == "Br8667_mid"]


fnl_dat <- colData(spe) |> data.frame()


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





# Find Dominate Spots -----------------------------------------------------

dom_thres <- 0.5
res <- "broad"
deconv <- "tangram"
c_type_oi <- "excit"


deconv_count <- deconv_cell_counts |>
    data.frame() |>
    pull(paste(res, deconv, sep = "_"))

coi_perc <- fnl_dat |>
    dplyr::select(n_coi = paste(res, deconv, c_type_oi, sep = "_")) |>
    cbind(n_cell = deconv_count) |>
    mutate(
        perc_coi = n_coi / n_cell,
        include = perc_coi > dom_thres
    )
coi_perc <- coi_perc / deconv_count


coloc_spe <- calc_coloc(spe_dat, "EFNA5", "EPHA5", sample_id = "Br8667_mid")


tmp <- vis_coloc(coloc_spe[, which(coi_perc$include == TRUE)],
    "EFNA5", "EPHA5",
    sample_id = "Br8667_mid"
)
