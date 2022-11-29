library(tidyverse)
library(readr)
library(SpatialExperiment)
library(unglue)

deconvo_res_path <- paste0(
    "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/",
    "processed-data/spot_deconvo/05-shared_utilities/nonIF/",
    "results_raw_layer.csv"
)


tangram_res <- read.table(deconvo_res_path,sep = ",", header  = TRUE) |>
    filter(deconvo_tool == "tangram") |>
    mutate(new_key = paste(barcode, sample_id, sep = "_"))

n_tangram <- nrow(tangram_res)
# Load spe of all samples ------------------------------------------------------

# The clustering raw results have inconsistent naming convention,
# Cause problem when merging
# https://github.com/LieberInstitute/spatialDLPFC/issues/133
# cluster_fld_path <- paste(
#     "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC",
#     "processed-data", "rdata", "spe", "clustering_results", sep = "/"
# )
#
# sp9_path <- paste(cluster_fld_path,
#                   "bayesSpace_harmony_9",
#                   "clusters.csv", sep = "/")
# sp9 <- read.csv(sp9_path, header = TRUE)
#
#
# sp9_path <- paste(cluster_fld_path,
#                   "bayesSpace_harmony_9",
#                   "clusters.csv", sep = "/")
# sp9 <- read.csv(sp9_path, header = TRUE) |>
#     dplyr::rename(sp9 = cluster)
#
# # Code to debug sp9 cluster loading
# # tmp <- sp9 |> unglue_unnest(col = key,
# #                      patterns = "{spot_id}_{Br_num}_{pos}")
# #
# # tmp |> filter(Br_num == "Br2720")
#
# sp16_path <- paste(cluster_fld_path,
#                   "bayesSpace_harmony_9",
#                   "clusters.csv", sep = "/")
# sp16 <- read.csv(sp16_path, header = TRUE)|>
#     dplyr::rename(sp16 = cluster)

load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")

n_sample <- nrow(spe)

spe_col <- colData(spe) |> data.frame(
    check.names = FALSE
) |>
    mutate(sp9 = factor(paste("Sp9D", bayesSpace_harmony_9, sep="")),
           sp16 = factor(paste("Sp16D", bayesSpace_harmony_16, sep=""))
    ) |>
    dplyr::select(key:array_col,
                  sum_umi:col,
                  sp9, sp16) |>
    rownames_to_column("bc_tmp") |>
    unglue_unnest(col = key,
                  pattern = "{key_bc}_{key_sampleID}_{key_sec=ant|mid|post}{key_tail}",
                  remove = FALSE) |>
    unglue_unnest(col=bc_tmp,
                  pattern = "{bc_trim}{bc_tail=\\.\\d*|$}",
                  remove = FALSE)

tmp <- unglue::unglue_data(rownames(colData(spe)),
                           pattern = "{bc_trim}{bc_tail=\\.\\d*|$}") |>
    filter(bc_tail!="")

fnl_col <- spe_col |>
    # Confirming the information are matched for each spots
    filter(bc_trim == key_bc,
           paste(key_sampleID, key_sec, sep = "_") == sample_id) |>
    dplyr::rename(barcode = bc_trim) |>
    dplyr::select(-starts_with("key_"),
                  -starts_with("bc_")) |>
    mutate(new_key = paste(barcode, sample_id, sep="_"))

# Dimensionality staill matches after data wrangling
nrow(fnl_col) == ncol(spe)


# rm(spe)

fnl_dat <- tangram_res |>
    inner_join(
        fnl_col,
        by = c("barcode", "sample_id")
    )

nrow(fnl_dat) == ncol(spe)

saveRDS(fnl_dat,
        file = "~/cell_comp_full_dat.rds")




