
setwd("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/")

library(tidyverse)
library(SingleCellExperiment)
# library(readr)
library(SpatialExperiment)
library(unglue)
library(here)



# Commond Line Receiver ---------------------------------------------------
# Receive Arguments from Command Line
# args=(commandArgs(TRUE))
# if(length(args)==0){
#     print("No arguments supplied.")
# }else{
#     for(i in 1:length(args)){
#         eval(parse(text=args[[i]]))
#     }
# }

## NOTE: reserved for when command line arguments fails
# deconv_method <- c("tangram", "cell2location", "SPOTlight")
# deconv_method <- deconv_method[1]


# Data Prep ---------------------------------------------------------------

deconvo_res_path <- here(
    "processed-data/spot_deconvo/05-shared_utilities/nonIF/",
    "results_raw_layer.csv"
)


deconv_res <- read.table(deconvo_res_path,sep = ",", header  = TRUE)


# Wide format data
deconv_res_wide <- deconv_res |>
    rowwise() |>
    mutate(n_cell = sum(c_across(Astro:OPC))) |>
    ungroup() |>
    pivot_wider(
        names_from = deconvo_tool,
        values_from = c(Astro:OPC, n_cell)
    )

# method_dat <- read.table(deconvo_res_path,sep = ",", header  = TRUE) |>
#     filter(deconvo_tool == deconv_method) |>
#     mutate(new_key = paste(barcode, sample_id, sep = "_"))


load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")

n_sample <- nrow(spe)

spe_col <- colData(spe) |> data.frame(
    check.names = FALSE
) |>
    mutate(sp9 = factor(bayesSpace_harmony_9, levels = 1:9, labels = paste("Sp9D", 1:9, sep="")),
           sp16 = factor(bayesSpace_harmony_16, levels = 1:16, labels = paste("Sp16D", 1:16, sep=""))
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

# Some sanity check
# tmp <- unglue::unglue_data(rownames(colData(spe)),
#                            pattern = "{bc_trim}{bc_tail=\\.\\d*|$}") |>
#     filter(bc_tail!="")

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

fnl_dat <- deconv_res_wide |>
    inner_join(
        fnl_col,
        by = c("barcode", "sample_id")
    )

nrow(fnl_dat) == ncol(spe)

#NOTE: count should be the cell count from vistoseg
# see https://jhu-genomics.slack.com/archives/D04BRDW7Q8Y/p1669149653581219
fnl_dat <- fnl_dat |>
    dplyr::rename(n_cell_vs = count)

saveRDS(fnl_dat,
        file = here("processed-data/15_cell_composition/cell_comp_full_dat.rds")
)




