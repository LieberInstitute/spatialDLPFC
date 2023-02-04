library("here")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")

spe_in <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
)

anno_samples <- c("Br6522_ant", "Br6522_mid", "Br8667_post")

anno_wrinkle_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF",
    "spatialLIBD_ManualAnnotation_{sample_id}_wrinkle.csv"
)

anno_layers_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF",
    "spatialLIBD_ManualAnnotation_{sample_id}_layers.csv"
)

spe <- readRDS(spe_in)

anno_list <- list()

for (sample_id in anno_samples) {
    #   Adjust paths for this sample
    this_anno_wrinkle_path <- sub(
        "\\{sample_id\\}", sample_id, anno_wrinkle_path
    )
    this_anno_layers_path <- sub("\\{sample_id\\}", sample_id, anno_layers_path)

    #   Read in wrinkle annotation and use unique and informative colnames
    anno_wrinkle <- this_anno_wrinkle_path |>
        read.csv() |>
        as_tibble() |>
        rename(wrinkle_type = ManualAnnotation, barcode = spot_name)

    #   Read in layer annotation and use unique and informative colnames
    anno_layers <- this_anno_layers_path |>
        read.csv() |>
        as_tibble() |>
        rename(manual_layer_label = ManualAnnotation, barcode = spot_name)

    anno_list[[sample_id]] <- anno_layers |>
        full_join(anno_wrinkle, by = c("sample_id", "barcode"))
}

#   Combine into a single tibble with all 3 samples
anno <- do.call(rbind, anno_list) |>
    as_tibble()

stopifnot(all(anno$barcode %in% colnames(spe)))

added_coldata <- tibble(
    "barcode" = colnames(spe), "sample_id" = spe$sample_id
) |>
    left_join(anno)

# table(is.na(added_coldata$manual_layer_label))
# table(is.na(added_coldata$wrinkle_type))

#   Add the wrinkle and layer annotation to the SPE's colData
colData(spe) <- cbind(
    colData(spe),
    added_coldata |> select(c("manual_layer_label", "wrinkle_type"))
)

#   Colnames started out sorted, so we'll mantain that organization
colData(spe) <- colData(spe)[, sort(colnames(colData(spe)))]
