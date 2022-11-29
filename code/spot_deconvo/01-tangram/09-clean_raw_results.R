library("here")
library("jaffelab")
library("tidyverse")
library("sessioninfo")

#   We'd like to share the tangram results such that the corresponding SPE
#   object can be loaded and the tangram results trivially added with
#   spatialLIBD::cluster_import. Write a CSV given the original combined
#   results, in the appropriate format

for (dataset in c("IF", "nonIF")) {
    for (cell_group in c("broad", "layer")) {
        raw_results_path <- here(
            "processed-data", "spot_deconvo", "05-shared_utilities", dataset,
            paste0("results_raw_", cell_group, ".csv")
        )

        clean_results_path <- here(
            "processed-data", "spot_deconvo", "01-tangram", dataset, cell_group,
            "raw_results", "clusters.csv"
        )

        dir.create(
            dirname(clean_results_path),
            showWarnings = FALSE, recursive = TRUE
        )

        #   Read in just the tangram results
        raw_results <- read.csv(raw_results_path) |>
            filter(deconvo_tool == "tangram")

        #   Add the 'key' column and remove unnecessary columns
        raw_results$key <- paste(
            raw_results$barcode, raw_results$sample_id,
            sep = "_"
        )
        raw_results <- raw_results[
            , -match(
                c("barcode", "sample_id", "deconvo_tool"), colnames(raw_results)
            )
        ]

        #   Write to disk
        write.csv(
            raw_results,
            file = clean_results_path, row.names = FALSE,
            quote = FALSE
        )
    }
}

session_info()
