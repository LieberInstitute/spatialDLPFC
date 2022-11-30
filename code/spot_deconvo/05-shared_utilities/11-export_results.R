library("here")
library("jaffelab")
library("tidyverse")
library("sessioninfo")

#   We'd like to share the deconvo results such that the corresponding SPE
#   object can be loaded and the results trivially added with
#   spatialLIBD::cluster_import. Write CSVs for each deconvo tool and dataset
#   (IF and nonIF) given the original combined results

deconvo_tools = c("01-tangram", "03-cell2location", "04-spotlight")
names(deconvo_tools) = c("tangram", "cell2location", "SPOTlight")
    
for (dataset in c("IF", "nonIF")) {
    for (cell_group in c("broad", "layer")) {
        raw_results_path <- here(
            "processed-data", "spot_deconvo", "05-shared_utilities", dataset,
            paste0("results_raw_", cell_group, ".csv")
        )
        
        for (deconvo_tool in names(deconvo_tools)) {
            clean_results_path <- here(
                "processed-data", "spot_deconvo", deconvo_tools[deconvo_tool],
                dataset, cell_group, "raw_results", "clusters.csv"
            )
            
            dir.create(
                dirname(clean_results_path),
                showWarnings = FALSE, recursive = TRUE
            )
            
            #   Read in the results for this deconvo tool, cell group, and
            #   dataset
            raw_results <- read.csv(raw_results_path) |>
                filter(deconvo_tool == {{ deconvo_tool }})
            
            #   Add the 'key' column and remove unnecessary columns
            raw_results$key <- paste(
                raw_results$barcode, raw_results$sample_id,
                sep = "_"
            )
            raw_results <- raw_results[
                , -match(
                    c("barcode", "sample_id", "deconvo_tool"),
                    colnames(raw_results)
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
}

session_info()
