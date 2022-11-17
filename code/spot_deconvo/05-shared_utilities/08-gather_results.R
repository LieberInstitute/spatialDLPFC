library("here")
library("jaffelab")
library("tidyverse")
library("sessioninfo")

cell_group <- "broad" # "broad" or "layer"
dataset <- "IF" # "IF" or "nonIF"

for (cell_group in c("broad", "layer")) {
    for (dataset in c("IF", "nonIF")) {
        print(paste("Gathering results for", cell_group, dataset, "..."))
        
        ########################################################################
        #   Variable definitions
        ########################################################################
        
        sample_ids_path <- here(
            "processed-data", "spot_deconvo", "05-shared_utilities", dataset,
            "sample_ids.txt"
        )
        
        deconvo_tools <- c("01-tangram", "03-cell2location", "04-spotlight")
        deconvo_tool_names = c('tangram', 'cell2location', 'SPOTlight')
        
        #   "Ground-truth" cell counts from cellpose + trained classification
        #   tree (IF only)
        actual_paths <- here(
            "processed-data", "spot_deconvo", "02-cellpose", "{sample_id}",
            "clusters.csv"
        )
        
        #   Cell counts estimated by each deconvolution tool
        observed_paths <- here(
            "processed-data", "spot_deconvo", "{deconvo_tool}", dataset,
            cell_group, "{sample_id}", "clusters.csv"
        )
        
        raw_results_path <- here(
            "processed-data", "spot_deconvo", "05-shared_utilities", dataset,
            paste0("results_raw_", cell_group, ".csv")
        )
        
        collapsed_results_path <- here(
            "processed-data", "spot_deconvo", "05-shared_utilities", dataset,
            paste0("results_collapsed_", cell_group, ".csv")
        )
        
        cell_types_actual <- c("astro", "micro", "neuron", "oligo", "other")
        if (cell_group == "broad") {
            cell_types <- c(
                "Astro", "EndoMural", "Excit", "Inhib", "Micro", "Oligo", "OPC"
            )
        } else {
            cell_types <- c(
                "Astro", "EndoMural", "Excit_L2_3", "Excit_L3", "Excit_L3_4_5",
                "Excit_L4", "Excit_L5", "Excit_L5_6", "Excit_L6", "Inhib",
                "Micro", "Oligo", "OPC"
            )
        }
        
        ########################################################################
        #   Read in and aggregate cell counts across samples and methods
        ########################################################################
        
        sample_ids <- readLines(sample_ids_path)
        
        actual_list <- list()
        observed_list <- list()
        index <- 1
        
        added_colnames <- c("barcode", "sample_id", "deconvo_tool", "obs_type")
        
        #   Read in counts for all cell types for all samples and deconvolution
        #   tools (as well as the "ground-truth", for IF)
        for (sample_id in sample_ids) {
            if (dataset == "IF") {
                #   Read in the "ground-truth" counts for this sample
                actual_path <- sub("\\{sample_id\\}", sample_id, actual_paths)
                actual_df_small <- read.csv(actual_path)
                
                #   Make sure we have just counts for each cell type, barcode,
                #   and sample ID variables-- nothing else. Order column names
                actual_df_small$barcode <- ss(actual_df_small$key, "_", 1)
                actual_df_small$sample_id <- sample_id
                actual_df_small$obs_type <- "actual"
            }
            
            for (deconvo_tool in deconvo_tools) {
                #   Read in estimated cell counts for this deconvo tool and
                #   sample
                observed_path <- sub(
                    "\\{sample_id\\}", sample_id, observed_paths
                )
                observed_path <- sub(
                    "\\{deconvo_tool\\}", deconvo_tool, observed_path
                )
                observed_df_small <- read.csv(observed_path)
                colnames(observed_df_small) = gsub(
                    '\\.', '_', colnames(observed_df_small)
                )
                
                #   Make sure column names are consistent and include only info
                #   about barcode, sample_id, deconvo tool, and cell-type counts
                observed_df_small$barcode <- ss(observed_df_small$key, "_", 1)
                observed_df_small$sample_id <- sample_id
                observed_df_small$deconvo_tool <- deconvo_tool
                observed_df_small$obs_type <- "observed"
                observed_df_small <- observed_df_small[
                    , c(added_colnames, cell_types)
                ]
                
                if (dataset == "IF") {
                    actual_df_small$deconvo_tool <- deconvo_tool
                    
                    actual_list[[index]] <- actual_df_small[
                        , c(added_colnames, cell_types_actual)
                    ]
                }
                
                observed_list[[index]] <- observed_df_small
                index <- index + 1
            }
        }
        
        #   Form data frame containing cell-type counts for all spots, samples,
        #   and deconvolution tools
        observed_df <- as_tibble(do.call(rbind, observed_list))
        
        #   Write raw results to CSV
        write.csv(
            observed_df[, -match('obs_type', colnames(observed_df))],
            file = raw_results_path, quote = FALSE, row.names = FALSE
        )
        
        #   For IF, collapse cell types for comparison against the "ground
        #   truth". Then merge those ground-truth counts for export to a single
        #   CSV
        if (dataset == "IF") {
            actual_df <- as_tibble(do.call(rbind, actual_list))
            
            colnames(observed_df) <- tolower(colnames(observed_df))
            
            if (cell_group == "broad") {
                observed_df <- observed_df %>%
                    mutate(
                        "neuron" = excit + inhib,
                        "oligo" = oligo + opc,
                        "other" = endomural
                    ) %>%
                    select(all_of(c(added_colnames, cell_types_actual)))
            } else {
                observed_df = observed_df %>%
                    rowwise() %>%
                    mutate(
                        "neuron" = sum(
                            c_across(starts_with(c('excit_', 'inhib')))
                        ),
                        "oligo" = oligo + opc,
                        "other" = endomural
                    ) %>%
                    ungroup() %>%
                    select(all_of(c(added_colnames, cell_types_actual)))
            }
            
            stopifnot(all(colnames(observed_df) == colnames(actual_df)))
            
            #   Write collapsed results to CSV
            write.csv(
                rbind(observed_df, actual_df), file = collapsed_results_path,
                quote = FALSE, row.names = FALSE
            )
        }
    }
}

session_info()
