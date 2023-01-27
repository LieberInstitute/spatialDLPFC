#   Split deconvolution results by sample and add coordinates of each spot, in
#   preparation for importing the results as features to loopybrowser.com
#   using the python API

library("here")
library("jaffelab")
library("tidyverse")
library("sessioninfo")

dataset = "IF"
cell_group = "layer"

deconvo_tools <- c("tangram", "cell2location", "SPOTlight")

sample_ids_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", dataset,
    "sample_ids.txt"
)

collapsed_results_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", dataset,
    paste0("results_collapsed_", cell_group, ".csv")
)

out_dir = here("processed-data", "spot_deconvo", "05-shared_utilities", dataset)

spot_path = here(
    'processed-data', '01_spaceranger_IF', '{sample_id}', 'outs', 'spatial',
    'tissue_positions_list.csv'
)

cell_types = c("astro", "micro", "neuron", "oligo", "other")

sample_ids = readLines(sample_ids_path)

collapsed_results = read.csv(collapsed_results_path) |> as_tibble()

for (sample_id in sample_ids) {
    #   Read in the tissue positions file to get coordinates for each barcode
    #   on the full-resolution image
    tissue_positions = sub('\\{sample_id\\}', sample_id, spot_path) |>
        read.csv(header = FALSE) |>
        as_tibble() |>
        select(c(1, 5, 6))

    colnames(tissue_positions) = c("barcode", "x", "y")
    
    #   Sanity check: all barcodes from the results should be in the spaceranger
    #   output
    collapsed_barcodes = collapsed_results |>
        filter(sample_id == {{ sample_id }}) |>
        pull(barcode)
    
    stopifnot(all(collapsed_barcodes %in% (tissue_positions |> pull(barcode))))
    
    loopy_results = collapsed_results |>
        filter(
            sample_id == {{ sample_id }}, obs_type == "observed"
        ) |>
        left_join(tissue_positions, by = "barcode") |>
        pivot_longer(
            cols = all_of(cell_types), names_to = "cell_type",
            values_to = "count"
        ) |>
        select(c("x", "y", "deconvo_tool", "cell_type", "count"))
        
    write.csv(loopy_results, file.path(out_dir, paste0("loopy_", sample_id, ".csv")))
}

session_info()
