#   For the Visium-IF data, add cell-type counts to the SpatialExperiment
#   object and save a new copy. This includes counts from each deconvolution
#   tool as well as the "ground-truth" counts from cellpose + the trained
#   DecisionTreeClassifier.

library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")

cell_groups <- c("broad", "layer")
deconvo_tools <- c("tangram", "cell2location", "SPOTlight")

sample_ids <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "sample_ids.txt"
)

spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)

spe_IF_out <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF", "spe.rds"
)

#   CSV of all IF counts at one resolution estimated by all methods 
raw_results_paths <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "results_raw_{cell_group}.csv"
)

#   CSV including CART-calculated counts
collapsed_results_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "results_collapsed_broad.csv"
)

################################################################################
#   Add spot deconvo results to the SPE object
################################################################################

spe <- readRDS(spe_IF_in)

added_coldata = tibble(
    'barcode' = colnames(spe), 'sample_id' = spe$sample_id
)

extra_colnames =  c('barcode', 'sample_id', 'deconvo_tool')

for (cell_group in cell_groups) {
    #   Read in cell counts
    raw_results = sub("\\{cell_group\\}", cell_group, raw_results_paths) |>
        read.csv() |>
        as_tibble()
    
    #   Individually add cell counts for each method
    for (deconvo_tool in deconvo_tools) {
        these_results = raw_results |>
            filter(deconvo_tool == {{ deconvo_tool }})
        
        #   Take the columns that represent cell types and append the cell group
        #   and deconvo tool to the colnames
        added_indices = match(extra_colnames, colnames(these_results))
        colnames(these_results)[- added_indices] = tolower(
            paste(
                cell_group, deconvo_tool,
                colnames(these_results)[- added_indices], sep = "_"
            )
        )
        
        #   Add columns to colData(spe)
        added_coldata_temp = added_coldata |>
            left_join(these_results, by = c('barcode', 'sample_id')) |>
            select(- all_of(extra_colnames))
        
        if(any(is.na(added_coldata_temp))) {
            stop("Some cell counts were NA after merging with SPE.")
        }
        
        colData(spe) = cbind(colData(spe), added_coldata_temp)
    }
}

#   Read in CART counts
these_results = read.csv(collapsed_results_path) |>
    as_tibble() |>
    filter(
        obs_type == "actual",
        deconvo_tool == deconvo_tools[1]
    ) |>
    select(- obs_type)

#   Append 'cart' to cell-type column names
added_indices = match(extra_colnames, colnames(these_results))
colnames(these_results)[- added_indices] = paste(
    'cart', colnames(these_results)[- added_indices], sep = "_"
)

#   Add columns to colData(spe)
added_coldata = added_coldata |>
    left_join(these_results, by = c('barcode', 'sample_id')) |>
    select(- all_of(extra_colnames))

if(any(is.na(added_coldata))) {
    stop("Some cell counts were NA after merging with SPE.")
}

colData(spe) = cbind(colData(spe), added_coldata)

################################################################################
#   Add up-to-date cell counts from cellpose and clarify existing ones
################################################################################

#   Rename outdated VistoSeg counts for clarity
spe$VistoSeg_count_deprecated = spe$counts
spe$counts = NULL

#   Add cellpose counts by adding all cell types for CART counts at each spot
spe$cellpose_count = colData(spe) |>
    as_tibble() |>
    select(starts_with('cart_')) |>
    rowSums()

#   Save a new SPE object
saveRDS(spe, spe_IF_out)

session_info()
