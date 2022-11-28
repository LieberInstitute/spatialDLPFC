#   For the Visium-IF data, add cell-type counts to the SpatialExperiment
#   object and save a new copy. This includes counts from each deconvolution
#   tool as well as the "ground-truth" counts from cellpose + the trained
#   DecisionTreeClassifier.

library("here")
library("jaffelab")
library("spatialLIBD")
library("sessioninfo")

cell_groups <- c("broad", "layer")

sample_ids <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "sample_ids.txt"
)

deconvo_tools <- c("01-tangram", "03-cell2location", "04-spotlight")

spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)

spe_IF_out <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF", "spe.rds"
)

#   "Ground-truth" cell counts from cellpose + trained classification tree
actual_paths <- here(
    "processed-data", "spot_deconvo", "02-cellpose", "{sample_id}", "clusters.csv"
)

#   Cell counts estimated by each deconvolution tool
observed_paths <- here(
    "processed-data", "spot_deconvo", "{deconvo_tool}", "IF", "{cell_group}",
    "{sample_id}", "clusters.csv"
)

spe <- readRDS(spe_IF_in)

#   For each sample, resolution, and deconvolution tool (and ground-truth CART
#   counts), read in cell-type counts and add each as a column to colData(spe)
for (sample_id in unique(spe$sample_id)) {
    #   Read in the "ground-truth" counts for this sample
    actual_path <- sub("\\{sample_id\\}", sample_id, actual_paths)
    a <- read.csv(actual_path)
    stopifnot(all(a$key %in% spe$key))

    #   Use informative column names, and don't duplicate the 'key' column
    actual_key <- a$key
    a$counts <- a$n_cells
    a$n_cells <- NULL
    a$key <- NULL
    colnames(a) <- paste0("CART_", colnames(a))

    #   Add CART counts to the object
    for (cname in colnames(a)) {
        if (sample_id == unique(spe$sample_id)[1]) {
            spe[[cname]] <- NA
        }

        spe[[cname]][match(actual_key, spe$key)] <- a[[cname]]
    }

    for (cell_group in cell_groups) {
        for (deconvo_tool in deconvo_tools) {
            short_deconvo_tool <- ss(deconvo_tool, "-", 2)

            #   Fill in the path with the actual values instead of the
            #   placeholders
            observed_path <- sub("\\{sample_id\\}", sample_id, observed_paths)
            observed_path <- sub("\\{deconvo_tool\\}", deconvo_tool, observed_path)
            observed_path <- sub("\\{cell_group\\}", cell_group, observed_path)

            #   Read in cell-type counts for this deconvo tool
            a <- read.csv(observed_path)
            stopifnot(all(a$key %in% spe$key))

            #   Use unique column names to add to the SPE object
            observed_key <- a$key
            a$key <- NULL
            colnames(a) <- paste0(
                short_deconvo_tool, "_", cell_group, "_",
                gsub("\\.", "_", colnames(a))
            )

            #   Add counts to the object
            for (cname in colnames(a)) {
                if (sample_id == unique(spe$sample_id)[1]) {
                    spe[[cname]] <- NA
                }

                spe[[cname]][match(observed_key, spe$key)] <- a[[cname]]
            }
        }
    }
}

saveRDS(spe, spe_IF_out)

session_info()
