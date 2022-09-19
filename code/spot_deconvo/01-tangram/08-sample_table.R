#   Spot-deconvolution code for this project processes many samples, but a
#   number of different naming conventions are used for different files/
#   objects in the project (e.g. "Br2720_ant", "DLPFC_Br2720_ant_2", and
#   "DLPFC_Br2720_ant_manual_alignment" all refer to the same sample). This
#   script generates a table mapping between two sets of ID conventions that
#   will be needed for the spot deconvolution code.

library("here")
library("SpatialExperiment")

spe_in <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
)

spaceranger_dir <- here("processed-data", "rerun_spaceranger")
out_path <- here("processed-data", "spot_deconvo", "nonIF_ID_table.csv")

load(spe_in)

#   Read in both sets of ID naming conventions
old_spaceranger_ids <- list.files(spaceranger_dir)
spe_ids <- unique(spe$sample_id)

#   Order the second set of naming conventions along the first
new_spaceranger_ids <- c()
for (spe_id in spe_ids) {
    index <- grep(spe_id, old_spaceranger_ids)
    stopifnot(length(index) == 1)

    new_spaceranger_ids <- c(new_spaceranger_ids, old_spaceranger_ids[index])
}

#   Write a table mapping both sets of conventions
a <- data.frame("short_id" = spe_ids, "long_id" = new_spaceranger_ids)
write.csv(a, out_path, row.names = FALSE, quote = FALSE)
