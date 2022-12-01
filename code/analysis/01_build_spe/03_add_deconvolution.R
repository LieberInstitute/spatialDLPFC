library("spatialLIBD")
library("here")
library("sessioninfo")


## Takes about 4 min
Sys.time()
load(
    here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters.Rdata"
    ),
    verbose = TRUE
)
Sys.time()

## Check the size in GB
Sys.time()
lobstr::obj_size(spe) ## takes about 2 min
Sys.time()
# 6.92 GB

## Import spot deconvolution results
for (resolution in c("broad", "layer")) {
    for (deconvo in c("01-tangram", "03-cell2location", "04-spotlight")) {
        message(Sys.time(), " importing: ", resolution, " for ", deconvo)
        spe <- cluster_import(
            spe,
            here(
                "processed-data",
                "spot_deconvo",
                deconvo,
                "nonIF",
                resolution,
                "raw_results"
            ),
            prefix = paste0(resolution, "_", gsub(".+-", "", deconvo), "_")
        )
    }
}

## Fix some variables
vars <- colnames(colData(spe))
colnames(colData(spe))[grep("^X10x", vars)] <-
    gsub("X10x", "SpaceRanger_10x", vars[grep("^X10x", vars)])
colnames(colData(spe))[grep("^bayes", vars)] <-
    gsub("bayes", "Bayes", vars[grep("^bayes", vars)])
colnames(colData(spe))[grep("^region$", vars)] <-
    gsub("region", "position", vars[grep("^region$", vars)])

## Check the size in GB
lobstr::obj_size(spe)
# 6.96 GB

## Delete original VistoSeg counts (prior to updates in VistoSeg)
spe$VistoSeg_count_deprecated <- spe$count
spe$count <- NULL

## Import updated VistoSeg cell counts
nonIF_id_path <-
    here("processed-data", "spot_deconvo", "nonIF_ID_table.csv")
nonIF_counts_path <- here(
    "processed-data",
    "rerun_spaceranger",
    "{sample_id}",
    "outs",
    "spatial",
    "tissue_spot_counts.csv"
)

## Import code from
## https://github.com/LieberInstitute/spatialDLPFC/blob/444ce5cb408b90aa59066379df7a0e30cf0d2447/code/spot_deconvo/05-shared_utilities/01-r_to_python.R#L130-L167
id_table <- read.csv(nonIF_id_path)

## Initialize the new variables
spe$VistoSeg_count <- NA
spe$VistoSeg_percent <- NA

for (sample_id in unique(spe$sample_id)) {
    message(Sys.time(), " processing sample id ", sample_id)
    #   Correctly determine the path for the cell counts for this sample, then
    #   read in
    long_id <-
        id_table[match(sample_id, id_table$short_id), "long_id"]
    this_path <-
        sub("{sample_id}", long_id, nonIF_counts_path, fixed = TRUE)
    cell_counts <- read.csv(this_path)

    #   All spots in the object should have counts
    stopifnot(all(colnames(spe[, spe$sample_id == sample_id]) %in%
        cell_counts$barcode))

    #   Line up the rows of 'cell_counts' with the sample-subsetted SPE object
    cell_counts <- cell_counts[match(
        colnames(spe[, spe$sample_id == sample_id]),
        cell_counts$barcode
    ), ]

    #   Add this sample's counts to the SPE object
    spe$VistoSeg_count[spe$sample_id == sample_id] <-
        cell_counts$Nmask_dark_blue

    ## Also add the percent of the spot covered, which can be useful
    ## for detecting neuropil spots
    spe$VistoSeg_percent[spe$sample_id == sample_id] <-
        cell_counts$Pmask_dark_blue
}

#   Ensure counts were read in for all spots in the object
if (any(is.na(spe$VistoSeg_count))) {
    stop("Did not find cell counts for all non-IF spots.")
}
stopifnot(all(!is.na(spe$VistoSeg_percent)))


## Check the size in GB
Sys.time()
lobstr::obj_size(spe) ## takes about 2 min
Sys.time()
# 6.96 GB

## Save for later use, takes about 11 min
Sys.time()
saveRDS(
    spe,
    file = here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
    )
)
Sys.time()


## Start removing pieces for the shiny apps
imgData(spe) <- imgData(spe)[imgData(spe)$image_id == "lowres", ]
## Check the size in GB
lobstr::obj_size(spe)
# 4.59 GB

## Drop the counts which take quite a bit of space
counts(spe) <- NULL
## Check the size in GB
lobstr::obj_size(spe)
# 2.45 GB


## Save for use in the spatialLIBD shiny apps, takes about 3 min
Sys.time()
saveRDS(
    spe,
    file = here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_subset_for_spatialLIBD.rds"
    )
)
Sys.time()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
