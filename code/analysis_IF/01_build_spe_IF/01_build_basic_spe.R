library("SpatialExperiment")
library("spatialLIBD")
library("here")
library("rtracklayer")
library("lobstr")
library("sessioninfo")

## Create output directories
dir_rdata <- here::here("processed-data", "rdata", "spe_IF", "01_build_spe_IF")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Define some info for the samples
sample_info <- data.frame(
    sample_id = c(
        "Br2720_Ant_IF",
        "Br6432_Ant_IF",
        "Br6522_Ant_IF",
        "Br8667_Post_IF"
    )
)
sample_info$subject <- sample_info$sample_id
sample_info$sample_path <-
    file.path(
        here::here("processed-data", "01_spaceranger_IF"),
        sample_info$sample_id,
        "outs"
    )
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
## https://github.com/LieberInstitute/spatialDLPFC/blob/main/raw-data/sample_info/Visium_IF_DLPFC_MasterExcel_01262022.xlsx
donor_info <- data.frame(
    subject = c("Br2720_Ant_IF", "Br6432_Ant_IF", "Br6522_Ant_IF", "Br8667_Post_IF"),
    age = c(48.22, 48.88, 33.38, 37.33),
    sex = c("F", "M", "M", "F"),
    race = "CAUC",
    pmi = c(25.5, 30.5, 31.5, 13.5),
    diagnosis = c("Control", "Control", "Control", "Control"),
    rin = c(7.5, 7.4, 7.7, 6.9)
)

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)


## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE
)
Sys.time()
# [1] "2021-12-01 14:07:32 EST"
# 2021-12-01 14:07:32 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2021-12-01 14:09:36 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2021-12-01 14:09:42 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2021-12-01 14:09:45 rtracklayer::import: reading the reference GTF file
# 2021-12-01 14:10:44 adding gene information to the SPE object
# 2021-12-01 14:10:44 adding information used by spatialLIBD
# [1] "2021-12-01 14:10:52 EST"

Sys.time()

# ## Add segmentation images created by Madhavi Tippani
# image_types <- c("Abeta", "Abeta_seg", "pTau", "pTau_seg", "DAPI", "DAPI_seg", "merge", "merge_seg")
#
# for (img in image_types) {
#   for (res in c("lowres")) {
#     spe <- add_images(
#       spe = spe,
#       image_dir = here::here("raw-data", "Images", "VisiumIF", "VistoSeg"),
#       image_pattern = paste0(img, "_", res),
#       image_id_current = res
#     )
#   }
# }

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

## Read in cell counts and segmentation results
segmentations_list <-
    lapply(sample_info$sample_id, function(sampleid) {
        file <-
            here(
                "processed-data",
                "01_spaceranger_IF",
                sampleid,
                "outs",
                "spatial",
                "tissue_spot_counts.csv"
            )
        if (!file.exists(file)) {
            return(NULL)
        }
        x <- read.csv(file)
        x$key <- paste0(x$barcode, "_", sampleid)
        return(x)
    })
## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <-
    Reduce(function(...) {
        merge(..., all = TRUE)
    }, segmentations_list[lengths(segmentations_list) > 0])

## Add the information
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <-
    segmentations[segmentation_match, -which(
        colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
    )]
colData(spe) <- cbind(colData(spe), segmentation_info)

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 8748
length(no_expr) / nrow(spe) * 100
# [1] 23.90099
spe <- spe[-no_expr, ]

## For visualizing this later with spatialLIBD
spe$overlaps_tissue <-
    factor(ifelse(spe$in_tissue, "in", "out"))

## Save with and without dropping spots outside of the tissue
spe_raw <- spe
saveRDS(spe_raw, file.path(dir_rdata, "spe_raw.rds"))

saveRDS(spe_raw, file.path(dir_rdata, "spe_raw.rds"))

## Size in Gb
lobstr::obj_size(spe_raw)
# 1.651702

## Now drop the spots outside the tissue
spe <- spe_raw[, spe_raw$in_tissue]
dim(spe)
# [1] 27853 38287

## Now drop the spots outside the tissue
spe <- spe_raw[, spe_raw$in_tissue]
dim(spe)
# [1] 27853 38287
## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
    message("removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}

lobstr::obj_size(spe)
# 1.534376

saveRDS(spe, file.path(dir_rdata, "spe.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
