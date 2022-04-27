library("SpatialExperiment")
library("spatialLIBD")
library("here")
library("rtracklayer")
library("lobstr")
library("sessioninfo")

## Create output directories
dir_rdata <- here::here("processed-data", "spe_IF","01_build_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Define some info for the samples
sample_info <- data.frame(
  sample_id = c(
    "Br8667_Post_IF",
    "Br2720_Ant_IF",
    "Br6522_Ant_IF",
    "Br6432_Ant_IF",
  )
)
sample_info$sample_id <- gsub(".*_", "", sample_info$sample_id)
sample_info$sample_path <-
  file.path(
    here::here("processed-data", "spaceranger"),
    sample_info$sample_id,
    "outs"
  )
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
## https://github.com/LieberInstitute/spatialDLPFC/blob/main/raw-data/sample_info/Visium_IF_DLPFC_MasterExcel_01262022.xlsx
donor_info <- data.frame(
  subject = c("Br8667_Post", "Br2720_Ant", "Br6522_Ant", "Br6432_Ant"),
  age = c(37.33, 48.22, 33.38, 48.88),
  sex = c("F", "F", "M", "M"),
  race = "CAUC",
  pmi = c(13.5, 25.5, 31.5, 30.5),
  diagnosis = c("Control", "Control", "Control", "Control"),
  rin = c(6.9, 7.5, 7.7, 7.4)
  # BCrating = c("Def AD", "Def AD", "Prob AD", "No AP"),
  # braak = c("B3", "B3", "B3", "B2"),
  # cerad = c("C3", "C3", "C3", "C0")
)

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)


## Build basic SPE
Sys.time()
spe_wholegenome <- read10xVisiumWrapper(
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
spe_targeted <- read10xVisiumWrapper(
  gsub(
    "spaceranger",
    "spaceranger_targeted",
    sample_info$sample_path
  ),
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres"),
  load = TRUE
)
Sys.time()
# [1] "2021-12-01 14:11:47 EST"
# 2021-12-01 14:11:48 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2021-12-01 14:12:17 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2021-12-01 14:12:24 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2021-12-01 14:12:26 rtracklayer::import: reading the reference GTF file
# 2021-12-01 14:13:17 adding gene information to the SPE object
# 2021-12-01 14:13:18 adding information used by spatialLIBD
# [1] "2021-12-01 14:13:19 EST"

# ## Add images created by Madhavi Tippani
# image_types <- c("Abeta", "Abeta_seg", "pTau", "pTau_seg", "DAPI", "DAPI_seg", "merge", "merge_seg")
# 
# for (img in image_types) {
#   for (res in c("lowres")) {
#     spe_wholegenome <- add_images(
#       spe = spe_wholegenome,
#       image_dir = here("processed-data", "Images", "spatialLIBD_images"),
#       image_pattern = paste0(img, "_", res),
#       image_id_current = res
#     )
#   }
# }

## Update the imData in the targeted sequencing (do we need this part?)
imgData(spe_targeted) <- imgData(spe_wholegenome)

## This is the case since we didn't use the --target-panel option when
## running spaceranger as described at
## https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count
stopifnot(identical(rowData(spe_wholegenome), rowData(spe_targeted)))
# stopifnot(identical(colData(spe_wholegenome), colData(spe_targeted)))
## The above is no longer true since we are reading in the clustering results
## from SpaceRanger which are different between the regular Visium and the
## targeted sequencing Visium.

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
spe_wholegenome <- add_design(spe_wholegenome)
spe_targeted <- add_design(spe_targeted) #Do we need this line?

## Read in cell counts and segmentation results
segmentations_list <-
  lapply(sample_info$sample_id, function(sampleid) {
    file <-
      here(
        "processed-data",
        "spaceranger",
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
segmentation_match <- match(spe_wholegenome$key, segmentations$key)
segmentation_info <-
  segmentations[segmentation_match, -which(
    colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
  )]
colData(spe_wholegenome) <- cbind(colData(spe_wholegenome), segmentation_info)
colData(spe_targeted) <- cbind(colData(spe_targeted), segmentation_info) #Maybe we do not need this line?!

## Remove genes with no data
no_expr <- which(rowSums(counts(spe_wholegenome)) == 0)
length(no_expr)
# [1] 8748
length(no_expr) / nrow(spe_wholegenome) * 100
# [1] 23.90099
spe_wholegenome <- spe_wholegenome[-no_expr, ]

no_expr <- which(rowSums(counts(spe_targeted)) == 0)
length(no_expr)
# [1] 13116
length(no_expr) / nrow(spe_targeted) * 100
# [1] 35.83509
spe_targeted <- spe_targeted[-no_expr, ]

## For visualizing this later with spatialLIBD (maybe do need this either because it is for targeted sequencing?)
spe_targeted$overlaps_tissue <-
  spe_wholegenome$overlaps_tissue <-
  factor(ifelse(spe_wholegenome$in_tissue, "in", "out"))

## Save with and without dropping spots outside of the tissue
spe_raw_wholegenome <- spe_wholegenome
spe_raw_targeted <- spe_targeted

saveRDS(spe_raw_wholegenome, file.path(dir_rdata, "spe_raw_wholegenome.rds"))
#saveRDS(spe_raw_targeted, file = file.path(dir_rdata, "spe_raw_targeted.rds"))

## Size in Gb
lobstr::obj_size(spe_raw_wholegenome) / 1024^3
# 1.651702
$lobstr::obj_size(spe_raw_targeted) / 1024^3
# 1.235324

## Now drop the spots outside the tissue
spe_wholegenome <- spe_raw_wholegenome[, spe_raw_wholegenome$in_tissue]
dim(spe_wholegenome)
# [1] 27853 38287
## Remove spots without counts
if (any(colSums(counts(spe_wholegenome)) == 0)) {
  message("removing spots without counts for spe_wholegenome")
  spe_wholegenome <- spe_wholegenome[, -which(colSums(counts(spe_wholegenome)) == 0)]
  dim(spe_wholegenome)
}

# spe_targeted <- spe_raw_targeted[, spe_raw_targeted$in_tissue]
# dim(spe_targeted)
# [1] 23485 38287

## Remove spots without counts
if (any(colSums(counts(spe_targeted)) == 0)) {
  message("removing spots without counts for spe_targeted")
  spe_targeted <-
    spe_targeted[, -which(colSums(counts(spe_targeted)) == 0)]
  dim(spe_targeted)
}

lobstr::obj_size(spe_wholegenome) / 1024^3
# 1.534376
#lobstr::obj_size(spe_targeted) / 1024^3
# 1.198796

saveRDS(spe_wholegenome, file.path(dir_rdata, "spe_wholegenome.rds"))
#saveRDS(spe_targeted, file = file.path(dir_rdata, "spe_targeted.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()