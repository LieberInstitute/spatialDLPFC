## This script requires R 4.1
# module load conda_R/4.1.x

## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")

## vis
library("spatialLIBD")


## Define some info for the samples
sample_info <- data.frame(
  sample_id = c(
    "DLPFC_Br2743_ant_manual_alignment",
    "DLPFC_Br2743_mid_manual_alignment_extra_reads",
    "DLPFC_Br2743_post_manual_alignment",
    "DLPFC_Br3942_ant_manual_alignment",
    "DLPFC_Br3942_mid_manual_alignment",
    "DLPFC_Br3942_post_manual_alignment",
    "DLPFC_Br6423_ant_manual_alignment_extra_reads",
    "DLPFC_Br6423_mid_manual_alignment",
    "DLPFC_Br6423_post_extra_reads",
    "DLPFC_Br8492_ant_manual_alignment",
    "DLPFC_Br8492_mid_manual_alignment_extra_reads",
    "DLPFC_Br8492_post_manual_alignment",
    "Round2/DLPFC_Br2720_ant_manual_alignment",
    "Round2/DLPFC_Br2720_mid_manual_alignment",
    "Round2/DLPFC_Br2720_post_extra_reads",
    "Round2/DLPFC_Br6432_ant_manual_alignment",
    "Round2/DLPFC_Br6432_mid_manual_alignment",
    "Round2/DLPFC_Br6432_post_manual_alignment",
    "Round3/DLPFC_Br6471_ant_manual_alignment_all",
    "Round3/DLPFC_Br6471_mid_manual_alignment_all",
    "Round3/DLPFC_Br6471_post_manual_alignment_all",
    "Round3/DLPFC_Br6522_ant_manual_alignment_all",
    "Round3/DLPFC_Br6522_mid_manual_alignment_all",
    "Round3/DLPFC_Br6522_post_manual_alignment_all",
    "Round3/DLPFC_Br8325_ant_manual_alignment_all",
    "Round3/DLPFC_Br8325_mid_manual_alignment_all",
    "Round3/DLPFC_Br8325_post_manual_alignment_all",
    "Round3/DLPFC_Br8667_ant_extra_reads",
    "Round3/DLPFC_Br8667_mid_manual_alignment_all",
    "Round3/DLPFC_Br8667_post_manual_alignment_all",
    "Round4/DLPFC_Br2720_ant_2",
    "Round4/DLPFC_Br6432_ant_2",
    "Round4/DLPFC_Br8325_mid_2"
  ),
  subjects = c(rep(
    c("Br2743", "Br3942", "Br6423", "Br8492", "Br2720","Br6432","Br6471","Br6522","Br8325","Br8667"),
    each = 3
  ),"Br2720","Br6432","Br8325"),
  regions = c(rep(
    c("anterior", "middle", "posterior"),
    10
  ),"anterior","anterior","middle"),
  sex = c(rep(
    c("M", "M", "M", "F","F","M","M","M","F","F"),
    each = 3
  ),"F","M","F"),
  age = c(rep(
    c(61.54, 47.53, 51.73, 53.40,48.22,48.88,55.46,33.39,57.62,37.33),
    each = 3
  ),48.22,48.88,57.62),
  diagnosis = "Control"
)
sample_info$sample_path <- file.path(
  here::here("processed-data", "NextSeq"),
  sample_info$sample_id,
  "outs"
)
stopifnot(all(file.exists(sample_info$sample_path)))

#clean up sample_id
sample_info$sample_id <- gsub("_all|_extra_reads|DLPFC_|_manual_alignment", "", basename(sample_info$sample_id))

## Build SPE object
Sys.time()
spe_wrapper <- read10xVisiumWrapper(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  data = "raw",
  images = c("lowres", "hires", "detected", "aligned"),
  load = TRUE,
  reference_gtf = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
)
Sys.time()

# 2021-11-23 15:05:05 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2021-11-23 15:12:23 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2021-11-23 15:13:09 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2021-11-23 15:13:15 rtracklayer::import: reading the reference GTF file
# 2021-11-23 15:14:03 adding gene information to the SPE object
# 2021-11-23 15:14:04 adding information used by spatialLIBD

## Add the experimental information
### Equivalent to:
## colData(spe_wrapper)$subject <- ...
spe_wrapper$subject <- sample_info$subjects[match(spe_wrapper$sample_id, sample_info$sample_id)]
spe_wrapper$region <- sample_info$regions[match(spe_wrapper$sample_id, sample_info$sample_id)]
spe_wrapper$sex <- sample_info$sex[match(spe_wrapper$sample_id, sample_info$sample_id)]
spe_wrapper$age <- sample_info$age[match(spe_wrapper$sample_id, sample_info$sample_id)]
spe_wrapper$diagnosis <- sample_info$diagnosis[match(spe_wrapper$sample_id, sample_info$sample_id)]

## Remove genes with no data
no_expr <- which(rowSums(counts(spe_wrapper)) == 0)
length(no_expr)
# [1] 7468
length(no_expr) / nrow(spe_wrapper) * 100
# [1] 20.40381

#remove genes that are not expressed in any spots
spe_wrapper <- spe_wrapper[-no_expr, ]

spe_wrapper <- spe_wrapper[, which(spatialData(spe_wrapper)$in_tissue=="TRUE")]

## Remove spots without counts
spe_wrapper <- spe_wrapper[, -which(colSums(counts(spe_wrapper)) == 0)]

## Read in cell counts and segmentation results
segmentations_list <- lapply(sample_info$sample_id, function(sampleid) {
  current<-sample_info$sample_path[sample_info$sample_id==sampleid]
  file <- file.path(current, "spatial", "tissue_spot_counts.csv")
  if(!file.exists(file)) return(NULL)
  x <- read.csv(file)
  x$key <- paste0(x$barcode, "_", sampleid)
  return(x)
})
## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <- Reduce(function(...) merge(..., all = TRUE), segmentations_list[lengths(segmentations_list) > 0])

## Add the information
segmentation_match <- match(spe_wrapper$key, segmentations$key)
segmentation_info <- segmentations[segmentation_match, - which(colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")), drop=FALSE]
colData(spe_wrapper) <- cbind(colData(spe_wrapper), segmentation_info)




