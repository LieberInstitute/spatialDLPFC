#  Create the 'manifest.tsv' file, the "main" file required for upload to
#  synapse

library("sessioninfo")
library("here")
library("readxl")
library("SpatialExperiment")
library("jaffelab")

template_path <- here(
    "code", "synapse_upload", "templates", "template_manifest.xlsx"
)

spe_path <- here(
    "processed-data", "rdata", "spe", "01_build_spe", "spe_filtered_final.Rdata"
)

fastq_naming_path <- here(
    "processed-data", "synapse_upload", "01-prepare_fastq",
    "fastq_renaming_scheme.csv"
)

write_path <- here(
    "processed-data", "synapse_upload", "02-metadata_files", "manifest.tsv"
)

###############################################################################
#  Load and preprocess phenotype data and metadata
###############################################################################

fastq_naming <- read.csv(fastq_naming_path)

#   Load phenotype data and take the first row for each sample ID
load(spe_path, verbose = TRUE)
pd <- colData(spe)
pd <- pd[match(unique(pd$sample_id), pd$sample_id), ]

#   Sample IDs have slight differences between spaceranger data and the SPE
#   object. We'll use the former as the "correct IDs"
pd$sample_id_correct <- sapply(
    pd$sample_id,
    function(id) {
        #   Find the 1 sample ID in the spaceranger table corresponding to this
        #   SPE sample ID
        index <- grep(id, unique(fastq_naming$sample_id))
        stopifnot(length(index) == 1)

        return(unique(fastq_naming$sample_id)[index])
    }
)

paths <- c(
    here("processed-data", "synapse_upload", "02-metadata_files", "assay.csv"),
    here(
        "processed-data", "synapse_upload", "02-metadata_files",
        "biospecimen.csv"
    ),
    here(
        "processed-data", "synapse_upload", "02-metadata_files",
        "individual.csv"
    ),
    fastq_naming$new_path
)

###############################################################################
#  Populate data frame that will become TSV
###############################################################################

num_metadata <- 3
num_fastq <- length(fastq_naming$new_path)

#  Populate a data frame
meta_df <- data.frame(
    "path" = paths,
    "parent" = c(
        rep("syn32383331", num_metadata), rep("syn32383329", num_fastq)
    ),
    "individualID" = c(
        rep(NA, num_metadata),
        pd$subject[
            match(fastq_naming$sample_id, pd$sample_id_correct)
        ]
    ),
    "specimenID" = c(rep(NA, num_metadata), fastq_naming$sample_id),
    "isMultiIndividual" = c(
        rep(NA, num_metadata), rep(FALSE, num_fastq)
    ),
    "isMultiSpecimen" = c(
        rep(NA, num_metadata), rep(FALSE, num_fastq)
    ),
    "assay" = NA, # TODO: need this!
    "libraryID" = NA,
    "fileFormat" = c(
        rep("csv", num_metadata), rep("fastq", num_fastq)
    ),
    "consortium" = "PEC",
    "study" = "LIBD_U01_DLPFC",
    "grant" = NA, # TODO: need this!
    "resourceType" = "experimentalData",
    "dataType" = NA,
    "dataSubtype" = c(
        rep("metadata", num_metadata), rep("raw", num_fastq)
    ),
    "metadataType" = c(
        "assay", "biospecimen", "individual",
        rep(NA, num_fastq)
    ),
    "analysisType" = NA
)

#   Ensure we included all the required columns
template_names <- colnames(read_excel(template_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.table(meta_df, write_path, row.names = FALSE, sep = "\t")

session_info()
