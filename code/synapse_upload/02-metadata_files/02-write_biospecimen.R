#   Create the 'biospecimen.csv' file, 1 of 3 required metadata files for upload
#   to synapse

library("sessioninfo")
library("here")
library("readxl")
library("SpatialExperiment")

template_path <- here(
    "code", "synapse_upload", "templates", "template_biospecimen.xlsx"
)

spe_path <- here(
    "processed-data", "rdata", "spe", "01_build_spe", "spe_filtered_final.Rdata"
)

fastq_naming_path <- here(
    "processed-data", "synapse_upload", "01-prepare_fastq",
    "fastq_renaming_scheme.csv"
)

write_path <- here(
    "processed-data", "synapse_upload", "02-metadata_files", "biospecimen.csv"
)

###############################################################################
#  Load and preprocess phenotype data
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

###############################################################################
#  Populate biospecimen data frame
###############################################################################

meta_df <- data.frame(
    "individualID" = pd$subject,
    "specimenID" = pd$sample_id_correct,
    "organ" = "brain",
    "organWeight" = NA,
    "organRIN" = NA,
    "tissue" = "dorsolateral prefrontal cortex",
    "isPostMortem" = TRUE,
    "BrodmannArea" = NA,
    "nucleicAcidSource" = NA,
    "cellType" = NA,
    "reprogrammedCellType" = NA,
    "terminalDifferentiationPoint" = NA,
    "passage" = NA,
    "samplingAge" = pd$age,
    "samplingAgeUnits" = "years",
    "sampleStatus" = NA
)

#  Ensure we included all the required columns
template_names <- colnames(read_excel(template_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.csv(meta_df, write_path, row.names = FALSE)

session_info()
