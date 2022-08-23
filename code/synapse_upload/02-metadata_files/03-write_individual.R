#   Create the 'individual.csv' file, 1 of 3 required metadata files for upload
#   to synapse

library("sessioninfo")
library("here")
library("readxl")
library("SpatialExperiment")

template_path <- here(
    "code", "synapse_upload", "templates", "template_individual_human.xlsx"
)

spe_path <- here(
    "processed-data", "rdata", "spe", "01_build_spe", "spe_filtered_final.Rdata"
)

fastq_naming_path <- here(
    "processed-data", "synapse_upload", "01-prepare_fastq",
    "fastq_renaming_scheme.csv"
)

write_path <- here(
    "processed-data", "synapse_upload", "02-metadata_files", "individual.csv"
)

###############################################################################
#  Load and preprocess phenotype data
###############################################################################

fastq_naming <- read.csv(fastq_naming_path)

#   Load phenotype data and take the first row for each donor
load(spe_path, verbose = TRUE)
pd <- colData(spe)
pd <- pd[match(unique(pd$subject), pd$subject), ]

#   Tweak phenotype data to fit Synapse conventions
pd$sex[pd$sex == "M"] <- "male"
pd$sex[pd$sex == "F"] <- "female"
pd$diagnosis[pd$diagnosis == "Control"] <- "control"

###############################################################################
#  Populate "individual" data frame
###############################################################################

meta_df <- data.frame(
    "individualID" = pd$subject,
    "individualIdSource" = "LIBD",
    "species" = "Human",
    "reportedGender" = pd$sex,
    "sexChromosome" = NA,
    "race" = NA,
    "ethnicity" = NA,
    "genotypeInferredAncestry" = NA,
    "familialRelationship" = NA,
    "IQ" = NA,
    "BMI" = NA,
    "primaryDiagnosis" = pd$diagnosis,
    "primaryDiagnosisDetail" = NA,
    "otherDiagnosis" = NA,
    "otherDiagnosisDetail" = NA,
    "familyHistory" = NA,
    "familyHistoryDetails" = NA,
    "ageOnset" = NA,
    "neuropathDescription" = NA,
    "dementia" = NA,
    "CDR" = NA,
    "Braak" = NA,
    "otherMedicalDx" = NA,
    "otherMedicalDetail" = NA,
    "ageDeath" = pd$age,
    "ageDeathUnits" = "years",
    "causeDeath" = NA,
    "mannerDeath" = NA,
    "postmortemTox" = NA,
    "postmortemToxDetails" = NA,
    "postmortemToxSource" = NA,
    "medRecordTox" = NA,
    "PMICertain" = NA,
    "PMI" = NA,
    "pH" = NA
)

#  Ensure we included all the required columns
template_names <- colnames(read_excel(template_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.csv(meta_df, write_path, row.names = FALSE)

session_info()
