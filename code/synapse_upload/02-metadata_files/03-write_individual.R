#   Create the 'individual.csv' file, 1 of 3 required metadata files for upload
#   to synapse

library("sessioninfo")
library("here")
library("readxl")
library("SpatialExperiment")
library("tidyverse")

template_path <- here(
    "code", "synapse_upload", "templates", "template_individual_human.xlsx"
)

spe_path <- here(
    "processed-data", "rdata", "spe", "01_build_spe", "spe_filtered_final.Rdata"
)

write_path_synapse <- here(
    "processed-data", "synapse_upload", "02-metadata_files", "individual.csv"
)

write_path_paper <- here(
    "processed-data", "synapse_upload", "02-metadata_files", "demographics.csv"
)

###############################################################################
#  Load and preprocess phenotype data
###############################################################################

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

write.csv(meta_df, write_path_synapse, row.names = FALSE)

#   Take the non-NA columns with their original colnames and write to CSV
#   (for use as a supplementary table in the paper). Add the brain region
#   (DLPFC for all donors)
pd |>
    as_tibble() |>
    select(c("subject", "sex", "diagnosis", "age")) |>
    mutate("brain_region" = "DLPFC") |>
    write.csv(file = write_path_paper, quote = FALSE, row.names = FALSE)

session_info()
