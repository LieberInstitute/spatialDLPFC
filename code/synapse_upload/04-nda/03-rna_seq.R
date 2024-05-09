#   Prepare "data expected" CSV for all RNA-seq data for this study, which
#   consists of 3 components:
#       - n = 30 Visium FASTQs
#       - n = 4 Visium-SPG FASTQs
#       - n = 19 snRNA-seq FASTQs

library(tidyverse)
library(here)
library(readxl)
library(sessioninfo)

he_sample_info_path <- here(
    "raw-data", "sample_info", "Visium_dlpfc_mastersheet.xlsx"
)
fastq_mapping_path <- here(
    "raw-data", "FASTQ_Globus", "fastq_mapping.csv"
)
pd_path <- here("processed-data", "synapse_upload", "04-nda", "pheno_data.csv")
out_dir = here('processed-data', 'synapse_upload', '04-nda')

#   All columns in the data structure (even those we have no info for)
col_names <- c(
    "subjectkey", "src_subject_id", "interview_date", "interview_age", "sex",
    "experiment_id", "cellid", "samplesubtype", "libraryid", "gen_software",
    "softwareversion", "referenceset", "otherreferenceset", "librarybatch",
    "sequencingbatch", "libraryselection", "libraryconstructionprotocol",
    "otherlibraryconstructprotocol", "librarysource", "otherlibrarysource",
    "readlength", "librarylayout", "totalreads", "numbercells",
    "readstrandorigin", "isstranded", "libraryversion", "validbarcodereads",
    "mediangenes", "medianumis", "gen_rin", "rnabatch", "ribozero_batch",
    "data_file1", "data_file1_type", "hcdi_tissue",
    "dlpfc_rna_isola_prepoperator", "flowcell_batch", "flowcell_lane_a",
    "flowcell_lane_b", "flowcell_name", "hemisphere", "rat280",
    "sample_id_biorepository", "psych_enc_exclude_reason", "flowcell_2",
    "flowcell_given_to_core", "flowcell_id", "flowcell_name_2", "study",
    "brodmann_area", "psych_enc_exclude", "ercc_added", "librarykit",
    "librarytype", "mappedreads_multimapped", "mappedreads_primary",
    "nucleicacidsource", "readlength_max", "readlength_min", "rnaseqid",
    "rrnarate", "samplebarcode", "sequencingplatform", "tissuestate",
    "celltype",
    "externalreference",
    "filename",
    "library_prep_batch",
    "platform",
    "assay",
    "hcdi_organ",
    "ethnicity",
    "psych_enc_datatype",
    "rna_type",
    "sequencing_assay",
    "submission_file_name",
    "file_status",
    "data_file5_type",
    "data_file5",
    "visium_protocol_version"
)

################################################################################
#   Gather basic sample info: one row per FASTQ
################################################################################

sample_info <- read_csv(fastq_mapping_path, show_col_types = FALSE) |>
    select(sample_id, donor, fastq_globus, assay) |>
    rename(
        data_file1 = fastq_globus,
        src_subject_id = sample_id,
        assay_internal = assay
    ) |>
    mutate(data_file1_type = "FASTQ") |>
    left_join(
        read_csv(pd_path, show_col_types = FALSE),
        by = "donor"
    ) |>
    select(-donor)

################################################################################
#   Add additional fields as appropriate for the assay: Visium, Visium-SPG, or
#   snRNA-seq
################################################################################

sample_info |>
    mutate(
        experiment_id = "LIBD spatial DLPFC",
        samplesubtype = ifelse(assay_internal == "snRNA-seq", 2, 3),
        referenceset = 2, # GrCh38
        libraryconstructionprotocol = 31, # Visium, which we added
        librarysource = ifelse(assay_internal == "snRNA-seq", 2, 1),
        librarylayout = 2, # paired-end
        hcdi_tissue = 16, # DLPFC
        psych_enc_exclude_reason = "Not excluded",
        brodmann_area = 46,
        psych_enc_exclude = 0, # Not excluded from study
        rnaseqid = src_subject_id,
        filename = data_file1,
        platform = case_when(
            assay_internal == "snRNA-seq" ~ "Illumina Novaseq 6000",
            assay_internal == "Visium" ~ "10x Genomics Visium",
            assay_internal == "Visium-SPG" ~ "10x Genomics Visium-SPG"
        ),
        assay = ifelse(assay_internal == "snRNA-seq", 13, 5),
        hcdi_organ = 6, # brain
        submission_file_name = data_file1,
        file_status = 1, # raw data
        visium_protocol_version = ifelse(assay_internal == "snRNA-seq", NA, "V1")
    ) |>
    select(any_of(col_names)) |>
    write_csv(file.path(out_dir, 'rna_seq.csv'))

session_info()
