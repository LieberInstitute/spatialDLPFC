library(here)
library(tidyverse)
library(sessioninfo)

out_dir = here('processed-data', 'synapse_upload', '04-nda')
rna_seq_path = here('processed-data', 'synapse_upload', '04-nda', 'rna_seq.csv')

col_names <- c(
    "subjectkey", "src_subject_id", "interview_date", "interview_age", "sex",
    "race", 'phenotype', 'phenotype_description','twins_study', 'sibling_study',
    'family_study', 'sample_taken', 'sample_id_original', 'sample_description',
    'biorepository', 'patient_id_biorepository', 'sample_id_biorepository'
)

meta_df = read_csv(rna_seq_path, show_col_types = FALSE, skip = 1) |>
    group_by(src_subject_id) |>
    slice_head(n = 1) |>
    ungroup() |>
    mutate(
        race = "White",
        phenotype = "Control",
        phenotype_description = "Neurotypical control",
        twins_study = "No",
        sibling_study = "No",
        family_study = "No",
        sample_taken = "Yes",
        sample_id_original = src_subject_id,
        sample_description = "brain",
        biorepository = NA,
        patient_id_biorepository = NA,
        sample_id_biorepository = NA
    ) |>
    select(all_of(col_names))

out_path = file.path(out_dir, 'genomics_subject.csv')

#   Mimic the submission template from NDA, so this "CSV" can be directly
#   validated with the validator without any reformatting
write_csv(meta_df, out_path)
formatted_info = c(
    paste0('genomics_subject,02', paste(rep(',', ncol(meta_df) - 2), collapse = "")),
    readLines(out_path)
)
writeLines(formatted_info, out_path)

session_info()
