#   Gather sample info to form two tibbles:
#
#   - One with 34 rows (all H&E and IF samples) and columns for sample ID,
#     GUID, donor, sex, age, interveiw_date, image file, alignment file, fastq
#     files, spaceranger JSON, and loupe version
#   - One with 10 rows (one per donor) with sample ID, GUID, donor, sex, age,
#     and interview_date
#
#   This info will be borrowed by other R scripts in this directory

library(here)
library(tidyverse)
library(sessioninfo)

he_sample_info_path = here(
    "processed-data", "rdata", "spe", "01_build_spe", "spe_sample_info.csv"
)
he_sr_params_path = here(
    'code', 'spaceranger', 'spaceranger_parameters.txt'
)
if_sr_params_path = here(
    'code', 'spaceranger_IF', 'spaceranger_IF_parameters.txt'
)
out_dir = here('processed-data', 'synapse_upload', '04-nda')

guids = c(
    "NDAR_INVBM879VT6", "NDAR_INVZN585KGQ", "NDAR_INVEJ479WEA",
    "NDAR_INVCH371AU1", "NDAR_INVGP756DP1", "NDAR_INVCE963WP1",
    "NDAR_INVPY786FH2", "NDAR_INVRD044HHV", "NDAR_INVUM090ZYX",
    "NDAR_INVTG948WW3"
)

################################################################################
#   H&E samples
################################################################################

he_info = read_csv(he_sample_info_path, show_col_types = FALSE) |>
    rename(src_subject_id = subjects) |>
    mutate(
        interview_age = as.integer(round(age * 12, 0)),
        sample_id_old = str_replace(sample_id, '_2$', '')
    ) |>
    select(sample_id_old, src_subject_id, sex, interview_age)

all_json_paths = file.path(
    list.files(
        here("processed-data", "rerun_spaceranger"), full.names = TRUE,
        pattern = '^DLPFC'
    ),
    "outs", "spatial", "scalefactors_json.json"
)
stopifnot(all(file.exists(all_json_paths)))

he_info = read_delim(
        he_sr_params_path, col_names = FALSE, delim = " ",
        show_col_types = FALSE
    ) |>
    mutate(
        sample_id_old = str_extract(X1, 'Br[0-9]{4}_(ant|mid|post)'),
        sample_id = sprintf('%s_%s', X2, X3),
        loupe_version = "5.1.0",
        fastq_files = str_replace_all(X6, ',', ';'),
        #   Date of Visium sequencing
        interview_date = '06/25/2020',
        stain = "H&E",
        spaceranger_json = sapply(
            sample_id_old, function(x) all_json_paths[grep(x, all_json_paths)]
        )
    ) |>
    rename(image_file = X4, alignment_file = X5) |>
    select(
        sample_id_old, sample_id, loupe_version, fastq_files, image_file,
        alignment_file, interview_date, stain, spaceranger_json
    ) |>
    #   Add phenotype data
    left_join(he_info, by = 'sample_id_old') |>
    select(-sample_id_old) |>
    #   Add GUIDs
    left_join(
        tibble(
            src_subject_id = unique(he_info$src_subject_id),
            subjectkey = guids
        ),
        by = 'src_subject_id'
    )

################################################################################
#   Phenotype data that applies to many NDA data structures
################################################################################

#   Just one row per donor
pd = he_info |>
    group_by(src_subject_id) |>
    slice_head(n = 1) |>
    ungroup() |>
    select(src_subject_id, subjectkey, interview_age, interview_date, sex)

write_csv(pd, file.path(out_dir, "pheno_data.csv"))

################################################################################
#   IF samples
################################################################################

if_info = read_tsv(
        if_sr_params_path, col_names = FALSE, show_col_types = FALSE
    ) |>
    mutate(
        src_subject_id = str_extract(X1, '^Br[0-9]{4}'),
        sample_id = sprintf('%s_%s', X2, X3),
        loupe_version = "6.0.0",
        fastq_files = sapply(
            X6,
            function(x) paste(list.files(x, full.names = TRUE), collapse = ';')
        ),
        stain = "IF",
        spaceranger_json = here(
            "processed-data", "01_spaceranger_IF", X1, "outs", "spatial",
            "scalefactors_json.json"
        )
    ) |>
    rename(image_file = X4, alignment_file = X5) |>
    select(
        src_subject_id, sample_id, loupe_version, fastq_files, image_file,
        alignment_file, stain, spaceranger_json
    ) |>
    left_join(pd, by = 'src_subject_id')

################################################################################
#   Combine and export
################################################################################

#   At this point the same info should be collected for H&E and IF, and all
#   values should be non-NULL
stopifnot(!any(is.na(he_info)))
stopifnot(!any(is.na(if_info)))
stopifnot(identical(sort(colnames(if_info)), sort(colnames(he_info))))

if_info = if_info[,colnames(he_info)]
sample_info = rbind(he_info, if_info)

#   All file paths should exist
stopifnot(all(file.exists(sample_info$spaceranger_json)))
stopifnot(all(file.exists(sample_info$image_file)))
stopifnot(all(file.exists(sample_info$alignment_file)))

write_csv(sample_info, file.path(out_dir, "imaging_sample_info.csv"))

session_info()
