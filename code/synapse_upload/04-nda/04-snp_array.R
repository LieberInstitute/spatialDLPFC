library(tidyverse)
library(here)
library(readxl)
library(sessioninfo)

geno_dir = '/dcs05/lieber/liebercentral/libdGenotype_LIBD001/BrainGenotyping/subsets/sc_n10_Nick'
pd_path <- here("processed-data", "synapse_upload", "04-nda", "pheno_data.csv")
out_dir = here('processed-data', 'synapse_upload', '04-nda')

#   The set of columns we have information for in the SNP Array data structure
col_names <- c(
    "subjectkey", "src_subject_id", "interview_date", "interview_age", "sex",
    "experiment_id", "referenceset", "data_file1", "data_file1_type",
    "genotyping_chip_type", "study", "psych_enc_exclude",
    "psych_enc_exclude_reason", "platform", "assay"
)

#   Map of donors to sample IDs used for genotyping data
id_map = read_tsv(
        file.path(geno_dir, 'sc_n10.brnum2genoID.tab'), show_col_types = FALSE,
        col_names = FALSE
    ) |>
    rename(src_subject_id = X1, sample_id = X2)

#   General donor-level fields we need
pd = read_csv(pd_path, show_col_types = FALSE)

#   Gather, define, and export the necessary metadata to CSV
meta_df = tibble(
        data_file1 = list.files(
            file.path(geno_dir, 'idats'), full.names = TRUE
        ),
        data_file1_type = "Intensity data file"
    ) |>
    mutate(sample_id = str_extract(data_file1, '[0-9]{12}_R[0-9]{2}C01')) |>
    left_join(id_map, by = 'sample_id') |>
    left_join(pd, by = 'src_subject_id') |>
    mutate(
        experiment_id = 2607,
        referenceset = 3, # GRCh37
        genotyping_chip_type = ifelse(
            src_subject_id == "Br8325",
            "Infinium Omni2.5-8 v1.4",
            "Infinium Omni2.5-8 v1.5"
        ),
        study = "LIBD spatial DLPFC",
        psych_enc_exclude = 0, # not excluded
        psych_enc_exclude_reason = "Not excluded",
        platform = "Illumina genotyping microarray kit",
        assay = 14 # snpArray
    ) |>
    select(all_of(col_names))

writeLines(meta_df$data_file1, file.path(out_dir, 'snp_array_upload_list.txt'))

meta_df = meta_df |>
    mutate(data_file1 = basename(data_file1))

out_path = file.path(out_dir, 'snp_array.csv')

#   Mimic the submission template from NDA, so this "CSV" can be directly
#   validated with the validator without any reformatting
write_csv(meta_df, out_path)
formatted_info = c(
    paste0('snp_array,01', paste(rep(',', ncol(meta_df) - 2), collapse = "")),
    readLines(out_path)
)
writeLines(formatted_info, out_path)

session_info()
