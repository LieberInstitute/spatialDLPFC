library(here)
library(readxl)
library(tidyverse)

visium_info_path = here('raw-data', 'sample_info','Visium_dlpfc_mastersheet.xlsx')
visium_dir_new = here('raw-data', 'FASTQ_Globus', 'Visium')

dir.create(visium_dir_new, recursive = TRUE, showWarnings = FALSE)

message('Determining paths of old vs. new FASTQs for Visium...')

#   Construct a tibble with one row per old Visium (non-IF) FASTQ file
fastq_df = tibble(
    fastq_old = list.files(
        here('raw-data', 'FASTQ_renamed'), pattern = '\\.fastq\\.gz$', 
        full.names = TRUE
    )
)

#   Extract just the brain number and position
fastq_df$sample_id_old = fastq_df$fastq_old |>
    basename() |>
    str_extract('Br[0-9]{4}_(ant|mid|post)')

#   Read in sample info for Visium (non-IF); the old sample ID is
#   [brain]_[position], while the new is [slide]_[array]
visium_info = read_excel(visium_info_path) |>
    filter(sequenced == 'yes') |>
    rename(sample_id_old = `sample name`) |>
    mutate(sample_id_new = paste(`slide#`, `array number`, sep = '_'))

fastq_df = fastq_df |>
    #   Add sample info columns (notably 'sample_id_new') to fastq_df
    left_join(visium_info, by = 'sample_id_old') |>
    #   Determine a filepath for the new FASTQ files using sample_id_new
    mutate(
        fastq_new = file.path(
            visium_dir_new,
            paste0(
                sample_id_new,
                basename(fastq_old) |> str_extract('_R[12]_[0-9]*'),
                '.fastq.gz'
            )
        )
    ) |>
    #   Only interested in FASTQ file paths: old and new
    select(c(fastq_old, fastq_new))

message('Copying over Visium FASTQs to new location...')
file.copy(fastq_df$fastq_old, fastq_df$fastq_new)
