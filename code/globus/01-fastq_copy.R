library(here)
library(readxl)
library(tidyverse)
library(jaffelab)
library(sessioninfo)

visium_info_path = here('raw-data', 'sample_info','Visium_dlpfc_mastersheet.xlsx')
visium_spg_info_path = here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)

visium_dir_new = here('raw-data', 'FASTQ_Globus', 'Visium')
visium_spg_dir_new = here('raw-data', 'FASTQ_Globus', 'Visium_SPG')
sn_dir_new = here('raw-data', 'FASTQ_Globus', 'snRNA-seq')
sn_dir_old = '/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/raw-data/FASTQ_renamed'

out_table_path = here('raw-data', 'FASTQ_Globus', 'fastq_mapping.csv')

perform_copy = FALSE

dir.create(visium_dir_new, recursive = TRUE, showWarnings = FALSE)
dir.create(visium_spg_dir_new, recursive = FALSE, showWarnings = FALSE)
dir.create(sn_dir_new, recursive = FALSE, showWarnings = FALSE)

################################################################################
#   Visium
################################################################################

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
    mutate(
        sample_id_new = paste(`slide#`, `array number`, sep = '_'),
        donor = ss(sample_id_old, '_')
    )

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
    #   Interested in FASTQ paths, sample ID, and donor
    select(c(fastq_old, fastq_new, sample_id_new, donor))

#   Final sanity checks before copying
stopifnot(all(file.exists(fastq_df$fastq_old)))
stopifnot(length(unique(fastq_df$fastq_new)) == length(fastq_df$fastq_new))

if (perform_copy) {
    message('Copying over Visium FASTQs to new location...')
    file.copy(fastq_df$fastq_old, fastq_df$fastq_new)
}

#   Maintain a tibble keeping track of 3 versions of FASTQ file paths: original
#   from sequencing, symlinks for synapse upload, and FASTQs shared through
#   Globus
fastq_map = fastq_df |>
    rename(
        fastq_globus = fastq_new, fastq_synapse = fastq_old,
        sample_id = sample_id_new
    ) |>
    mutate(fastq_sequencing = normalizePath(fastq_synapse), assay = "Visium")

################################################################################
#   Visium-SPG
################################################################################

message('Determining paths of old vs. new FASTQs for Visium-SPG...')

visium_spg_info = read_excel(visium_spg_info_path) |>
    filter(!is.na(`Sample #`)) |>
    mutate(
        position = ss(Br_Region, '_', 2),
        sample_id_old = paste0('Br', BrNumbr, '_', position, '_IF'),
        sample_id_new = paste(`Slide SN #`, `Array #`, sep = '_'),
        donor = paste0('Br', BrNumbr)
    )

fastq_df = tibble(
    fastq_old = list.files(
        here('raw-data', 'FASTQ', 'Visium_IF'), pattern = '\\.fastq\\.gz$', 
        recursive = TRUE, full.names = TRUE
    )
)
fastq_df$sample_id_old = fastq_df$fastq_old |>
    str_extract('Br[0-9]{4}_(Ant|Post)_IF')

fastq_df = fastq_df |>
    #   Add sample info columns (notably 'sample_id_new') to fastq_df
    left_join(visium_spg_info, by = 'sample_id_old') |>
    #   Determine a filepath for the new FASTQ files using sample_id_new
    mutate(
        fastq_new = file.path(
            visium_spg_dir_new,
            paste0(
                sample_id_new,
                basename(fastq_old) |> str_extract('_R[12]'),
                '.fastq.gz'
            )
        )
    ) |>
    #   Interested in FASTQ paths, sample ID, and donor
    select(c(fastq_old, fastq_new, sample_id_new, donor))

#   Actually, a single sample sometimes has more than one FASTQ per mate.
#   Deduplicate filenames accordingly (just add a suffix '_1', '_2', etc)
for (fastq_path in unique(fastq_df$fastq_new)) {
    num_copies = length(which(fastq_df$fastq_new == fastq_path))
    base_path = ss(
        basename(fastq_df$fastq_new[fastq_df$fastq_new == fastq_path]), '\\.'
    )

    fastq_df$fastq_new[fastq_df$fastq_new == fastq_path] = file.path(
        visium_spg_dir_new,
        paste0(base_path, '_', 1:num_copies, '.fastq.gz')
    )
}

#   Final sanity checks before copying
stopifnot(all(file.exists(fastq_df$fastq_old)))
stopifnot(length(unique(fastq_df$fastq_new)) == length(fastq_df$fastq_new))

if (perform_copy) {
    message('Copying over Visium-SPG FASTQs to new location...')
    file.copy(fastq_df$fastq_old, fastq_df$fastq_new)
}

#   Add Visium-SPG FASTQs to the big table of FASTQs
fastq_map = fastq_df |>
    rename(
        fastq_globus = fastq_new, sample_id = sample_id_new
    ) |>
    mutate(
        fastq_synapse = fastq_globus,
        fastq_sequencing = normalizePath(fastq_old),
        assay = "Visium-SPG"
    ) |>
    select(!fastq_old) |>
    rbind(fastq_map)

################################################################################
#   snRNA-seq
################################################################################

#   The old FASTQ paths already use the right naming convention and just need
#   to be moved to the new parent directory
fastq_old = list.files(
    sn_dir_old, pattern = '\\.fastq\\.gz$', full.names = TRUE
)
fastq_new = file.path(sn_dir_new, basename(fastq_old))

#   Final sanity checks before copying
stopifnot(all(file.exists(fastq_old)))
stopifnot(length(unique(fastq_new)) == length(fastq_new))

if (perform_copy) {
    message('Copying over snRNA-seq FASTQs to new location...')
    file.copy(fastq_old, fastq_new)
}

#   Add snRNA-seq FASTQs to the table of all FASTQs
fastq_map = tibble(
    fastq_synapse = fastq_old,
    fastq_globus = fastq_new
) |>
    mutate(
        donor = ss(basename(fastq_new), '_'),
        sample_id = paste(donor, ss(basename(fastq_new), '_', 2), sep = '_'),
        assay = 'snRNA-seq',
        fastq_sequencing = normalizePath(fastq_synapse)
    ) |>
    rbind(fastq_map) |>
    relocate(sample_id, donor, starts_with('fastq'), assay)

write.csv(fastq_map, out_table_path, quote = FALSE, row.names = FALSE)

session_info()
