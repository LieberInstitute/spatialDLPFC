library(tidyverse)
library(here)
library(readxl)

he_sample_info_path = here(
    'processed-data', 'rdata', 'spe', '01_build_spe', 'spe_sample_info.csv'
)
he_sample_info_2_path = here(
    'raw-data', 'sample_info', 'Visium_dlpfc_mastersheet.xlsx'
)
if_sample_info_path = here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)

if_ids = sprintf('V10B01-087_%s1', c('A', 'B', 'C', 'D'))
if_image_paths = here('raw-data', 'Images', 'VisiumIF', 'VistoSeg', '%s.tif')

col_names = c(
    "subjectkey", "src_subject_id", "interview_date", "interview_age", "sex",
    "comments_misc", "image_file", "image_description", "scan_type",
    "scan_object", "image_file_format", "image_modality",
    "transformation_performed", "transformation_type", "image_history",
    "image_num_dimensions", "image_extent1", "image_extent2", "image_unit1",
    "image_unit2", "image_resolution1", "image_resolution2",
    "emission_wavelength", "objective_magnification", "objective_na",
    "immersion", "exposure_time", "stain", "stain_details",
    "pipeline_stage", "deconvolved", "type_of_microscopy"
)

guids = c(
    "NDAR_INVBM879VT6", "NDAR_INVZN585KGQ", "NDAR_INVEJ479WEA",
    "NDAR_INVCH371AU1", "NDAR_INVGP756DP1", "NDAR_INVCE963WP1",
    "NDAR_INVPY786FH2", "NDAR_INVRD044HHV", "NDAR_INVUM090ZYX",
    "NDAR_INVTG948WW3"
)

################################################################################
#   Read in and preprocess sample info: gather age, sex, ID, and image path
#   for all H&E and IF samples
################################################################################

#-------------------------------------------------------------------------------
#   H&E images
#-------------------------------------------------------------------------------

sample_info = read_csv(he_sample_info_path, show_col_types = FALSE) |>
    rename(interview_age = age) |>
    select(sample_id, sex, interview_age) |>
    mutate(donor = str_extract(sample_id, '^Br[0-9]{4}'))

sample_info_2 = read_excel(he_sample_info_2_path) |>
    head(n = 30) |>
    rename(sample_id = `sample name`) |>
    mutate(
        src_subject_id = sprintf('%s_%s', `slide#`, `array number`),
        image_file = str_replace(
            `image file path`,
            '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC',
            here('raw-data')
        )
    ) |>
    select(sample_id, src_subject_id, image_file)

sample_info = left_join(sample_info, sample_info_2, by = 'sample_id') |>
    select(src_subject_id, donor, sex, interview_age, image_file)

#   Create a version with just one donor per row and immediate phenotype data
#   only
pd = sample_info |>
    group_by(donor) |>
    slice_head(n = 1) |>
    ungroup() |>
    select(donor, sex, interview_age)

#-------------------------------------------------------------------------------
#   IF images
#-------------------------------------------------------------------------------

sample_info_if = read_excel(if_sample_info_path) |>
    head(n = 4) |>
    mutate(
        src_subject_id = sprintf('%s_%s', `Slide SN #`, `Array #`),
        donor = paste0('Br', BrNumbr),
        image_file = sprintf(if_image_paths, src_subject_id)
    ) |>
    left_join(pd, by = 'donor') |>
    select(src_subject_id, donor, sex, interview_age, image_file)

meta_if = tibble(
        donor = if_ids,
        image_file = sprintf(if_image_paths, if_ids),
        image_description = "Leica CS2",
        scan_type = "microscopy",
        scan_object = "Post-mortem",
        image_file_format = "TIFF",
        image_modality = "microscopy",
        transformation_performed = "No",
        transformation_type = NA,
        image_history = NA,
        image_num_dimensions = 2,
        emission_wavelength = "300-700",
        objective_magnification = "40x",
        objective_na = 0.75,
        immersion = 0,
        exposure_time = NA,
        stain = "H&E",
        pipeline_stage = 1,
        deconvolved = 0,
        type_of_microscopy = "BF",
        #   image_extent1, image_extent2, "image_unit1", "image_unit2",
        #   "image_resolution1", "image_resolution2", stain_details
    )

meta = rbind(meta_if, meta_he) |>
    left_join(sample_info, by = 'donor') |>
    select(-donor)
    