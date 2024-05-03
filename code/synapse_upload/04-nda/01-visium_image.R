library(tidyverse)
library(here)

sample_info_path = here(
    'processed-data', 'rdata', 'spe', '01_build_spe', 'spe_sample_info.csv'
)
if_info_path = here(
    'processed-data', 'rdata', 'spe_IF', '01_build_spe_IF',
    'spe_IF_sample_info.csv'
)
if_ids = sprintf('V10B01-087_%s1', c('A', 'B', 'C', 'D'))
if_image_paths = here('raw-data', 'Images', 'VisiumIF', 'VistoSeg', '%s.tif')

col_names = c(
    "subjectkey", "src_subject_id", "interview_date", "interview_age", "sex",
    "comments_misc", "image_file", "image_description", "scan_type",
    "scan_object", "image_file_format", "data_file2", "data_file2_type",
    "image_modality", "transformation_performed", "transformation_type",
    "image_history", "image_num_dimensions", "image_extent1", "image_extent2",
    "image_unit1", "image_unit2", "image_resolution1", "image_resolution2",
    "manifest", "emission_wavelength", "objective_magnification",
    "objective_na", "immersion", "exposure_time", "stain", "stain_details",
    "pipeline_stage", "deconvolved", "decon_software", "decon_method",
    "psf_type", "psf_file", "decon_snr", "decon_iterations",
    "micro_temmplate_name", "in_stack", "decon_template_name", "stack",
    "slices", "slice_number", "slice_thickness", "type_of_microscopy",
    "excitation_wavelength"
)

guids = c(
    "NDAR_INVBM879VT6", "NDAR_INVZN585KGQ", "NDAR_INVEJ479WEA",
    "NDAR_INVCH371AU1", "NDAR_INVGP756DP1", "NDAR_INVCE963WP1",
    "NDAR_INVPY786FH2", "NDAR_INVRD044HHV", "NDAR_INVUM090ZYX",
    "NDAR_INVTG948WW3"
)

sample_info = read_csv(sample_info_path, show_col_types = FALSE) |>
    rename(src_subject_id = sample_id, interview_age = age) |>
    select(src_subject_id, sex, interview_age) |>
    mutate(donor = str_extract(src_subject_id, '^Br[0-9]{4}'))

meta_if = tibble(
        donor = if_ids,
        image_file = sprintf(if_image_paths, if_ids),
        image_description = "Leica CS2",
        scan_type = "microscopy",
        scan_object = "Post-mortem",
        image_file_format = "TIFF",
    ) |>
    left_join(sample_info, by = 'donor')