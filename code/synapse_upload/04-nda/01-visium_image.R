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
#   Functions
################################################################################

get_image_metadata = function(image_paths) {
    for (image_path in image_paths) {
        readTIFF(image_path, payload = FALSE) |>
            as_tibble() |>
            rename(image_extent1 = width, image_extent2 = length) |>

    }
}

################################################################################
#   Read in and preprocess sample info: gather age, sex, ID, image path, and
#   spaceranger JSON for all H&E and IF samples
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
        ),
        spaceranger_json = file.path(
            `spaceranger file path`, 'outs', 'spatial', 'scalefactors_json.json'
        )
    ) |>
    select(sample_id, src_subject_id, image_file, spaceranger_json)

sample_info = left_join(sample_info, sample_info_2, by = 'sample_id') |>
    select(
        src_subject_id, donor, sex, interview_age, image_file, spaceranger_json
    ) |>
    #   For internally distinguishing between IF and H&E samples
    mutate(interview_date = '06/25/2020', image_type = "H&E")

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
        image_file = sprintf(if_image_paths, src_subject_id),
        sample_id = sprintf(
            '%s_%s_IF', donor, str_extract(Br_Region, '(Ant|Post)')
        ),
        spaceranger_json = here(
            'processed-data', '01_spaceranger_IF', sample_id, 'outs', 'spatial',
            'scalefactors_json.json'
        )
    ) |>
    left_join(pd, by = 'donor') |>
    select(
        src_subject_id, donor, sex, interview_age, image_file, spaceranger_json
    ) |>
    #   For internally distinguishing between IF and H&E samples
    mutate(interview_date = '06/25/2022', image_type = "IF")

#-------------------------------------------------------------------------------
#   Combine and add guids
#-------------------------------------------------------------------------------

sample_info = rbind(sample_info, sample_info_if) |>
    #   Add in guid associated with each donor
    left_join(tibble(donor = pd$donor, subjectkey = guids), by = 'donor') |>
    select(-donor)

################################################################################
#   Fill in main metadata expected for the data structure
################################################################################

meta_if = tibble(
        image_file = sample_info$image_file,
        image_description = ifelse(
            sample_info$image_type == "H&E",
            "Leica CS2", "Vectra Polaris + inForm unmixing"
        ),
        scan_type = "microscopy",
        scan_object = "Post-mortem",
        image_file_format = "TIFF",
        image_modality = "microscopy",
        transformation_performed = ifelse(
            sample_info$image_type == "H&E", "No", "Yes"
        ),
        transformation_type = ifelse(
            sample_info$image_type == "H&E", NA, "spectral unmixing"
        ),
        image_history = ifelse(
            sample_info$image_type == "H&E",
            NA, "single channel selected after spectral unmixing"
        ),
        image_num_dimensions = 2,
        emission_wavelength = ifelse(
            sample_info$image_type == "H&E", "300-700", "600"
        ),
        objective_magnification = ifelse(
            sample_info$image_type == "H&E", "40x", "20x"
        ),
        objective_na = 0.75,
        immersion = 0,
        exposure_time = ifelse(
            sample_info$image_type == "H&E",
            NA, 0.0021 # "DAPI: 0.0021; Opal 520: 0.143; Opal 570: 0.330; Opal 620: 0.2; Opal 690: 1.07; Autofluorescence: 0.1"
        ),
        stain = sample_info$image_type,
        stain_details = ifelse(
            sample_info$image_type == "H&E",
            TODO,
            "Immunofluorescence (IF) staining was conducted according to the manufacturer’s instruction (catalog no.CG000312 Rev C, 10x Genomics). In brief, tissue blocks were cryosectioned at 10μm thickness and tissue sections were collected on a Visium Spatial Gene Expression Slide (catalog no. 2000233, 10x Genomics). Tissue was then fixed in pre-chilled methanol, treated with BSA-containing blocking buffer, and incubated for 30 minutes at room temperature with primary antibodies against the four cell-type marker proteins, NeuN for neurons, TMEM119 for microglia, GFAP for astrocytes, and OLIG2 for oligodendrocytes. For primary antibodies we used mouse anti-NeuN antibody conjugated to Alexa 488 (Sigma Aldrich, Cat# MAB377X, 1:100), rabbit anti-TMEM119 antibody (Sigma Aldrich, Cat# HPA051870, 1:20), rat anti-GFAP antibody (Thermofisher, Cat# 13-0300, 1:100), and goat anti-OLIG2 antibody (R&D systems, Cat# AF2418, 1:20). Following a total of 5 washes, secondary antibodies were applied for 30 minutes at room temperature. Detailed product information of the secondary antibodies is provided as follows: donkey anti-rabbit IgG conjugated to Alexa 555 (Thermofisher, Cat# A-31572, 1:300), donkey anti-rat IgG conjugated to Alexa 594 (Thermofisher, Cat# A-21209, 1:600), and donkey anti-goat IgG conjugated to Alexa 647 (Thermofisher, Cat# A-21447, 1:400). DAPI (Thermofisher, Cat# D1306, 1:3000, Final 1.67 μg/ml) was used for nuclear counterstaining."
        ),
        pipeline_stage = ifelse(
            sample_info$image_type == "H&E", 1, 2
        ),
        deconvolved = 0,
        type_of_microscopy = ifelse(
            sample_info$image_type == "H&E", "BF", "WFF"
        )
        #   image_extent1, image_extent2, "image_unit1", "image_unit2",
        #   "image_resolution1", "image_resolution2"
    )

meta = rbind(meta_if, meta_he) |>
    left_join(sample_info, by = 'donor') |>
    select(-donor)
    