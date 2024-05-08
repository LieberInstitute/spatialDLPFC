library(tidyverse)
library(here)
library(readxl)
library(rjson)
library(sessioninfo)
library(tiff)

he_sample_info_path = here(
    'processed-data', 'rdata', 'spe', '01_build_spe', 'spe_sample_info.csv'
)
he_sample_info_2_path = here(
    'raw-data', 'sample_info', 'Visium_dlpfc_mastersheet.xlsx'
)
if_sample_info_path = here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)
image_map_path = here(
    "processed-data", 'synapse_upload', '04-nda', 'image_mapping.csv'
)
if_image_paths = here('raw-data', 'Images', 'VisiumIF', 'VistoSeg', '%s.tif')
out_dir = here('processed-data', 'synapse_upload', '04-nda')

col_names = c(
    "subjectkey", "src_subject_id", "interview_date", "interview_age", "sex",
    "image_file", "image_description", "scan_type",
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

#   Given a tibble 'sample_info' with columns 'spaceranger_json' and
#   'image_file', return a copy of 'sample_info' with additional columns
#   'image_extent1', 'image_extent2' 'image_resolution1', and
#   'image_resolution1' as required in the VisiumImage NDA data structure
add_image_metadata = function(sample_info) {
    meta_list = list()
    for (i in 1:nrow(sample_info)) {
        #   From: https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs
        #   "The number of pixels that span the diameter of a theoretical 65µm
        #   spot in the original, full-resolution image."
        #   Here, we're computing the image density in micrometers per pixel
        um_per_px = 65 / fromJSON(file = sample_info$spaceranger_json[i])[['spot_diameter_fullres']]

        #   Grab image dimensions from the TIFF's metadata and compute um extent
        meta_list[[i]] = readTIFF(sample_info$image_file[i], payload = FALSE) |>
            as_tibble() |>
            #   Dimension size in pixels
            rename(image_resolution1 = width, image_resolution2 = length) |>
            #   Dimension size in um
            mutate(
                image_extent1 = image_resolution1 * um_per_px,
                image_extent2 = image_resolution2 * um_per_px
            ) |>
            select(matches('^image_'))
    }
    return(cbind(sample_info, do.call(rbind, meta_list)))
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
    mutate(
        sample_id = str_replace(sample_id, '_2$', ''),
        donor = str_extract(sample_id, '^Br[0-9]{4}')
    )

sample_info_2 = read_excel(he_sample_info_2_path) |>
    head(n = 30) |>
    rename(sample_id = `sample name`) |>
    mutate(
        src_subject_id = sprintf('%s_%s', `slide#`, `array number`),
        #   Before compressing images in 05-compress_images.R, the uncompressed
        #   paths were used
        # image_file = str_replace(
        #     `image file path`,
        #     '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC',
        #     here('raw-data')
        # ),
        spaceranger_json = str_replace(
            file.path(
                `spaceranger file path`, 'outs', 'spatial',
                'scalefactors_json.json'
            ),
            '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC',
            here('processed-data')
        )
    ) |>
    select(sample_id, src_subject_id, spaceranger_json)

sample_info = left_join(sample_info, sample_info_2, by = 'sample_id') |>
    select(
        src_subject_id, donor, sex, interview_age, spaceranger_json
    ) |>
    #   Add path to compressed images
    left_join(
        read_csv(image_map_path, show_col_types = FALSE) |>
            select(-original_path),
        by = 'src_subject_id'
    ) |>
    rename(image_file = compressed_path) |>
    #   For internally distinguishing between IF and H&E samples
    mutate(interview_date = '06/25/2020', stain = "H&E")

#-------------------------------------------------------------------------------
#   Phenotype data that applies to many NDA data structures
#-------------------------------------------------------------------------------

#   Just one row per donor
pd = tibble(
        donor = unique(sample_info$donor),
        subjectkey = guids
    ) |>
    left_join(sample_info, by = 'donor', multiple = "any") |>
    select(donor, subjectkey, interview_date, interview_age, sex)

write_csv(pd, file.path(out_dir, "pheno_data.csv"))

#-------------------------------------------------------------------------------
#   H&E: fix bad spaceranger JSON paths
#-------------------------------------------------------------------------------

#   Fix misnamed sample IDs
temp = str_replace(sample_info$spaceranger_json, 'manual_alignment', 'extra_reads')
sample_info$spaceranger_json = ifelse(
    !file.exists(sample_info$spaceranger_json) & file.exists(temp),
    temp, sample_info$spaceranger_json
)

#   Fix misplaced spaceranger outputs
temp = str_replace(
    sample_info$spaceranger_json,
    'NextSeq(/Round[234])?',
    'NextSeq_not_used'
)
sample_info$spaceranger_json = ifelse(
    !file.exists(sample_info$spaceranger_json) & file.exists(temp),
    temp, sample_info$spaceranger_json
)

stopifnot(all(file.exists(sample_info$spaceranger_json)))

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
    mutate(interview_date = '06/25/2022', stain = "IF")

#-------------------------------------------------------------------------------
#   Combine, add guids, and add image dimension info
#-------------------------------------------------------------------------------

sample_info = rbind(sample_info, sample_info_if) |>
    #   Add in guid associated with each donor
    left_join(pd |> select(donor, subjectkey), by = 'donor')

sample_info = add_image_metadata(sample_info)

################################################################################
#   Fill in main metadata expected for the data structure
################################################################################

sample_info = sample_info |>
    mutate(
        image_description = ifelse(
            stain == "H&E", "Leica CS2", "Vectra Polaris + inForm unmixing"
        ),
        scan_type = "microscopy",
        scan_object = "Post-mortem",
        image_file_format = "TIFF",
        image_modality = "microscopy",
        transformation_performed = ifelse(
            stain == "H&E", "No", "Yes"
        ),
        transformation_type = ifelse(
            stain == "H&E", NA, "spectral unmixing"
        ),
        image_history = ifelse(
            stain == "H&E",
            NA, "single channel selected after spectral unmixing"
        ),
        image_num_dimensions = 2,
        image_unit1 = "Micrometers",
        image_unit2 = "Micrometers",
        emission_wavelength = ifelse(
            stain == "H&E", "300-700", "600"
        ),
        objective_magnification = ifelse(
            stain == "H&E", "40x", "20x"
        ),
        objective_na = 0.75,
        immersion = 0,
        exposure_time = ifelse(
            stain == "H&E", NA, 0.0021 # "DAPI: 0.0021; Opal 520: 0.143; Opal 570: 0.330; Opal 620: 0.2; Opal 690: 1.07; Autofluorescence: 0.1"
        ),
        stain_details = ifelse(
            stain == "H&E",
            "H&E staining was conducted according to the manufacturer's instructions (protocol CG000160, Rev B, 10x Genomics)",
            "Immunofluorescence (IF) staining was conducted according to the manufacturer’s instruction (catalog no.CG000312 Rev C, 10x Genomics)"
        ),
        pipeline_stage = ifelse(
            stain == "H&E", 1, 2
        ),
        deconvolved = 0,
        type_of_microscopy = ifelse(
            stain == "H&E", "BF", "WFF"
        )
    ) |>
    select(all_of(col_names))

write_csv(sample_info, file.path(out_dir, 'VisiumImage.csv'))

session_info()
