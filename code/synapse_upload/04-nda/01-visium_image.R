library(tidyverse)
library(here)
library(readxl)
library(rjson)
library(sessioninfo)
library(tiff)

sample_info_path = here(
   'processed-data', 'synapse_upload', '04-nda', "imaging_sample_info.csv"
)
image_map_path = here(
    "processed-data", 'synapse_upload', '04-nda', 'image_mapping.csv'
)
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
        um_per_px = 65 /
            fromJSON(file = sample_info$spaceranger_json[i])[['spot_diameter_fullres']]

        #   Grab image dimensions from the TIFF's metadata and compute um extent
        meta_list[[i]] = readTIFF(sample_info$image_file[i], payload = FALSE) |>
            as_tibble() |>
            #   Dimension size in pixels
            rename(image_resolution1 = width, image_resolution2 = length) |>
            #   Dimension size in um
            mutate(
                image_extent1 = as.integer(round(image_resolution1 * um_per_px, 0)),
                image_extent2 = as.integer(round(image_resolution2 * um_per_px))
            ) |>
            select(matches('^image_'))
    }
    return(as_tibble(cbind(sample_info, do.call(rbind, meta_list))))
}

################################################################################
#   Read in sample info and add image dimension info
################################################################################

sample_info = read_csv(sample_info_path, show_col_types = FALSE) |>
    add_image_metadata() |>
    #   Replace image path with the compressed version
    select(-image_file) |>
    left_join(
        read_csv(image_map_path, show_col_types = FALSE) |>
            select(-original_path),
        by = 'src_subject_id'
    ) |>
    rename(image_file = compressed_path)

writeLines(
    sample_info$image_file, file.path(out_dir, "visium_image_upload_list.txt")
)

#   NDA validator expects short filename
sample_info = sample_info |> mutate(image_file = basename(image_file))

################################################################################
#   Fill in main metadata expected for the data structure
################################################################################

sample_info = sample_info |>
    mutate(
        interview_date = ifelse(stain == "H&E", '06/25/2020', '06/25/2022'),
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
            stain == "H&E", "40", "20"
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

#   Mimic the submission template from NDA, so this "CSV" can be directly
#   validated with the validator without any reformatting
out_path = file.path(out_dir, 'visium_image.csv')
write_csv(sample_info, out_path)
formatted_info = c(
    paste0('visiumimage,01', paste(rep(',', ncol(sample_info) - 2), collapse = "")),
    readLines(out_path)
)
writeLines(formatted_info, out_path)

session_info()
