#   In the final workflow, this script is actually run before 01-visium_image.R

library(tidyverse)
library(here)
library(sessioninfo)

image_meta_path = here(
   'processed-data', 'synapse_upload', '04-nda', "imaging_sample_info.csv"
)
image_out_dir = here(
    "processed-data", 'synapse_upload', '04-nda', 'compressed_images'
)
map_path_out = here(
    "processed-data", 'synapse_upload', '04-nda', 'image_mapping.csv'
)

dir.create(image_out_dir, showWarnings = FALSE)

compress_images = function(in_paths, out_paths) {
    stopifnot(all(file.exists(in_paths)))
    stopifnot(all(file.exists(dirname(out_paths))))
    
    for (i in 1:length(in_paths)) {
        command = sprintf('gzip -c %s > %s', in_paths[i], out_paths[i])
        message(sprintf("Running %s...", command))
        system(command)
    }
}

#   Tibble with sample ID, and path to uncompressed and compressed images
image_map = read_csv(image_meta_path, show_col_types = FALSE) |>
    select(sample_id, image_file) |>
    rename(original_path = image_file) |>
    mutate(
        compressed_path = file.path(
            image_out_dir, paste0(sample_id, '.tif.gz')
        )
    )

write_csv(image_map, map_path_out)

# compress_images(image_map$original_path, image_map$compressed_path)

session_info()
