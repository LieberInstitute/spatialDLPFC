#   While NDA looks in the metadata CSV files to determine files to upload, the
#   CSV files only contain the file basenames, and files are located in many
#   different directories. This script creates symbolic links, located in a
#   single directory, pointing to the files to upload.

library(here)
library(tidyverse)
library(sessioninfo)

dest_dir = here('processed-data', 'synapse_upload', '04-nda', 'to_upload')

upload_list_paths = here(
    'processed-data', 'synapse_upload', '04-nda',
    c(
        'rna_seq_upload_list.txt', 'snp_array_upload_list.txt',
        'visium_image_upload_list.txt'
    )
)

dir.create(dest_dir, showWarnings = FALSE)

src_paths = do.call(c, lapply(upload_list_paths, readLines))
dest_paths = file.path(dest_dir, basename(src_paths))

#   Files are all in the same directory, and therefore must have unique
#   basenames (also source paths should exist)
stopifnot(!any(duplicated(dest_paths)))
stopifnot(all(file.exists(src_paths)))

#   Create symlinks
all(file.symlink(src_paths, dest_paths))

session_info()
