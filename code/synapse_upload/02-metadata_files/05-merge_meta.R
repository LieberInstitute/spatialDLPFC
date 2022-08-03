#   This synapse project involves uploading data from two different repos (this
#   one and 'DLPFC_snRNAseq'). It ended up being simpler to just create two sets
#   of metadata files for synapse (one in each repo), but for the actual upload,
#   only a single set (with the exception of the need for 2 'assays.csv' files)
#   is needed. This script merges the two sets into the one needed for upload.

library('here')
library("sessioninfo")

spatial_dir = here("processed-data", "synapse_upload", "02-metadata_files")
sc_dir = '/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/04_synapse_upload'
out_dir = here(
    "processed-data", "synapse_upload", "02-metadata_files", "combined"
)

#   Note: assays files and individual file are already copied in the shell
#   script wrapper

################################################################################
#   Biospecimen
################################################################################

spatial_meta = read.csv(file.path(spatial_dir, 'biospecimen.csv'))
sc_meta = read.csv(file.path(sc_dir, 'biospecimen.csv'))
meta = rbind(spatial_meta, sc_meta)
write.csv(meta, file.path(out_dir, 'biospecimen.csv'), row.names = FALSE)

################################################################################
#   Manifest
################################################################################

spatial_meta = read.table(file.path(spatial_dir, 'manifest.tsv'), header = TRUE)
sc_meta = read.table(file.path(sc_dir, 'manifest.tsv'), header = TRUE)

#   Fix paths the assay metadata files
spatial_meta[basename(spatial_meta$path) == 'assay.csv', 'path'] = file.path(
    out_dir, 'assay_TODO.csv'
)
sc_meta[basename(sc_meta$path) == 'assay.csv', 'path'] = file.path(
    out_dir, 'assay_scrnaSeq.csv'
)

#   Use the path to the merged individual and biospecimen files only
spatial_meta[
    basename(spatial_meta$path) %in% c('individual.csv', 'biospecimen.csv'),
    'path'
] = 
    file.path(
        out_dir, c('individual.csv', 'biospecimen.csv')
    )
sc_meta = sc_meta[
    - which(basename(sc_meta$path) %in% c('individual.csv', 'biospecimen.csv')),
]

meta = rbind(spatial_meta, sc_meta)
write.table(
    meta, file.path(out_dir, 'manifest.tsv'), row.names = FALSE, sep = "\t"
)

session_info()
