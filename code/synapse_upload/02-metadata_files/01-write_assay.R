#  Create the 'assay.csv' file, 1 of 3 required metadata files for upload to
#  synapse

library('sessioninfo')
library('here')
library('readxl')

template_path = here(
    'code', 'synapse_upload', 'templates', 'template_assay_rnaSeq.xlsx'
)

fastq_naming_path = here(
    'processed-data', 'synapse_upload', '01-prepare_fastq',
    'fastq_renaming_scheme.csv'
)

write_path = here(
    'processed-data', 'synapse_upload', '02-metadata_files', 'assay.csv'
)

#  Ensure we included all the required columns
template_names = colnames(read_excel(template_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.csv(meta_df, write_path, row.names=FALSE)

session_info()
