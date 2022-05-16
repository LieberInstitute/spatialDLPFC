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

session_info()
