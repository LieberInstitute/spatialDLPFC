#   FASTQ files for this dataset have unclear names (such as names that don't
#   include sample ID) and sometimes highly misleading names (such as names
#   including a different sample ID, due to ID-swap issues). This script creates
#   symbolic links where FASTQs have a consistent naming convention, always
#   including the correct sample ID.

library('here')

sr_table_path = here('code', 'spaceranger', 'spaceranger_parameters.txt')
dest_dir = here('code', 'synapse_upload', 'FASTQ')

sr_table = read.table(sr_table_path)
ids = sr_table[,1]

#   For each row (sample), check if:
#       1. every directory in the 6th column contains the sample ID (in col 1)
#       2. the FASTQs in those directories contain the sample ID
#   Trying to check if information about the sample ID is contained in FASTQ
#   full file paths. I manually found that at least one sample has FASTQ file
#   paths lacking the needed information

for (i in 1:nrow(sr_table)) {
    contains_id = FALSE
    id = sr_table[i, 1]
    print(paste0('Checking sample ID ', id, ':'))
    fastq_dirs = strsplit(sr_table[i, 6], ',')[[1]]
    contains_id = any(grepl(id, basename(fastq_dirs)))
    if (! contains_id) {
        print('Did not find ID in FASTQ dirnames!')
    }
}
