#   Enumerate and quantify the size of:
#       - files currently uploaded to synapse
#       - spatial FASTQs not yet uploaded
#       - spatial images not yet uploaded
#   as of 2024/04/30

library(here)
library(readxl)
library(tidyverse)

sn_man_path = '/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/04_synapse_upload/manifest.tsv'
sample_info_path = here(
    'raw-data', 'sample_info', 'Visium_dlpfc_mastersheet.xlsx'
)
old_repo_path = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC'
n30_fastq_path = here('raw-data', 'FASTQ_Globus', 'Visium')
n19_fastq_path = here('raw-data', 'FASTQ_Globus', 'snRNA-seq')
n4_fastq_path = here('raw-data', 'FASTQ_Globus', 'Visium_SPG')
n4_image_path = here('raw-data', 'Images', 'VisiumIF', 'VistoSeg')
out_dir = here('processed-data', 'synapse_upload', '01-prepare_fastq')

#   Given a vector of file paths, return a numeric vector of file sizes in GB
get_file_size = function(fns) {
    sapply(
        fns,
        function(x) {
            as.numeric(
                system(paste('du -k', x, '| cut -f 1'), intern = TRUE)
            ) / 1e6
        }
    )
}

#   30 H&E images
sample_info = read_excel(sample_info_path) |>
    head(n = 30) |> 
    mutate(
        path = str_replace(`image file path`, old_repo_path, here('raw-data')),
        source = "n30_HE_images",
        size_gb = get_file_size(path)
    ) |>
    select(path, source, size_gb)

#   30 samples of spatial FASTQs
sample_info_temp = tibble(path = list.files(n30_fastq_path, full.names = TRUE)) |>
    mutate(source = "n30_fastqs", size_gb = get_file_size(path))
sample_info = rbind(sample_info, sample_info_temp)

#   snRNA-seq FASTQs
sample_info_temp = tibble(path = list.files(n19_fastq_path, full.names = TRUE)) |>
    mutate(source = "snRNAseq_fastqs", size_gb = get_file_size(path))
sample_info = rbind(sample_info, sample_info_temp)

#   4 samples of spatial FASTQs
sample_info_temp = tibble(path = list.files(n4_fastq_path, full.names = TRUE)) |>
    mutate(source = "n4_fastqs", size_gb = get_file_size(path))
sample_info = rbind(sample_info, sample_info_temp)

#   4 IF images
sample_info_temp = tibble(
        path = list.files(n4_image_path, pattern = '\\.tif$', full.names = TRUE)
    ) |>
    mutate(source = "n4_IF_images", size_gb = get_file_size(path))
sample_info = rbind(sample_info, sample_info_temp)

#   Stats by source
sample_info |>
    group_by(source) |>
    summarize(
        max_size_gb = max(size_gb),
        med_size_gb = median(size_gb),
        total_size_gb = sum(size_gb)
    ) |>
    write_csv(file.path(out_dir, 'to_upload_stats_2024_04_30.csv'))

write_csv(sample_info, file.path(out_dir, 'to_upload_full_2024_04_30.csv'))
