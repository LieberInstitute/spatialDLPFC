#   FASTQ files for this dataset have unclear names (such as names that don't
#   include sample ID) and sometimes highly misleading names (such as names
#   including a different sample ID, due to ID-swap issues). This script creates
#   symbolic links where FASTQs have a consistent naming convention, always
#   including the correct sample ID.

library("here")
library("jaffelab")
library("sessioninfo")
library("tidyverse")

sr_table_path <- here("code", "spaceranger", "spaceranger_parameters.txt")
dest_dir <- here("raw-data", "FASTQ_renamed")
df_out_path <- here(
    "processed-data", "synapse_upload", "01-prepare_fastq",
    "fastq_renaming_scheme.csv"
)

dir.create(dest_dir, showWarnings = FALSE)

sr_table <- read.table(sr_table_path)

ids <- sr_table[, 1]
fastq_old <- lapply(strsplit(sr_table[, 6], ","), list.files, full.names = TRUE)

#   The new naming convention for the FASTQ files to upload to synapse will be:
#
#       [sample ID]_[MATE]_[INDEX].fastq.gz
#
#   where [INDEX] is assigned randomly. For example, if there are 3 pairs of
#   FASTQ files for a given sample, one of the pairs will be randomly given the
#   index 1, and so on for index 2 and 3.

fastq_new <- lapply(
    1:nrow(sr_table),
    function(i) {
        num_files <- length(fastq_old[[i]])
        stopifnot(num_files %% 2 == 0) # FASTQ files should come in pairs
        num_pairs <- num_files / 2

        #   Here we take advantage of how 'list.files' orders files, and the
        #   existing naming convention for fastq_old.
        #   Specifically, mates 1 and 2 for a given pair will be adjacent, and
        #   pairs for lanes 1 and 2 will be adjacent.
        new_names <- paste0(
            ids[i],
            rep(c("_R1_", "_R2_"), times = num_pairs),
            rep(1:num_pairs, each = 2),
            ".fastq.gz"
        )

        return(file.path(dest_dir, new_names))
    }
)

#   Create a table with columns ID, old filename, new filename. It will be in
#   "long format" in the sense that an ID will be repeated some multiple of 2
#   times. Save this table so we have a reference for how the old filenames
#   relate to the new ones

ids_long <- ss(basename(unlist(fastq_new)), "_R[12]_")

file_df <- data.frame(
    "sample_id" = ids_long,
    "old_path" = unlist(fastq_old),
    "new_path" = unlist(fastq_new)
)

# write.csv(file_df, file = df_out_path, quote = FALSE, row.names = FALSE)

#   Now symbolically link the old FASTQs to the new destination with the new
#   naming scheme
# all(file.symlink(unlist(fastq_old), unlist(fastq_new)))

#   The column of directories containing FASTQ files also accepts a
#   comma-separated list of FASTQ files. We'll overwrite the column with the
#   new FASTQs in the latter format
sr_table[,6] = fastq_new |>
    sapply(function(x) do.call(paste, as.list(x))) |>
    str_replace_all(' ', ',')

write.table(
    sr_table, file = sr_table_path, quote = FALSE, row.names = FALSE,
    col.names = FALSE
)

session_info()
