library('here')
library('ggplot2')
library('jaffelab')
library('tidyverse')

cell_group = "broad" # "broad" or "layer"

sample_ids = here(
    "processed-data", "spot_deconvo", "05-shared_utilities", "IF",
    "sample_ids.txt"
)

deconvo_tools = c('01-tangram', '03-cell2location', '04-spotlight')

#   "Ground-truth" cell counts from cellpose + trained classification tree
actual_paths = here(
    "processed-data", "spot_deconvo", "02-cellpose", "{sample_id}", "clusters.csv"
)

#   Cell counts estimated by each deconvolution tool
observed_paths <- here(
    "processed-data", "spot_deconvo", "{deconvo_tool}", "IF", cell_group,
    "{sample_id}", "clusters.csv"
)

cell_types_actual = c('micro', 'neuron', 'oligo', 'other')
if (cell_group == "broad") {
    cell_types = c('Astro', 'Excit', 'Inhib', 'Micro', 'Oligo', 'OPC')
} else {
    cell_types = c() # TODO
}

sample_ids = readLines(sample_ids)

actual_list = list()
observed_list = list()
index = 1

added_colnames = c("barcode", "sample_id", "deconvo_tool", "obs_type")

for (sample_id in sample_ids) {
    #   Read in the "ground-truth" counts for this sample
    actual_path = sub("\\{sample_id\\}", sample_id, actual_paths)
    actual_df_small = read.csv(actual_path)
    
    #   Make sure we have just counts for each cell type, barcode, and
    #   sample ID variables-- nothing else. Order column names
    actual_df_small$barcode = ss(actual_df_small$key, '_', 1)
    actual_df_small$sample_id = sample_id
    actual_df_small$deconvo_tool = "none"
    actual_df_small$obs_type = "actual"
    actual_list[[index]] = actual_df_small[
        , c(added_colnames, cell_types_actual)
    ]
    
    for (deconvo_tool in deconvo_tools) {
        #   Read in estimated cell counts for this deconvo tool and sample
        observed_path = sub("\\{sample_id\\}", sample_id, observed_paths)
        observed_path = sub("\\{deconvo_tool\\}", deconvo_tool, observed_path)
        observed_df_small = read.csv(observed_path)
        
        #   Make sure column names are consistent and include only info about
        #   barcode, sample_id, deconvo tool, and cell-type counts
        observed_df_small$barcode = ss(observed_df_small$key, '_', 1)
        observed_df_small$sample_id = sample_id
        observed_df_small$deconvo_tool = deconvo_tool
        observed_df_small$obs_type = "observed"
        observed_df_small = observed_df_small[
            , c(added_colnames, cell_types)
        ]
        
        observed_list[[index]] = observed_df_small
        index = index + 1
    }
}

#   Form data frames containing cell-type counts for all spots, samples, and
#   deconvolution tools
actual_df = as_tibble(do.call(rbind, actual_list))
observed_df = as_tibble(do.call(rbind, observed_list))

#   Combine cell types as appropriate for comparison against the relatively
#   narrow types in the ground-truth
if (cell_group == "broad") {
    colnames(observed_df) = tolower(colnames(observed_df))
    
    observed_df = observed_df %>%
        mutate("neuron" = excit + inhib) %>%
        mutate("other" = astro + opc) %>%
        select(all_of(c(added_colnames, cell_types_actual)))
    
    stopifnot(all(colnames(observed_df) == colnames(actual_df)))
} else {
    #   TODO
}

#   TODO: (might use reshape2::cast) reshape full_df such that instead of the
#   'obs_type' column and 'count', we have 'count_actual' and 'count_observed'
#   (and thus no need for the 'obs_type' column)
full_df = rbind(observed_df, actual_df) %>%
    melt(
        id.vars = added_colnames, variable.name = "cell_type",
        value.name = "count"
    )
