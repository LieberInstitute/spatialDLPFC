library('here')
library('ggplot2')

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

sample_ids = readLines(sample_ids)

actual_props_list = list()
for (sample_id in sample_ids) {
    #   Read in the "ground-truth" counts for this sample
    actual_path = sub("{sample_id}", sample_id, actual_paths)
    actual_df = read.csv(actual_path)
    
    #   Grab only the columns with counts for each cell type
    actual_df = actual_df[
        , -match(c("key", "sample", "cell_count"), colnames(actual_df))
    ]
    
    #   Make sure colnames are kept so we can verify cell types are in order!
    #   Extract the proportion of each cell type in all spots combined
    actual_counts = colSums(actual_df)
    actual_props_list[[sample_id]] = actual_counts / sum(actual_counts)
    
    # observed_list = list()
    # for (deconvo_tool in deconvo_tools) {
    #     #   Read in estimated cell counts for this deconvo tool and sample
    #     observed_path = sub("{sample_id}", sample_id, observed_paths)
    #     observed_path = sub("{deconvo_tool}", deconvo_tool, observed_path)
    #     observed_df = read.csv(observed_path)
    # }
    # 
    # observed_df = do.call(rbind, observed_list)
}

#   For ggplot, we want a dataframe with columns (for the first scatterplot in #99):
#   prop_cell_type1, prop_cell_type2, ... , sample_id, deconvo_tool,
#   observed count, actual count
#
#   Even better would be a dataframe with the additional columns:
#   cell group, spot barcode
#
#   This would be really large though-- not memory efficient
