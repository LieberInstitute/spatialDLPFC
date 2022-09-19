## for all the results where I subsetted the data for the two clusters of interest

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})

# load spe object
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")

# load results from previous script
res_list <- readRDS(here::here("processed-data", "rdata", "spe", "13_nnSVG", "pairwise", "res_list_subset16_4.rds"))

samples <- names(res_list)

# remove samples that had too few spots to run nnSVG
# remove <- c()
# for(i in seq_along(samples)){
#   if(is.null(dim(res_list[[samples[i]]]))){
#     remove <- append(remove,samples[i])
#   }
# }
# samples <- samples[! samples %in% remove]
# res_list <- res_list[samples]

num_sig_svgs <- data.frame(matrix(ncol = 3, nrow = 0))
# Number of significant SVGs per capture area

for (i in seq_along(samples)) {
    # if(is.null(dim(res_list[[samples[i]]]))){
    #   next
    # }
    temp <- as.data.frame(table(res_list[[samples[i]]]$padj <= 0.05))
    temp$sample_id <- samples[i]
    num_sig_svgs <- rbind(num_sig_svgs, temp)
}

saveRDS(num_sig_svgs, file = here::here("processed-data", "rdata", "spe", "13_nnSVG", "pairwise", "num_sig_svgs_16_4.rds"))


# Create vector of samples for nnSVG on whole tissue. Just rename previous vector.
sample_ids <- samples

# ---------------
# combine results
# ---------------

# sum gene ranks across sample-parts to generate overall ranking

# number of genes that passed filtering for each sample-part
sapply(res_list, nrow)

# match results from each sample-part and store in correct rows
res_ranks <- matrix(NA, nrow = nrow(spe), ncol = length(sample_ids))
rownames(res_ranks) <- rownames(spe)
colnames(res_ranks) <- sample_ids

for (s in seq_along(sample_ids)) {
    stopifnot(colnames(res_ranks)[s] == sample_ids[s])
    stopifnot(colnames(res_ranks)[s] == names(res_list)[s])

    rownames_s <- rownames(res_list[[s]])
    res_ranks[rownames_s, s] <- res_list[[s]][, "rank"]
}

# keep only genes that were not filtered out in all samples
res_ranks <- na.omit(res_ranks)

# calculate average ranks
avg_ranks <- sort(rowMeans(res_ranks))

# summary table
df_summary <- data.frame(
    gene_id = names(avg_ranks),
    gene_name = rowData(spe)[names(avg_ranks), "gene_name"],
    gene_type = rowData(spe)[names(avg_ranks), "gene_type"],
    avg_rank = unname(avg_ranks),
    row.names = names(avg_ranks)
)

head(df_summary, 20)

# Export summary as .csv file
write.csv(df_summary, file = here::here("processed-data", "rdata", "spe", "13_nnSVG", "pairwise", "avg_rank_16_4.csv"), row.names = FALSE)

# Plot top 20 SVGs
SVGs <- df_summary$gene_name[1:20]

# Locate the marker genes
SVG_search <- rowData(spe)$gene_search[match(SVGs, rowData(spe)$gene_name)]

for (i in SVG_search) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "13_nnSVG", "pairwise", "subset_16_4", paste0(gsub("; ", "_", i), ".pdf")),
        assayname = "logcounts",
        minCount = 0,
        cont_colors = c(
            "aquamarine4",
            "springgreen", "goldenrod", "red"
        ),
        alpha = 0.5,
        sample_order = unique(spe$sample_id),
        point_size = 2
    )
}

#################################
# load results from previous script
#################################
res_list <- readRDS(here::here("processed-data", "rdata", "spe", "13_nnSVG", "pairwise", "res_list_subset13_7.rds"))

samples <- names(res_list)

# remove samples that had too few spots to run nnSVG
# remove <- c()
# for(i in seq_along(samples)){
#   if(is.null(dim(res_list[[samples[i]]]))){
#     remove <- append(remove,samples[i])
#   }
# }
# samples <- samples[! samples %in% remove]
# res_list <- res_list[samples]

num_sig_svgs <- data.frame(matrix(ncol = 3, nrow = 0))
# Number of significant SVGs per capture area

for (i in seq_along(samples)) {
    # if(is.null(dim(res_list[[samples[i]]]))){
    #   next
    # }
    temp <- as.data.frame(table(res_list[[samples[i]]]$padj <= 0.05))
    temp$sample_id <- samples[i]
    num_sig_svgs <- rbind(num_sig_svgs, temp)
}

saveRDS(num_sig_svgs, file = here::here("processed-data", "rdata", "spe", "13_nnSVG", "pairwise", "num_sig_svgs_13_7.rds"))


# Create vector of samples for nnSVG on whole tissue. Just rename previous vector.
sample_ids <- samples

# ---------------
# combine results
# ---------------

# sum gene ranks across sample-parts to generate overall ranking

# number of genes that passed filtering for each sample-part
sapply(res_list, nrow)

# match results from each sample-part and store in correct rows
res_ranks <- matrix(NA, nrow = nrow(spe), ncol = length(sample_ids))
rownames(res_ranks) <- rownames(spe)
colnames(res_ranks) <- sample_ids

for (s in seq_along(sample_ids)) {
    stopifnot(colnames(res_ranks)[s] == sample_ids[s])
    stopifnot(colnames(res_ranks)[s] == names(res_list)[s])

    rownames_s <- rownames(res_list[[s]])
    res_ranks[rownames_s, s] <- res_list[[s]][, "rank"]
}

# keep only genes that were not filtered out in all samples
res_ranks <- na.omit(res_ranks)

# calculate average ranks
avg_ranks <- sort(rowMeans(res_ranks))

# summary table
df_summary <- data.frame(
    gene_id = names(avg_ranks),
    gene_name = rowData(spe)[names(avg_ranks), "gene_name"],
    gene_type = rowData(spe)[names(avg_ranks), "gene_type"],
    avg_rank = unname(avg_ranks),
    row.names = names(avg_ranks)
)

head(df_summary, 20)

# Export summary as .csv file
write.csv(df_summary, file = here::here("processed-data", "rdata", "spe", "13_nnSVG", "pairwise", "avg_rank_13_7.csv"), row.names = FALSE)

# Plot top 20 SVGs
SVGs <- df_summary$gene_name[1:20]

# Locate the marker genes
SVG_search <- rowData(spe)$gene_search[match(SVGs, rowData(spe)$gene_name)]

for (i in SVG_search) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "13_nnSVG", "pairwise", "subset_13_7", paste0(gsub("; ", "_", i), ".pdf")),
        assayname = "logcounts",
        minCount = 0,
        cont_colors = c(
            "aquamarine4",
            "springgreen", "goldenrod", "red"
        ),
        alpha = 0.5,
        sample_order = unique(spe$sample_id),
        point_size = 2
    )
}


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
