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
res_list <- readRDS(here::here("processed-data", "rdata","spe","13_nnSVG", "pairwise","trash","4_1.rds"))

#remove samples that had too few spots to run nnSVG
samples <- names(res_list)
remove <- c()
for(i in seq_along(samples)){
  if(is.null(dim(res_list[[samples[i]]]))){
    remove <- append(remove,samples[i])
  }
}

samples <- samples[! samples %in% remove]
res_list <- res_list[samples]

num_sig_svgs <- data.frame(matrix(ncol = 3, nrow = 0))
# Number of significant SVGs per capture area

for(i in seq_along(samples)){
  # if(is.null(dim(res_list[[samples[i]]]))){
  #   next
  # }
  temp <- as.data.frame(table(res_list[[samples[i]]]$padj<=0.05))
  temp$sample_id <- samples[i]
  num_sig_svgs <- rbind(num_sig_svgs,temp)
}

save(num_sig_svgs, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/13_nnSVG/pairwise/trash/num_sig_svgs_4_1.Rdata")

# [1] "Br2743_ant"
# 
# FALSE  TRUE 
# 2527    54 
# [1] "Br2743_mid"
# 
# FALSE  TRUE 
# 3048    52 
# [1] "Br2743_post"
# 
# FALSE  TRUE 
# 2490     1 
# [1] "Br3942_ant"
# 
# FALSE  TRUE 
# 3118   192 
# [1] "Br3942_mid"
# 
# FALSE  TRUE 
# 2578    80 
# [1] "Br3942_post"
# 
# FALSE  TRUE 
# 4559    12 
# [1] "Br6423_ant"
# 
# FALSE  TRUE 
# 3191    27 
# [1] "Br6423_mid"
# 
# FALSE  TRUE 
# 1842    27 
# [1] "Br6423_post"
# 
# FALSE 
# 3870 
# [1] "Br8492_ant"
# 
# FALSE  TRUE 
# 1752     3 
# [1] "Br8492_mid"
# 
# FALSE  TRUE 
# 1989    22 
# [1] "Br8492_post"
# 
# FALSE  TRUE 
# 1545    20 
# [1] "Br2720_ant"
# 
# FALSE  TRUE 
# 441     1 
# [1] "Br2720_post"
# 
# FALSE  TRUE 
# 1582    11 
# [1] "Br6432_ant"
# 
# FALSE  TRUE 
# 1822   139 
# [1] "Br6432_mid"
# 
# FALSE  TRUE 
# 4119    12 
# [1] "Br6432_post"
# 
# FALSE  TRUE 
# 4476    90 
# [1] "Br6471_ant"
# 
# FALSE  TRUE 
# 2363   422 
# [1] "Br6471_mid"
# 
# FALSE  TRUE 
# 2748   476 
# [1] "Br6471_post"
# 
# FALSE  TRUE 
# 1888   353 
# [1] "Br6522_ant"
# 
# FALSE  TRUE 
# 3232   311 
# [1] "Br6522_mid"
# 
# FALSE  TRUE 
# 4481   597 
# [1] "Br6522_post"
# 
# FALSE  TRUE 
# 1700   494 
# [1] "Br8325_ant"
# 
# FALSE  TRUE 
# 1836    90 
# [1] "Br8325_mid"
# 
# FALSE  TRUE 
# 2450     4 
# [1] "Br8325_post"
# 
# FALSE  TRUE 
# 1255   134 
# [1] "Br8667_ant"
# 
# FALSE  TRUE 
# 1858   619 
# [1] "Br8667_mid"
# 
# FALSE  TRUE 
# 3405    74 
# [1] "Br8667_post"
# 
# FALSE  TRUE 
# 1847   194 


# Create vector of samples for nnSVG on whole tissue
#sample_ids <- num_sig_svgs$sample_id[num_sig_svgs$Var1==TRUE]
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

# directory to save whole tissue results
dir_outputs <- here("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/13_nnSVG/pairwise/trash")
fn_out <- file.path(dir_outputs, "DG_nnSVG_avgrank_4_1")

# Export summary as .csv file
write.csv(df_summary, fn_out, row.names = FALSE)

# Plot top 20 SVGs
SVGs <- df_summary$gene_name[1:20]

# Locate the marker genes
SVG_search <- rowData(spe)$gene_search[match(SVGs, rowData(spe)$gene_name)]

for (i in SVG_search) {
  vis_grid_gene(
    spe = spe,
    geneid = i,
    pdf = here::here("plots", "13_nnSVG", "pairwise", paste0(gsub("; ", "_", i), ".pdf")),
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