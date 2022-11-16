library(spatialLIBD)

# load parsed modeling results
k <- 9
load(file = here::here("processed-data", "rdata", "spe", "07_layer_differential_expression", paste0("parsed_modeling_results_k", k, ".Rdata")))
stats <- modeling_results$enrichment
rownames(stats) <- paste0(stats$gene, "_", stats$ensembl)
dim(stats)
# [1] 12225    29

stats <-
    stats[, grep("[f|t]_stat_", colnames(stats))]


df_enriched <- data.frame(matrix(nrow = 500, ncol = ncol(stats)))
for (i in 1:ncol(stats)) {
    idx <- order(stats[, i], decreasing = TRUE)[seq_len(500)]
    df_enriched[, i] <- rownames(stats[idx, ])
    colnames(df_enriched)[i] <- i
}

write.csv(df_enriched,
    file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/07_layer_differential_expression/top_500_enriched_genes_k9.csv",
    row.names = TRUE, col.names = TRUE
)

# load parsed modeling results
k <- 16
load(file = here::here("processed-data", "rdata", "spe", "07_layer_differential_expression", paste0("parsed_modeling_results_k", k, ".Rdata")))
stats <- modeling_results$enrichment
rownames(stats) <- paste0(stats$gene, "_", stats$ensembl)
dim(stats)
# [1] 12225    29

stats <-
    stats[, grep("[f|t]_stat_", colnames(stats))]


df_enriched <- data.frame(matrix(nrow = 500, ncol = ncol(stats)))
for (i in 1:ncol(stats)) {
    idx <- order(stats[, i], decreasing = TRUE)[seq_len(500)]
    df_enriched[, i] <- rownames(stats[idx, ])
    colnames(df_enriched)[i] <- i
}

write.csv(df_enriched,
    file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/07_layer_differential_expression/top_500_enriched_genes_k16.csv",
    row.names = TRUE, col.names = TRUE
)

#### top 25 genes from pairwise model k =16

library(spatialLIBD)
load(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/07_layer_differential_expression/parsed_modeling_results_k16.Rdata")

sig_genes <- sig_genes_extract(
    n = 25,
    modeling_results = modeling_results,
    model_type = names(modeling_results)[3],
    reverse = FALSE,
    sce_layer = fetch_data(type = "sce_layer")
)

# Pairwise: 5 vs. 9 (Layer 4?)
# Pairwise: 4 vs. 16 (Layer 5?)
# Pairwise: 7 vs. 13 (Layer "6A" vs. "6B")
# Pairwise: 12 vs. 13 (layer 6 apex vs base)
tests <- c("BayesSpace5-BayesSpace9", "BayesSpace4-BayesSpace16", "BayesSpace7-BayesSpace13", "BayesSpace12-BayesSpace13")

sig_genes <- sig_genes[which(sig_genes$test %in% tests), ]
write.csv(sig_genes,
    file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/07_layer_differential_expression/top_25_pairwise_k16_subset.csv",
    row.names = TRUE, col.names = TRUE
)
