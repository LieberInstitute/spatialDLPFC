library(spatialLIBD)

#load parsed modeling results
k = 9
load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("parsed_modeling_results_k",k,".Rdata")))
stats <-modeling_results$enrichment
rownames(stats) = paste0(stats$gene,"_",stats$ensembl)
dim(stats)
# [1] 12225    29 

stats <-
  stats[, grep("[f|t]_stat_", colnames(stats))]


df_enriched <- data.frame(matrix(nrow = 500,ncol = ncol(stats)))
 for(i in 1:ncol(stats)) {
   idx <- order(stats[,i],decreasing = TRUE)[seq_len(500)]
   df_enriched[,i] <- rownames(stats[idx,])
   colnames(df_enriched)[i] <- i
 }

write.csv(df_enriched,
          file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/08_layer_differential_expression/top_500_enriched_genes_k9.csv",
          row.names = TRUE,col.names = TRUE)

#load parsed modeling results
k = 16
load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("parsed_modeling_results_k",k,".Rdata")))
stats <-modeling_results$enrichment
rownames(stats) = paste0(stats$gene,"_",stats$ensembl)
dim(stats)
# [1] 12225    29 

stats <-
  stats[, grep("[f|t]_stat_", colnames(stats))]


df_enriched <- data.frame(matrix(nrow = 500,ncol = ncol(stats)))
for(i in 1:ncol(stats)) {
  idx <- order(stats[,i],decreasing = TRUE)[seq_len(500)]
  df_enriched[,i] <- rownames(stats[idx,])
  colnames(df_enriched)[i] <- i
}

write.csv(df_enriched,
          file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/08_layer_differential_expression/top_500_enriched_genes_k16.csv",
          row.names = TRUE,col.names = TRUE)