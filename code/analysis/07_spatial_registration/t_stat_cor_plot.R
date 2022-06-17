library(spatialLIBD)

k = 7
load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("parsed_modeling_results_k",k,".Rdata")))
my_enrich_tstats <-modeling_results$enrichment
rownames(my_enrich_tstats) = my_enrich_tstats$ensembl

ground_truth <- spatialLIBD::fetch_data("modeling_results")
ground_truth_enrich <-ground_truth$enrichment
rownames(ground_truth_enrich) = ground_truth_enrich$ensembl

tstats <-
  ground_truth_enrich[, grep("[f|t]_stat_", colnames(ground_truth_enrich))]
colnames(tstats) <-
  gsub("[f|t]_stat_", "", colnames(tstats))

top_n_index <- unique(as.vector(apply(tstats, 2, function(t) {
  order(t, decreasing = TRUE)[seq_len(100)]
})))

tstats <- tstats[top_n_index, , drop = FALSE]

mm<-match(rownames(tstats),rownames(my_enrich_tstats))

tstats <- tstats[!is.na(mm), ]
my_enrich_tstats <- my_enrich_tstats[mm[!is.na(mm)], ]

top_genes <- rownames(tstats)

#make single dataframe 

pdf(file = here::here("plots","07_spatial_registration","t_cor_k7_wm.pdf"))
plot(my_enrich_tstats$t_stat_7,tstats$WM)
dev.off()

#try to plot all genes
k = 7
load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("parsed_modeling_results_k",k,".Rdata")))
my_enrich_tstats <-modeling_results$enrichment
rownames(my_enrich_tstats) = my_enrich_tstats$ensembl

ground_truth <- spatialLIBD::fetch_data("modeling_results")
ground_truth_enrich <-ground_truth$enrichment
rownames(ground_truth_enrich) = ground_truth_enrich$ensembl

common_genes <- intersect(rownames(ground_truth_enrich),rownames(my_enrich_tstats))
ground_truth_enrich <-ground_truth_enrich[common_genes,]
my_enrich_tstats <- my_enrich_tstats[common_genes,]

pdf(file = here::here("plots","07_spatial_registration","t_cor_k7_wm_all_genes.pdf"))
plot(my_enrich_tstats$t_stat_7,tstats$WM)
dev.off()

