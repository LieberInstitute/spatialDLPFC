library(spatialLIBD)

load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("parsed_modeling_results_k",k,".Rdata")))
my_enrich_tstats <-modeling_results$enrichment
rownames(my_enrich_tstats) = my_enrich_tstats$ensembl

ground_truth <- spatialLIBD::fetch_data("modeling_results")
ground_truth_enrich <-ground_truth$enrichment
rownames(ground_truth_enrich) = ground_truth_enrich$ensembl

common_genes <- intersect(rownames(my_enrich_tstats),rownames(ground_truth_enrich))

ground_truth_enrich <- ground_truth_enrich[common_genes,]
my_enrich_tstats <- my_enrich_tstats[common_genes,]
#must add rownames and intersect 
pdf(file = "~/tcor.pdf")
plot(my_enrich_tstats$t_stat_1,ground_truth_enrich$t_stat_WM)
dev.off()
