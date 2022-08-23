library(spatialLIBD)
modeling_results <- fetch_data(type = "modeling_results")

sig_genes <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = names(modeling_results)[2],
    reverse = FALSE,
    sce_layer = fetch_data(type = "sce_layer")
)

save(sig_genes, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/sig_genes_manual_annotations.rda")
