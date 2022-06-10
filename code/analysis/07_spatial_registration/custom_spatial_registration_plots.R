load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("parsed_modeling_results_k",k,".Rdata")))

t0_contrasts_cell<-modeling_results$enrichment
rownames(t0_contrasts_cell) = t0_contrasts_cell$ensembl

t0_contrasts_cell <- t0_contrasts_cell[,! names(t0_contrasts_cell)%in%c("ensembl","gene")]


ground_truth <- spatialLIBD::fetch_data("modeling_results")

cor_stats_layer <- layer_stat_cor(
  t0_contrasts_cell,
  modeling_results = ground_truth,
  model_type = "enrichment",
  top_n = 100
)

### heatmap ### here can also use layer_stat_cor_plot() from spatialLIBD

##plot output directory
dir_plots <-
  here::here("plots", "07_spatial_registration")
dir.create(dir_plots, showWarnings = FALSE)

#http://research.libd.org/spatialLIBD/reference/layer_stat_cor_plot.html newer function for plotting



pdf(
  file = here::here(
    "plots",
    "07_spatial_registration",
    paste0(
      "dlpfc_pseudobulked_bayesSpace_vs_mannual_annotations_k",
      k,
      ".pdf"
    )
  ),
  width = 8
)
layer_stat_cor_plot_AS(
  cor_stats_layer,
  max = 1,
  cex = 2.5
)

dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()