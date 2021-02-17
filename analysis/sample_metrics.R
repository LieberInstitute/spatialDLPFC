sample_names <- c("DLPFC_Br2743_ant_manual_alignment", "DLPFC_Br2743_mid_manual_alignment","DLPFC_Br2743_post_manual_alignment","DLPFC_Br3942_ant_manual_alignment","DLPFC_Br3942_mid_manual_alignment","DLPFC_Br3942_post_manual_alignment","DLPFC_Br6423_ant_manual_alignment","DLPFC_Br6423_mid_manual_alignment","DLPFC_Br6423_post_manual_alignment","DLPFC_Br8492_ant_manual_alignment","DLPFC_Br8492_mid_manual_alignment","DLPFC_Br8492_post_manual_alignment")
dir_outputs <- "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/outputs/NextSeq"

df_metrics_all <-data.frame(ncol = 12, nrow = 23)
for (i in seq_along(sample_names)){
  sample_name <- sample_names[i]
  dir_csv <- file.path(dir_outputs, sample_name, "outs", "metrics_summary.csv")
  df_metrics <- read.csv(dir_csv, header = TRUE)
  df_metrics <- as.data.frame(df_metrics, header = TRUE)
  df_metrics <- t(df_metrics)
  print(dim(df_metrics))
  df_metrics_all <-cbind(df_metrics_all,df_metrics)
}

save(df_metrics_all, file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sample_metrics.RDS")
