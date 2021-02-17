sample_names <- c("DLPFC_Br2743_ant_manual_alignment", "DLPFC_Br2743_mid_manual_alignment","DLPFC_Br2743_post_manual_alignment","DLPFC_Br3942_ant_manual_alignment","DLPFC_Br3942_mid_manual_alignment","DLPFC_Br3942_post_manual_alignment","DLPFC_Br6423_ant_manual_alignment","DLPFC_Br6423_mid_manual_alignment","DLPFC_Br6423_post_manual_alignment","DLPFC_Br8492_ant_manual_alignment","DLPFC_Br8492_mid_manual_alignment","DLPFC_Br8492_post_manual_alignment")
dir_outputs <- here::here("outputs", "NextSeq")

df_metrics_all <- NULL
for (sample_name in sample_names){
  dir_csv <- file.path(dir_outputs, sample_name, "outs", "metrics_summary.csv")
  df_metrics <- read.csv(dir_csv, header = TRUE)
  df_metrics_all <- rbind(df_metrics_all,df_metrics)
}

# save
save(df_metrics_all, file = here::here("rdata", "spe", "sample_metrics.Rdata"))
