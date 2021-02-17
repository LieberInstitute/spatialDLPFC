sample_names <- c("DLPFC_Br2743_ant_manual_alignment", "DLPFC_Br2743_mid_manual_alignment","DLPFC_Br2743_post_manual_alignment","DLPFC_Br3942_ant_manual_alignment","DLPFC_Br3942_mid_manual_alignment","DLPFC_Br3942_post_manual_alignment","DLPFC_Br6423_ant_manual_alignment","DLPFC_Br6423_mid_manual_alignment","DLPFC_Br6423_post_manual_alignment","DLPFC_Br8492_ant_manual_alignment","DLPFC_Br8492_mid_manual_alignment","DLPFC_Br8492_post_manual_alignment")
dir_outputs <- here::here("outputs", "NextSeq")

df_metrics_all <- NULL
for (sample_name in sample_names){
  dir_csv <- file.path(dir_outputs, sample_name, "outs", "metrics_summary.csv")
  df_metrics <- read.csv(dir_csv, header = TRUE)
  df_metrics_all <- rbind(df_metrics_all,df_metrics)
}

# save
sample_metrics <- df_metrics_all
save(sample_metrics, file = here::here("rdata", "spe", "sample_metrics.Rdata"))
write.csv(sample_metrics, file = here::here("rdata", "spe", "sample_metrics.csv"))


## Read in older data
pilot_metrics <- read.table("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/visium_dlpfc_pilot_sample_metrics.tsv", header = TRUE, sep = "\t")

## Reshape the data
pilot_metrics <- as.data.frame(t(pilot_metrics))
colnames(pilot_metrics) <- pilot_metrics[1, ]
pilot_metrics <- pilot_metrics[-1, ]
rownames(pilot_metrics) <- gsub("X", "", rownames(pilot_metrics))

## Make it usable in R
for(i in c(
    "Estimated.Number.of.Spots",
    "Mean.Reads.per.Spot",
    "Median.Genes.per.Spot",
    "Number.of.Reads",
    "Total.Genes.Detected",
    "Median.UMI.Counts.per.Spot",
    "Mean.Cells.Per.Spot",
    "Proportion.0.Cells.Per.Spot",
    "Proportion.1.Cell.Per.Spot",
    "Age.Death"
)) {
    pilot_metrics[[i]] <- as.numeric(gsub(",", "", pilot_metrics[[i]]))
}
for(i in c(
    "Valid.Barcodes",
    "Sequencing.Saturation",
    "Q30.Bases.in.Barcode",
    "Q30.Bases.in.RNA.Read",
    "Q30.Bases.in.Sample.Index",
    "Q30.Bases.in.UMI",
    "Reads.Mapped.to.Genome",
    "Reads.Mapped.Confidently.to.Genome",
    "Reads.Mapped.Confidently.to.Intergenic.Regions",
    "Reads.Mapped.Confidently.to.Intronic.Regions",
    "Reads.Mapped.Confidently.to.Exonic.Regions",
    "Reads.Mapped.Confidently.to.Transcriptome",
    "Reads.Mapped.Antisense.to.Gene",
    "Fraction.Reads.in.Spots"
)) {
    pilot_metrics[[i]] <- as.numeric(gsub("%", "", pilot_metrics[[i]]))
}
summary(pilot_metrics)
save(pilot_metrics, file = here::here("rdata", "spe", "pilot_metrics.Rdata"))
