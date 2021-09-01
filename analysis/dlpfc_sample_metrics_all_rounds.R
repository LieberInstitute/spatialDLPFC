## Automatically style the code in this script:
styler::style_file(here::here("analysis", "dlpfc_sample_metrics_round2_3.R"),
                   transformers = biocthis::bioc_style()
)


library("ggplot2")
library("sessioninfo")


sample_names <-
  c(
    "DLPFC_Br2743_ant_manual_alignment",
    "DLPFC_Br2743_mid_manual_alignment_extra_reads",
    "DLPFC_Br2743_post_manual_alignment",
    "DLPFC_Br3942_ant_manual_alignment",
    "DLPFC_Br3942_mid_manual_alignment",
    "DLPFC_Br3942_post_manual_alignment",
    "DLPFC_Br6423_ant_manual_alignment_extra_reads",
    "DLPFC_Br6423_mid_manual_alignment",
    "DLPFC_Br6423_post_manual_alignment",
    "DLPFC_Br8492_ant_manual_alignment",
    "DLPFC_Br8492_mid_manual_alignment_extra_reads",
    "DLPFC_Br8492_post_manual_alignment",
    "Round2/DLPFC_Br2720_ant_manual_alignment",
    "Round2/DLPFC_Br2720_mid_manual_alignment",
    "Round2/DLPFC_Br2720_post_manual_alignment",
    "Round2/DLPFC_Br6432_ant_manual_alignment",
    "Round2/DLPFC_Br6432_mid_manual_alignment",
    "Round2/DLPFC_Br6432_post_manual_alignment",
    "Round2/DLPFC_Br6471_ant_manual_alignment_all",
    "Round2/DLPFC_Br6471_mid_manual_alignment_all",
    "Round2/DLPFC_Br6471_post_manual_alignment_all",
    "Round2/DLPFC_Br6522_ant_manual_alignment_all",
    "Round2/DLPFC_Br6522_mid_manual_alignment_all",
    "Round2/DLPFC_Br6522_post_manual_alignment_all",
    "Round2/DLPFC_Br8325_ant_manual_alignment_all",
    "Round2/DLPFC_Br8325_mid_manual_alignment_all",
    "Round2/DLPFC_Br8325_post_manual_alignment_all",
    "Round2/DLPFC_Br8667_ant_manual_alignment_all",
    "Round2/DLPFC_Br8667_mid_manual_alignment_all",
    "Round2/DLPFC_Br8667_post_manual_alignment_all"
    
  )
dir_outputs <- here::here("outputs", "NextSeq")

df_metrics_all <- NULL
for (sample_name in sample_names) {
  dir_csv <-
    file.path(dir_outputs, sample_name, "outs", "metrics_summary.csv")
  df_metrics <- read.csv(dir_csv, header = TRUE)
  df_metrics_all <- rbind(df_metrics_all, df_metrics)
}


# save
sample_metrics <- df_metrics_all
save(sample_metrics,
     file = here::here("rdata", "spe", "sample_metrics_all.Rdata")
)
write.csv(sample_metrics,
          file = here::here("rdata", "spe", "sample_metrics_all.csv")
)


## Read in older data
# round1_metrics <-
#   read.csv(
#     "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/rdata/spe/sample_metrics.csv",
#     header = TRUE)

# ## Reshape the data
# pilot_metrics <- as.data.frame(t(pilot_metrics))
# colnames(pilot_metrics) <- pilot_metrics[1, ]
# pilot_metrics <- pilot_metrics[-1, ]
# rownames(pilot_metrics) <- gsub("X", "", rownames(pilot_metrics))
# 
# ## Make it usable in R
# for (i in c(
#   "Estimated.Number.of.Spots",
#   "Mean.Reads.per.Spot",
#   "Median.Genes.per.Spot",
#   "Number.of.Reads",
#   "Total.Genes.Detected",
#   "Median.UMI.Counts.per.Spot",
#   "Mean.Cells.Per.Spot",
#   "Proportion.0.Cells.Per.Spot",
#   "Proportion.1.Cell.Per.Spot",
#   "Age.Death"
# )) {
#   pilot_metrics[[i]] <- as.numeric(gsub(",", "", pilot_metrics[[i]]))
# }
# for (i in c(
#   "Valid.Barcodes",
#   "Sequencing.Saturation",
#   "Q30.Bases.in.Barcode",
#   "Q30.Bases.in.RNA.Read",
#   "Q30.Bases.in.Sample.Index",
#   "Q30.Bases.in.UMI",
#   "Reads.Mapped.to.Genome",
#   "Reads.Mapped.Confidently.to.Genome",
#   "Reads.Mapped.Confidently.to.Intergenic.Regions",
#   "Reads.Mapped.Confidently.to.Intronic.Regions",
#   "Reads.Mapped.Confidently.to.Exonic.Regions",
#   "Reads.Mapped.Confidently.to.Transcriptome",
#   "Reads.Mapped.Antisense.to.Gene",
#   "Fraction.Reads.in.Spots"
# )) {
#   pilot_metrics[[i]] <- as.numeric(gsub("%", "", pilot_metrics[[i]]))
# }
# summary(pilot_metrics)
# save(pilot_metrics,
#      file = here::here("rdata", "spe", "pilot_metrics.Rdata")
# )
# 
# write.csv(pilot_metrics,
#           file = here::here("rdata", "spe", "pilot_metrics.csv")
# )


# rownames(sample_metrics) <- gsub("DLPFC_|_manual_alignment", "", sample_metrics$Sample.ID)
# shared_cols <-
#   intersect(colnames(sample_metrics), colnames(round1_metrics))
# shared_cols
#  [1] "Number.of.Reads"                                "Mean.Reads.per.Spot"
#  [3] "Median.Genes.per.Spot"                          "Median.UMI.Counts.per.Spot"
#  [5] "Valid.Barcodes"                                 "Sequencing.Saturation"
#  [7] "Q30.Bases.in.Barcode"                           "Q30.Bases.in.RNA.Read"
#  [9] "Q30.Bases.in.UMI"                               "Reads.Mapped.to.Genome"
# [11] "Reads.Mapped.Confidently.to.Genome"             "Reads.Mapped.Confidently.to.Intergenic.Regions"
# [13] "Reads.Mapped.Confidently.to.Intronic.Regions"   "Reads.Mapped.Confidently.to.Exonic.Regions"
# [15] "Reads.Mapped.Confidently.to.Transcriptome"      "Reads.Mapped.Antisense.to.Gene"
# [17] "Total.Genes.Detected"


## Document missing variables from each table
# colnames(pilot_metrics)[!colnames(pilot_metrics) %in% shared_cols]
# [1] "Estimated.Number.of.Spots"   "Q30.Bases.in.Sample.Index"   "Fraction.Reads.in.Spots"     "Mean.Cells.Per.Spot"
# [5] "Proportion.0.Cells.Per.Spot" "Proportion.1.Cell.Per.Spot"  "Brain.Number"                "Position"
# [9] "Replicate"                   "Age.Death"                   "Sex"                         "Primary.Diagnosis"
# colnames(sample_metrics)[!colnames(sample_metrics) %in% shared_cols]
# [1] "Sample.ID"                            "Number.of.Spots.Under.Tissue"         "Mean.Reads.Under.Tissue.per.Spot"
# [4] "Fraction.of.Spots.Under.Tissue"       "Valid.UMIs"                           "Fraction.Reads.in.Spots.Under.Tissue"

## Combine the shared metrics
# tmp <- sample_metrics[, shared_cols]
# regular_cols <-
#   which(
#     colnames(tmp) %in% c(
#       "Number.of.Reads",
#       "Mean.Reads.per.Spot",
#       "Median.Genes.per.Spot",
#       "Median.UMI.Counts.per.Spot",
#       "Total.Genes.Detected"
#     )
#   )
# tmp[, -regular_cols] <- tmp[, -regular_cols] * 100
# round1_metrics<-round1_metrics[,shared_cols]
# shared_metrics <- rbind(sample_metrics,
#                         round1_metrics)
# shared_metrics$study <- c(rep("new",14),rep("round1",12))
sample_metrics$study <- c(rep("round1",12),rep("round2",6),rep("round3",12))

save(sample_metrics,
     file = here::here("rdata", "spe", "sample_metrics_all.Rdata")
)
write.csv(sample_metrics,
          file = here::here("rdata", "spe", "sample_metrics_all.csv")
)

pdf(
  here::here("plots", "spaceranger_metrics_by_number_of_reads_072821_all.pdf"),
  useDingbats = FALSE,
  width = 10
)
ggplot(
  sample_metrics,
  aes(x = Mean.Reads.per.Spot, y = Number.of.Reads / 1e6)
) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 20) +
  facet_grid(~
               study)

ggplot(
  sample_metrics,
  aes(x = Median.Genes.per.Spot, y = Number.of.Reads / 1e6, color = study)
) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 20)
ggplot(
  sample_metrics,
  aes(x = Total.Genes.Detected, y = Number.of.Reads / 1e6, color = study)
) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 20)
ggplot(
  sample_metrics,
  aes(x = Median.UMI.Counts.per.Spot, y = Number.of.Reads / 1e6, color = study)
) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw(base_size = 20)

dev.off()


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()