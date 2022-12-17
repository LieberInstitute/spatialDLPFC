library("ggplot2")
library("ggpubr")
library("here")
library("sessioninfo")

## Load the ARI results
load(here("processed-data", "rdata", "pilot_dlpfc_data", "05_ARI", "pilot_ari_clustering_across.Rdata"), verbose = TRUE)

pdf(here::here("plots", "05_ARI", "ggboxplot_pilot_data_ARI_clustering_across.pdf"))
ggboxplot(ari.df.long,
    x = "method", y = "ari", color = "general_method",
    palette = viridis(3), add = "jitter", repel = TRUE,
    font.label = list(size = 10), legend = "none", ggtheme = theme_pubr(base_size = 20),
    ylab = "Adjusted Rand Index", xlab = "Clustering Method", size = 1
) +
    font("xy.text", size = 11) +
    font("xlab", size = 16) +
    font("ylab", size = 16)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
