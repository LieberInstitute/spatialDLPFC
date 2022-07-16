library(ggplot2)
library(here)
library(spatialLIBD)

load(file = here::here("processed-data", "rdata", "spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

# uncorrected umap
pdf(file = here::here("plots", "01_build_spe", "sfigu_harmony_UMAP_sample_id.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
    geom_point(size = 0.3, alpha = 0.5) +
    labs(color = "sample_id") +
    theme_bw() +
    ggtitle("UMAP of Uncorrected Data") +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    scale_colour_discrete(name = "Sample ID") +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 18), plot.title = element_text(size = 20), legend.text = element_text(size = 15))
dev.off()


# corrected umpas
pdf(file = here::here("plots", "01_build_spe", "sfigu_harmony_UMAP_harmony_sample_id.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
    geom_point(size = 0.3, alpha = 0.5) +
    labs(color = "sample_id") +
    theme_bw() +
    ggtitle("UMAP of Batch-Corrected Data") +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    scale_colour_discrete(name = "Sample ID") +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 18), plot.title = element_text(size = 20), legend.text = element_text(size = 15))
dev.off()
