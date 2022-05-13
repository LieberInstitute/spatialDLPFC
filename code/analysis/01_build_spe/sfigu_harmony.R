
#uncorrected umap
pdf(file=here::here("plots", "01_build_spe","UMAP_sample_id.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw() +
  ggtitle("UMAP of Uncorrected Data")
dev.off()


#corrected umpas
pdf(file=here::here("plots", "01_build_spe","UMAP_harmony_sample_id.pdf"))
ggplot(data.frame(reducedDim(spe, "UMAP.HARMONY")),
       aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw() +
  ggtitle("UMAP of Batch-Corrected Data")
dev.off()
