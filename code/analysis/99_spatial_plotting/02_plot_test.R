library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <-
  here(
    "plots",
    "99_spatial_plotting",
    "02_test_plots"
  )
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

sample_lims <- read.csv(here("processed-data","rdata", "spe","99_spatial_plotting","sample_xy_limits.csv"))

# Load full data
load(
  here(
    "processed-data",
    "rdata",
    "spe",
    "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
  )
)

rowData(spe)$gene_search <- rowData(spe)$gene_name
spe$bayesSpace_harmony_9 <- as.factor(spe$bayesSpace_harmony_9)
spe$bayesSpace_harmony_16 <- as.factor(spe$bayesSpace_harmony_16)

###
sl <- sample_lims[1,]

vis_gene_test <- vis_gene(
  spe = spe,
  point_size = 1.2,
  sampleid = "Br2720_ant",
  geneid = "CLDN5"
) + 
  coord_fixed() +
  geom_hline(yintercept = c(sl$x_min, sl$x_max)/100)


ggsave(vis_gene_test, filename = here(plot_dir, "vis_gene_ggsave.png"))
ggsave(vis_gene_test, filename = here(plot_dir, "vis_gene_ggsave.pdf"))


img <- SpatialExperiment::imgRaster(spe, sample_id = "Br2720_ant", image_id = "lowres")
dim(img)

# sample_id x_min x_max y_min y_max
# 1   Br2720_ant 15782 33439  2876 21502
