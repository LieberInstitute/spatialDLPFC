library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <-
  here(
    "plots",
    "99_spatial_plotting",
    "03_fig1_plots"
  )
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

frame_lims <- read.csv(here("processed-data","rdata", "spe","99_spatial_plotting","frame_limits.csv"))
frame_edge_lims <- read.csv(here("processed-data","rdata", "spe","99_spatial_plotting","frame_edge_limits.csv"))

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


## try custom functions
source("vis_gene_crop.R")

walk(c("SNAP25","MBP","PCP4"), function(gene){
  s <- "Br8667_mid"
  vis_gene_crop_plot <- vis_gene_crop(
    spe = spe,
    point_size = 2.2,
    frame_lim_df = frame_edge_lims,
    sampleid = s,
    geneid = gene,
    legend_overlap = TRUE
  ) 
  
  ggsave(vis_gene_crop_plot, filename = here(plot_dir, paste0("vis_gene_crop_",s,"-",gene,".png")))
  ggsave(vis_gene_crop_plot, filename = here(plot_dir, paste0("vis_gene_crop_",s,"-",gene,".pdf")))
  
})

# walk(c("SNAP25","MBP","PCP4"), function(gene){
#   s <- "Br8667_mid"
#   vis_gene_plot <- vis_gene(
#     spe = spe,
#     point_size = 1.8,
#     # frame_lim_df = frame_lims,
#     sampleid = s,
#     geneid = gene
#   ) +coord_fixed()
#   
#   ggsave(vis_gene_plot, filename = here(plot_dir, paste0("vis_gene_",s,"-",gene,".png")), height = 10, width = 10)
#   
# })

