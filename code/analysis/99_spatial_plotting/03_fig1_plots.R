library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")
library("patchwork")

plot_dir <-
  here(
    "plots",
    "99_spatial_plotting",
    "03_fig1_plots"
  )
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}


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
spe$bayesSpace_harmony_7 <- as.factor(spe$bayesSpace_harmony_7)
spe$bayesSpace_harmony_9 <- as.factor(spe$bayesSpace_harmony_9)
spe$bayesSpace_harmony_16 <- as.factor(spe$bayesSpace_harmony_16)

#### k7 dotplot ####
## prep colors
k_colors <- Polychrome::palette36.colors(28)
names(k_colors) <- c(1:28)

vis_clust_k7 <- vis_clus(
  spe = spe,
  viridis = FALSE,
  point_size = 2,
  colors = k_colors,
  sampleid = "Br8667_mid",
  clustervar = "bayesSpace_harmony_7"
) + labs(title = NULL)+ 
  theme(legend.position = c(0.9, 0.8) # top right
  )

ggsave(vis_clust_k7, filename = here(plot_dir, "vis_clust_Br8667_mid_k7.png"))
ggsave(vis_clust_k7, filename = here(plot_dir, "vis_clust_Br8667_mid_k7.pdf"))

#### plot example genes  ####
example_genes <- c("SNAP25","MBP","PCP4")
names(example_genes) <- example_genes
my_text_size <- 25

gene_plots <- map(example_genes, function(gene){
  s <- "Br8667_mid"
  vis_gene_plot <- vis_gene(
    spe = spe,
    viridis = FALSE,
    point_size = 2,
    cont_colors = viridisLite::plasma(21),
    # frame_lim_df = frame_lims,
    sampleid = s,
    geneid = gene
  ) +
    labs(title = gene) + 
    theme(plot.title = element_text(face = "italic", hjust = 0.5),
          legend.position = c(0.9, 0.8), # top right
          text = element_text(size = my_text_size)
          )
  
  ggsave(vis_gene_plot, filename = here(plot_dir, paste0("vis_gene_",s,"-",gene,".png")))
  ggsave(vis_gene_plot, filename = here(plot_dir, paste0("vis_gene_",s,"-",gene,".pdf")))
  return(vis_gene_plot)
})

## patchwork together

patchwork_genes <- gene_plots$SNAP25 + gene_plots$MBP + gene_plots$PCP4
ggsave(patchwork_genes, filename = here(plot_dir, "vis_gene_Br8667_mid-ALL.png"),width = 18)
ggsave(patchwork_genes, filename = here(plot_dir, "vis_gene_Br8667_mid-ALL.pdf"),width = 18)

# ## try custom functions
# source("vis_gene_crop.R")
# 
# frame_lims <- read.csv(here("processed-data","rdata", "spe","99_spatial_plotting","frame_limits.csv"))
# frame_edge_lims <- read.csv(here("processed-data","rdata", "spe","99_spatial_plotting","frame_edge_limits.csv"))
# 
# walk(c("SNAP25","MBP","PCP4"), function(gene){
#   s <- "Br8667_mid"
#   vis_gene_crop_plot <- vis_gene_crop(
#     spe = spe,
#     point_size = 2.2,
#     frame_lim_df = frame_edge_lims,
#     sampleid = s,
#     geneid = gene,
#     legend_overlap = TRUE
#   ) 
#   
#   ggsave(vis_gene_crop_plot, filename = here(plot_dir, paste0("vis_gene_crop_",s,"-",gene,".png")))
#   ggsave(vis_gene_crop_plot, filename = here(plot_dir, paste0("vis_gene_crop_",s,"-",gene,".pdf")))
#   
# })


