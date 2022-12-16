library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")
library("patchwork")

plot_dir <-
  here(
    "plots",
    "99_spatial_plotting",
    "04_fig2_plots"
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
spe$bayesSpace_harmony_9 <- as.factor(spe$bayesSpace_harmony_9)
spe$bayesSpace_harmony_16 <- as.factor(spe$bayesSpace_harmony_16)
spe$bayesSpace_harmony_28 <- as.factor(spe$bayesSpace_harmony_28)

#### Cluster Plots ####
samples <- c("Br8667_mid", "Br6522_ant", "Br6432_ant")
names(samples) <- samples

## prep colors
k_colors <- Polychrome::palette36.colors(28)
names(k_colors) <- c(1:28)

k_levels <- list(k09 = "bayesSpace_harmony_9", k16 ="bayesSpace_harmony_16", k28 ="bayesSpace_harmony_28")

cluster_plots <- map2(k_levels, names(k_levels), function(k, k_label){
  message(k_label)
  cluster_row_plots <- map(samples, function(s){
    vis_clus_plot <- vis_clus(
      spe = spe,
      viridis = FALSE,
      point_size = 2.2,
      colors = k_colors,
      sampleid = s,
      clustervar = k
    ) +
      theme(legend.position = "None", ## using heat maps label colors
            axis.title.x = element_blank()) 
    
    ## Add sample labels to top row
    if(k_label == "k09") {
      message("add title for k09")
      vis_clus_plot <- vis_clus_plot +
        labs(title = s) +
        theme(plot.title = element_text(hjust = 0.5))
    } else {
      vis_clus_plot <- vis_clus_plot +
        labs(title = "")
      # vis_clus_plot <- vis_clus_plot + theme(title=element_blank())
    }
    ## Add k labels to left 
    if(s == "Br8667_mid") vis_clus_plot <- vis_clus_plot + ylab(k_label)
    
    # ggsave(vis_clus_plot, filename = here(plot_dir, paste0("vis_clus_",s,"-",k,".png")))
    return(vis_clus_plot)
  })
  cluster_row <- Reduce("+", cluster_row_plots)
  # ggsave(cluster_row, filename = here(plot_dir, paste0("vis_clust_",k_label,"_row.png")), width = 18)
  return(cluster_row)
})

ggsave(cluster_plots$k09, filename = ggsave(here(plot_dir, "test09.png")), width = 18)
ggsave(cluster_plots$k16, filename = ggsave(here(plot_dir, "test16.png")), width = 18)

cluster_grid <- Reduce("/", cluster_plots)
ggsave(cluster_grid, filename = here(plot_dir, "vis_grid.png"), width = 18, height = 18)

## plot example genes 
example_genes <- c("SNAP25","MBP","PCP4")
names(example_genes) <- example_genes

gene_plots <- map(example_genes, function(gene){
  s <- "Br8667_mid"
  vis_gene_plot <- vis_gene(
    spe = spe,
    viridis = FALSE,
    point_size = 2.2,
    cont_colors = viridisLite::plasma(21),
    # frame_lim_df = frame_lims,
    sampleid = s,
    geneid = gene
  ) +
    labs(title = gene) + 
    theme(plot.title = element_text(face = "italic"),
          # legend.position = c(0.9, 0.02) # bottom right
          legend.position = c(0.9, 0.8) # top right
          )
  
  ggsave(vis_gene_plot, filename = here(plot_dir, paste0("vis_gene_",s,"-",gene,".png")))
  ggsave(vis_gene_plot, filename = here(plot_dir, paste0("vis_gene_",s,"-",gene,".pdf")))
  return(vis_gene_plot)
})

## patchwork together

patchwork_genes <- gene_plots$SNAP25 + gene_plots$MBP + gene_plots$PCP4
ggsave(patchwork_genes, filename = here(plot_dir, "vis_gene_Br8667_mid-ALL.png"),width = 18)
ggsave(patchwork_genes, filename = here(plot_dir, "vis_gene_Br8667_mid-ALL.pdf"),width = 18)


# })


