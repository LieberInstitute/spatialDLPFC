## This script requires R 4.1
# module load conda_R/devel

## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")

##analysis
library("scran")

## vis
library("spatialLIBD")
library("RColorBrewer")

#set up colors
nb.cols <- 25
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

## Load SPE data
load(here::here("processed-data", "rdata", "spe","01_build_spe", "spe_filtered_final.Rdata"), verbose = TRUE)

## Find marker genes
human_markers <-
  c(
    "SNAP25",
    "MBP",
    "MOBP",
    "PCP4",
    "RELN",
    "AQP4",
    "CUX2",
    "CCK",
    "HPCAL1",
    "NR4A2",
    "RORB"
  )

## Locate the marker genes
human_markers_search <- rowData(spe)$gene_search[match(human_markers, rowData(spe)$gene_name)]


for (i in human_markers_search) {
  vis_grid_gene(
    spe = spe,
    geneid = i,
    pdf = here::here("plots", "01a_marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
    assayname = "counts"
  )
}


#make known marker genes plot
regions <- unique(spe$region)
for (i in length(regions)){
  known_markers <-
    c(
      "SNAP25",
      "MOBP",
      "PCP4"
    )
  for(j in length(known_markers)){
    p <-vis_gene( #returns ggplot2 object
      spe,
      sampleid,
      geneid = known_markers[i],
      spatial = TRUE,
      assayname = "logcounts",
      minCount = 0,
      viridis = TRUE,
      image_id = "lowres",
      alpha = 1,
      cont_colors = if (viridis) viridisLite::viridis(21) else c("aquamarine4",
                                                                 "springgreen", "goldenrod", "red"),
      point_size = 1.25,
      ...
    )
    #append p to plots
  }
  
  pdf(file = here::here("plots", "01a_marker_genes",paste0("vis_genes_known_markers_sfig_",regions[i],".pdf")), height = 24, width = 36)
  print(cowplot::plot_grid(plotlist = plots))
  dev.off()
}



#use vis_gene
p <-vis_gene(
  spe,
  sampleid,
  geneid = known_markers[i],
  spatial = TRUE,
  assayname = "logcounts",
  minCount = 0,
  viridis = TRUE,
  image_id = "lowres",
  alpha = 1,
  cont_colors = if (viridis) viridisLite::viridis(21) else c("aquamarine4",
                                                             "springgreen", "goldenrod", "red"),
  point_size = 1.25,
  ...
)

pdf(file = here::here("plots", "01a_marker_genes",paste0("vis_genes_known_markers_sfig.pdf")), height = 24, width = 36)
print(cowplot::plot_grid(plotlist = p))
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()