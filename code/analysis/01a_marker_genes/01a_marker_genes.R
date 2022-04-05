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

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()