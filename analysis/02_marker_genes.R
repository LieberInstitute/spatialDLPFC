## Automatically style the code in this script:
styler::style_file(
    here::here("analysis", "02_marker_genes.R"),
    transformers = biocthis::bioc_style()
)

## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")

## vis
library("spatialLIBD")


## Load SPE raw data
load(here::here("rdata", "spe", "spe_raw.Rdata"), verbose = TRUE)

## Filter down to spots in tissue
spe <- spe_raw[, which(inTissue(spe_raw))]

## Find marker genes
human_markers <-
    c(
        "SNAP25",
        "MBP",
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
human_markers_search <- rowData(spe)$gene_search[match(human_markers, rowData(spe)$gene_id)]

## Create plots directory
dir.create(here::here("plots", "human_markers"), showWarnings = FALSE)
for (i in human_markers_search) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "human_markers", paste0(gsub("; ", "_", i), ".pdf")),
        assayname = "counts"
    )
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
