## This script requires R 4.2.x with spatialLIBD version 1.11.4
# module load conda_R/4.2.x

## utils
library("here")
library("sessioninfo")

## vis
library("spatialLIBD")
library("ggplot2")

## Load SPE data
load(here::here("processed-data", "rdata", "spe", "01_build_spe", "spe_filtered_final.Rdata"), verbose = TRUE)

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

## Make a grid plot for each marker
for (i in human_markers_search) {
    message(Sys.time(), " processing gene ", i)
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "01a_marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
        cont_colors = viridisLite::plasma(21)
    )
}


## Plot some new marker genes
new_markers <-
    c(
        "CLDN5",
        "C1QL2",
        "APOE",
        "MSX1",
        "SPARC"
    )
new_markers <- rowData(spe)$gene_search[match(new_markers, rowData(spe)$gene_name)]
for (i in new_markers) {
    message(Sys.time(), " processing gene ", i)
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "01a_marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
        cont_colors = viridisLite::plasma(21)
    )
}

## Individual plots for Br8667_mid on these new markers
pdf(file = here::here("plots", "01a_marker_genes", "vis_genes_new_markers_CLDN5.pdf"))
p <- vis_gene(
    spe,
    sampleid = "Br8667_mid",
    geneid = "CLDN5; ENSG00000184113",
    cont_colors = viridisLite::plasma(21)
)
print(p)
dev.off()

pdf(file = here::here("plots", "01a_marker_genes", "vis_genes_new_markers_MSX1.pdf"))
p <- vis_gene(
    spe,
    sampleid = "Br8667_mid",
    geneid = "MSX1; ENSG00000163132",
    cont_colors = viridisLite::plasma(21)
)
print(p)
dev.off()




## Get histology images by taking advantage of setting alpha = 0
plots_histology <- vis_grid_gene(
    spe,
    clustervar = "10x_graphclust",
    return_plots = TRUE,
    alpha = 0
)
plots_histology <- lapply(plots_histology, function(p) {
    p + ggtitle("")
})

## Obtain plots for a few key genes
key_genes <- c("MOBP", "SNAP25", "PCP4")
names(key_genes) <- key_genes
plots_genes <- lapply(key_genes, function(g) {
    message(Sys.time(), " processing gene ", g)

    ## Obtain the spatial plots
    p_list <- vis_grid_gene(
        spe,
        geneid = rowRanges(spe)$gene_search[rowRanges(spe)$gene_name == g],
        spatial = FALSE,
        assayname = "logcounts",
        cont_colors = viridisLite::plasma(21),
        return_plots = TRUE
    )

    ## Remove the title and the legend title, plus make the key size larger
    lapply(p_list, function(p) {
        p +
            labs(title = NULL) +
            theme(legend.key.size = unit(1.5, "cm"),
                legend.title = element_blank())
    })
})

## Obtain the sample information from the unique IDs
sample_info <- data.frame(
    sample_id = unique(spe$sample_id)
)
sample_info$subjects <- gsub("_.*", "", sample_info$sample_id)
sample_info$positions <- c("ant" = "anterior", "mid" = "middle", "post" = "posterior")[gsub(".*_", "", sample_info$sample_id)]
sample_info$row <- seq_len(nrow(sample_info))
donor_order <- unique(sample_info$subjects)

lapply(unique(sample_info$positions), function(position) {
    message(Sys.time(), " processing position ", position)
    position_subset <- subset(sample_info, positions == position)
    i <- match(donor_order, position_subset$subjects)
    plots_list <- c(
        plots_histology[position_subset$row[i]],
        plots_genes$SNAP25[position_subset$row[i]],
        plots_genes$MOBP[position_subset$row[i]],
        plots_genes$PCP4[position_subset$row[i]]
    )
    pdf(file = here::here("plots", "01a_marker_genes", paste0("vis_genes_known_markers_sfig_", position, ".pdf")), height = 8 * 10, width = 8 * 4)
    print(cowplot::plot_grid(plotlist = plots_list, ncol = 4, nrow = 10, byrow = FALSE))
    dev.off()
    return(NULL)
})

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
