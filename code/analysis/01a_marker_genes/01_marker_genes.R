## This script requires R 4.1
# module load conda_R/devel

## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")

## analysis
library("scran")

## vis
library("spatialLIBD")
library("RColorBrewer")

# set up colors
nb.cols <- 25
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

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


for (i in human_markers_search) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "01a_marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
        assayname = "counts"
    )
}

#
library(ggplot2)
sample_order <- unique(spe$sample_id) # may want to update sample_order
plots_histology <-
    lapply(sample_order, function(sampleid) {
        message(Sys.time(), " ", sampleid)
        spe_sub <- spe[, spe$sample_id == sampleid]
        sample_df <- as.data.frame(colData(spe_sub), optional = TRUE)

        ## Obtain the histology image
        img <- SpatialExperiment::imgRaster(spe_sub)

        ## Transform to a rasterGrob object
        grob <- grid::rasterGrob(img, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))

        ## Make a plot using geom_spatial
        p <- ggplot2::ggplot(
            sample_df,
            ggplot2::aes(
                x = pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe_sub),
                y = pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe_sub),
            )
        ) +
            geom_spatial(
                data = tibble::tibble(grob = list(grob)),
                ggplot2::aes(grob = grob),
                x = 0.5,
                y = 0.5
            ) +
            xlab("") +
            ylab("") #+
        # labs(title = sampleid,caption = NULL)
    })
names(plots_histology) <- sample_order

# plot MOBP
plots_mobp <- vis_grid_gene(
    spe,
    geneid = rowRanges(spe)$gene_search[rowRanges(spe)$gene_name == "MOBP"],
    spatial = FALSE,
    assayname = "logcounts",
    minCount = -1,
    image_id = "lowres",
    alpha = 1,
    point_size = 2.0,
    cont_colors = rev(viridisLite::viridis(21)),
    return_plots = TRUE,
    sample_order = sample_order
)
plots_mobp <- lapply(sample_order, function(sampleid) {
    p <- plots_mobp[[sampleid]]
    # p + labs(title = sampleid,caption = NULL)
    p +
        labs(title = NULL, caption = NULL) +
        theme(legend.key.size = unit(1.5, "cm"), legend.text = element_text(size = 30))
})
names(plots_mobp) <- sample_order

# plot SNAP25
plots_snap25 <- vis_grid_gene(
    spe,
    geneid = rowRanges(spe)$gene_search[rowRanges(spe)$gene_name == "SNAP25"],
    spatial = FALSE,
    assayname = "logcounts",
    minCount = -1,
    image_id = "lowres",
    alpha = 1,
    point_size = 2.0,
    cont_colors = rev(viridisLite::viridis(21)),
    return_plots = TRUE,
    sample_order = sample_order
)
plots_snap25 <- lapply(sample_order, function(sampleid) {
    p <- plots_snap25[[sampleid]]
    # p + labs(title = sampleid,caption = NULL)
    p +
        labs(title = NULL, caption = NULL) +
        theme(legend.key.size = unit(1.5, "cm"), legend.text = element_text(size = 30))
})
names(plots_snap25) <- sample_order

# plot PCP4
plots_pcp4 <- vis_grid_gene(
    spe,
    geneid = rowRanges(spe)$gene_search[rowRanges(spe)$gene_name == "PCP4"],
    spatial = FALSE,
    assayname = "logcounts",
    minCount = -1,
    image_id = "lowres",
    alpha = 1,
    point_size = 2.0,
    cont_colors = rev(viridisLite::viridis(21)),
    return_plots = TRUE,
    sample_order = sample_order
)
plots_pcp4 <- lapply(sample_order, function(sampleid) {
    p <- plots_pcp4[[sampleid]]
    # p + labs(title = sampleid,caption = NULL)
    p +
        labs(title = NULL, caption = NULL) +
        theme(legend.key.size = unit(1.5, "cm"), legend.text = element_text(size = 30))
})
names(plots_pcp4) <- sample_order

sample_info <- data.frame(
    sample_id = c(
        "DLPFC_Br2743_ant_manual_alignment",
        "DLPFC_Br2743_mid_manual_alignment_extra_reads",
        "DLPFC_Br2743_post_manual_alignment",
        "DLPFC_Br3942_ant_manual_alignment",
        "DLPFC_Br3942_mid_manual_alignment",
        "DLPFC_Br3942_post_manual_alignment",
        "DLPFC_Br6423_ant_manual_alignment_extra_reads",
        "DLPFC_Br6423_mid_manual_alignment",
        "DLPFC_Br6423_post_extra_reads",
        "DLPFC_Br8492_ant_manual_alignment",
        "DLPFC_Br8492_mid_manual_alignment_extra_reads",
        "DLPFC_Br8492_post_manual_alignment",
        "DLPFC_Br2720_ant_2",
        "DLPFC_Br2720_mid_manual_alignment",
        "DLPFC_Br2720_post_extra_reads",
        "DLPFC_Br6432_ant_2",
        "DLPFC_Br6432_mid_manual_alignment",
        "DLPFC_Br6432_post_manual_alignment",
        "DLPFC_Br6471_ant_manual_alignment_all",
        "DLPFC_Br6471_mid_manual_alignment_all",
        "DLPFC_Br6471_post_manual_alignment_all",
        "DLPFC_Br6522_ant_manual_alignment_all",
        "DLPFC_Br6522_mid_manual_alignment_all",
        "DLPFC_Br6522_post_manual_alignment_all",
        "DLPFC_Br8325_ant_manual_alignment_all",
        "DLPFC_Br8325_mid_2",
        "DLPFC_Br8325_post_manual_alignment_all",
        "DLPFC_Br8667_ant_extra_reads",
        "DLPFC_Br8667_mid_manual_alignment_all",
        "DLPFC_Br8667_post_manual_alignment_all"
    ),
    subjects = c(rep(
        c("Br2743", "Br3942", "Br6423", "Br8492", "Br2720", "Br6432", "Br6471", "Br6522", "Br8325", "Br8667"),
        each = 3
    )),
    regions = c(rep(
        c("anterior", "middle", "posterior"),
        10
    )),
    sex = c(rep(
        c("M", "M", "M", "F", "F", "M", "M", "M", "F", "F"),
        each = 3
    )),
    age = c(rep(
        c(61.54, 47.53, 51.73, 53.40, 48.22, 48.88, 55.46, 33.39, 57.62, 37.33),
        each = 3
    )),
    diagnosis = "Control"
)

# clean up sample_id
sample_info$sample_id <- gsub("_all|_extra_reads|DLPFC_|_manual_alignment|_2", "", basename(sample_info$sample_id))
sample_info$row <- seq_len(nrow(sample_info))
donor_order <- unique(sample_info$subjects)

lapply(unique(sample_info$regions), function(region) {
    message(Sys.time(), " ", region)
    region_subset <- subset(sample_info, regions == region)
    i <- match(donor_order, region_subset$subjects)
    plots_list <- c(
        plots_histology[region_subset$row[i]],
        plots_snap25[region_subset$row[i]],
        plots_mobp[region_subset$row[i]],
        plots_pcp4[region_subset$row[i]]
    )
    pdf(file = here::here("plots", "01a_marker_genes", paste0("vis_genes_known_markers_sfig_", region, ".pdf")), height = 8 * 10, width = 8 * 4)
    print(cowplot::plot_grid(plotlist = plots_list, ncol = 4, nrow = 10, byrow = FALSE))
    dev.off()
    return(NULL)
})

# make known marker genes plot
regions <- unique(spe$region)
for (i in length(regions)) {
    known_markers <-
        c(
            "SNAP25",
            "MOBP",
            "PCP4"
        )
    for (j in length(known_markers)) {
        p <- vis_gene( # returns ggplot2 object
            spe,
            sampleid,
            geneid = known_markers[i],
            spatial = TRUE,
            assayname = "logcounts",
            minCount = 0,
            viridis = TRUE,
            image_id = "lowres",
            alpha = 1,
            cont_colors = if (viridis) {
                viridisLite::viridis(21)
            } else {
                c(
                    "aquamarine4",
                    "springgreen", "goldenrod", "red"
                )
            },
            point_size = 1.25,
            ...
        )
        # append p to plots
    }

    pdf(file = here::here("plots", "01a_marker_genes", paste0("vis_genes_known_markers_sfig_", regions[i], ".pdf")), height = 24, width = 36)
    print(cowplot::plot_grid(plotlist = plots))
    dev.off()
}

# new_markers_search <- rowData(spe)$gene_search[match(new_markers, rowData(spe)$gene_name)]
new_markers <-
    c(
        "CLDN5",
        "C1QL2",
        "APOE",
        "MSX1",
        "SPARC"
    )
new_markers <- rowData(spe)$gene_search[match(new_markers, rowData(spe)$gene_name)]
samples <- unique(spe$sample_id)


for (i in new_markers) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "01a_marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
        assayname = "counts"
    )
}

for (i in length(new_markers)) {
    pdf(file = here::here("plots", "01a_marker_genes", paste0("vis_genes_new_markers_", ".pdf")))
    for (j in length(samples)) {
        p <- vis_gene( # returns ggplot2 object
            spe,
            sampleid = samples[j],
            geneid = new_markers[i],
            spatial = TRUE,
            assayname = "logcounts",
            minCount = 0,
            viridis = TRUE,
            image_id = "lowres",
            alpha = 1,
            point_size = 1.25
        )
        # append p to plots
        print(p)
    }
    dev.off()
}


pdf(file = here::here("plots", "01a_marker_genes", "vis_genes_new_markers_CLDN5.pdf"))
p <- vis_gene( # returns ggplot2 object
    spe,
    sampleid = "Br8667_mid",
    geneid = "CLDN5; ENSG00000184113",
    spatial = TRUE,
    assayname = "logcounts",
    minCount = 0,
    viridis = TRUE,
    image_id = "lowres",
    alpha = 1,
    point_size = 1.25
)
# append p to plots
print(p)
dev.off()

pdf(file = here::here("plots", "01a_marker_genes", "vis_genes_new_markers_MSX1.pdf"))
p <- vis_gene( # returns ggplot2 object
    spe,
    sampleid = "Br8667_mid",
    geneid = "MSX1; ENSG00000163132",
    spatial = TRUE,
    assayname = "logcounts",
    minCount = 0,
    viridis = TRUE,
    image_id = "lowres",
    alpha = 1,
    point_size = 1.25
)
# append p to plots
print(p)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
