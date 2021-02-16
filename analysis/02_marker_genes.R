## Automatically style the code in this script:
# styler::style_file("02_marker_genes.R", transformers = biocthis::bioc_style())

## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")

## vis
library("spatialLIBD")


## Load SPE data

## Filter down to spots in tissue


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

colors <- c("navy", "dodgerblue2")


## Plot marker genes
pdf("DLPFC/marker_genes.pdf", useDingbats = FALSE)

for (i in seq_along(sample_names)) {
    # select sample
    sce <- sce_list[[i]]

    for (j in seq_along(human_markers)) {
        # identify marker gene
        ix_marker <-
            which(toupper(rowData(sce)$gene_name) == toupper(human_markers[j]))
        stopifnot(length(ix_marker) == 1)
        colData(sce)$counts_marker <- counts(sce)[ix_marker, ]


        # plot UMI counts for marker gene

        p <- ggplot(
            as.data.frame(colData(sce)),
            aes(
                x = pxl_row_in_fullres,
                y = pxl_col_in_fullres,
                color = counts_marker
            )
        ) +
            geom_point(size = 1.0) +
            coord_fixed() +
            scale_y_reverse() +
            scale_color_gradient(low = "gray95", high = colors[1]) +
            ggtitle(paste0("UMI counts: ", human_markers[j], ": ", sample_names[i])) +
            labs(color = "counts") +
            theme_bw() +
            theme(
                panel.grid = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank()
            )

        print(p)
    }
}
dev.off()


# human_markers <- c("SNAP25", "MBP","PCP4", "RELN","AQP4","CUX2","CCK","HPCAL1")
#
# colors <- c("navy", "dodgerblue2")
pdf("DLPFC/marker_genes_by_gene.pdf", useDingbats = FALSE)
for (i in seq_along(human_markers)) {
    # select sample


    for (j in seq_along(sample_names)) {
        sce <- sce_list[[j]]

        # identify marker gene
        ix_marker <-
            which(toupper(rowData(sce)$gene_name) == toupper(human_markers[i]))
        stopifnot(length(ix_marker) == 1)
        colData(sce)$counts_marker <- counts(sce)[ix_marker, ]


        # plot UMI counts for marker gene

        p <- ggplot(
            as.data.frame(colData(sce)),
            aes(
                x = pxl_row_in_fullres,
                y = pxl_col_in_fullres,
                color = counts_marker
            )
        ) +
            geom_point(size = 1.0) +
            coord_fixed() +
            scale_y_reverse() +
            scale_color_gradient(low = "gray95", high = colors[1]) +
            ggtitle(paste0("UMI counts: ", human_markers[i], ": ", sample_names[j])) +
            labs(color = "counts") +
            theme_bw() +
            theme(
                panel.grid = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank()
            )

        print(p)
    }
}
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
