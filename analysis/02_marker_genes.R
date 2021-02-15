## Automatically style the code in this script:
# styler::style_file("02_marker_genes.R", transformers = biocthis::bioc_style())

## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")
library("rtracklayer")

## vis
library("spatialLIBD")

## analysis
library("scran")
library("scater")
library("BiocParallel")

## are they needed?
# library("tidyverse")
# library("ggplot2")


# remove spots w/ zero UMIs
load(file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_combined.rda")
remove <- which(colData(sce)$sum_umi == 0)
sce <- sce[, -remove]

pdf(
    "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/hist_sum_umi.pdf",
    useDingbats = FALSE
)
hist((colData(sce)$sum_umi), breaks = 200)
dev.off()

# quality control (scran) start here 1/25/21
qcstats <- perCellQCMetrics(sce)
qcfilter <- quickPerCellQC(qcstats)
colSums(as.matrix(qcfilter))

sce$scran_discard <-
    factor(qcfilter$discard, levels = c("TRUE", "FALSE")) # make boxplot of the sum_umis for the  scran_discard=TRUE spots
sce$scran_low_lib_size <-
    factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
sce$low_n_features <-
    factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))

summary(sce$scran_discard)
# TRUE FALSE
# 3055 46944

dim(assay(sce))
# [1] 36601 49999

# add reduced dimensions to sce, takes ~8min, eliminate spots with less than 10 UMIs before clusterings
sce_nonzero <- sce[, sce$sum_umi > 0]
set.seed(20191112)
Sys.time()
clusters_nonzero <- quickCluster(
    sce_nonzero,
    BPPARAM = MulticoreParam(4),
    block = sce_nonzero$sample_name,
    block.BPPARAM = MulticoreParam(4)
)
Sys.time()
save(clusters_nonzero, file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/clusters_nonzero.rda")
save(sce_nonzero, file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_nonzero.rda")
load(file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/clusters_nonzero.rda")

sce_more_than_10_umis <- sce[, sce$sum_umi > 10]
set.seed(20191112)
Sys.time()
clusters_more_than_10_umis <- quickCluster(
    sce_more_than_10_umis,
    BPPARAM = MulticoreParam(4),
    block = sce_more_than_10_umis$sample_name,
    block.BPPARAM = MulticoreParam(4)
)
Sys.time()
save(clusters_more_than_10_umis, file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/clusters_more_than_10_umis.rda")
save(sce_more_than_10_umis, file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_more_than_10_umis.rda")
load(file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_more_than_10_umis.rda")
load(file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/clusters_sce_more_than_10_umis.rda")

sce_more_than_100_umis <- sce[, sce$sum_umi > 100]
set.seed(20191112)
Sys.time()
clusters_more_than_100_umis <- quickCluster(
    sce_more_than_100_umis,
    BPPARAM = MulticoreParam(4),
    block = sce_more_than_100_umis$sample_name,
    block.BPPARAM = MulticoreParam(4)
)
Sys.time()
save(clusters_more_than_100_umis, file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/clusters_more_than_100_umis.rda")
save(sce_more_than_10_umis, file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_more_than_100_umis.rda")
load(file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_more_than_100_umis.rda")

Sys.time()
sce_more_than_10_umis <-
    computeSumFactors(sce_more_than_10_umis,
        clusters = clusters_more_than_10_umis,
        BPPARAM = MulticoreParam(4)
    )
Sys.time()

Sys.time()
sce_more_than_100_umis <-
    computeSumFactors(sce_more_than_100_umis,
        clusters = clusters_more_than_100_umis,
        BPPARAM = MulticoreParam(4)
    )
Sys.time()

sce_nonzero <- sce[, sce$sum_umi > 0]
options(error = recover)
sce_test <-
    computeSumFactors(sce_nonzero[, sce_nonzero$sum_umi > 10], clusters = clusters[sce_nonzero$sum_umi > 10], BPPARAM = MulticoreParam(4))

save(sce, file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_filtered_combined.rda")
load(file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_filtered_combined.rda")

summary(sizeFactors(sce))

sce <- logNormCounts(sce)

## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(sce,
    block = sce$sample_name,
    BPPARAM = MulticoreParam(4)
)
Sys.time()

pdf(
    "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/modelGeneVar.pdf",
    useDingbats = FALSE
)
mapply(function(block, blockname) {
    plot(
        block$mean,
        block$total,
        xlab = "Mean log-expression",
        ylab = "Variance",
        main = blockname
    )
    curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
    )
}, dec$per.block, names(dec$per.block))
dev.off()

top.hvgs <- getTopHVGs(dec, prop = 0.1)
length(top.hvgs)
# 2017

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
# 16831

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
# 15855

set.seed(20191112)
Sys.time()
sce <- runPCA(sce, subset_row = top.hvgs)
Sys.time()

reducedDimNames(sce)

summary(apply(reducedDim(sce, "PCA"), 2, sd))

summary(colMeans(reducedDim(sce, "PCA")))

Sys.time()
set.seed(20191206)
sce <-
    runTSNE(sce,
        dimred = "PCA",
        name = "TSNE_perplexity50",
        perplexity = 50
    )
Sys.time()

Sys.time()
set.seed(20191206)
sce <-
    runTSNE(sce,
        dimred = "PCA",
        name = "TSNE_perplexity5",
        perplexity = 5
    )
Sys.time()

Sys.time()
set.seed(20191206)
sce <-
    runTSNE(sce,
        dimred = "PCA",
        name = "TSNE_perplexity20",
        perplexity = 20
    )
Sys.time()

Sys.time()
set.seed(20191206)
sce <-
    runTSNE(sce,
        dimred = "PCA",
        name = "TSNE_perplexity80",
        perplexity = 80
    )
Sys.time()

Sys.time()
set.seed(20191206)
sce <- runUMAP(sce, dimred = "PCA", name = "UMAP_neighbors15")
Sys.time()

# stopped here and saved object on 2/3/21



spatialLIBD::check_sce(sce)


# save
save(sce, file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_filtered_combined_reduced_dim.rda")



#### plot log UMIs for scran_discard==TRUE
load(file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/sce_combined.rda")
# remove <- which(colData(sce)$scran_discard == TRUE)
# sce_discard <- sce[,remove]

x <-
    which(colData(sce)$sample_name == "DLPFC_Br2743_ant_manual_alignment")
sce_x <- sce[, x]
remove <- which(colData(sce_x)$scran_discard == "FALSE")
sce_x$sum_umi[remove] <- NA
sce_x$height <- 600
sce_x$width <- 504

colData(sce)$height <- NA
colData(sce)$width <- NA
for (i in seq_along(sample_names)){
  x = which(colData(sce)$sample_name == sample_names[i])
  colData(sce)$height[x] <- metadata(sce)$image$height[which(metadata(sce)$image$sample == sample_names[i])]
  colData(sce)$width[x] <- metadata(sce)$image$width[which(metadata(sce)$image$sample == sample_names[i])]
}

pdf(
    "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/scran_discard.pdf",
    useDingbats = FALSE
) # make pdf larger
sce_image_grid_gene(
    sce,
    geneid = "sum_umi",
    spatial = TRUE,
    return_plots = TRUE,
    minCount = -1
)
dev.off()

sce$discard <- qcfilter$discard

plots_discard <-
    sce_image_clus(sce,
        "DLPFC_Br2743_ant_manual_alignment",
        "discard",
        colors = c("light blue", "red")
    )
pdf(
    "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/analysis/scran_discard.pdf",
    height = 24,
    width = 36
)
plot_grid(plotlist = plots_discard)
dev.off()


for (i in seq_along(sample_names)) {
    # select sample
    sce <- sce_list[[i]]

    # identify mitochondrial genes
    is_mito <- grepl("(^MT-)|(^mt-)", rowData(sce)$gene_name)
    table(is_mito)
    rowData(sce)$gene_name[is_mito]

    # calculate QC metrics using scater package
    sce <- addPerCellQC(sce, subsets = list(mito = is_mito))

    colData(sce)

    # store
    sce_list[[i]] <- sce
}

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
