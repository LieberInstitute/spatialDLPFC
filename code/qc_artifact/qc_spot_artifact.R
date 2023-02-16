rm(list = ls())
.libPaths(c(.libPaths(), "/data/apps/extern/anaconda/envs/seurat/4.1.1/lib/R/library", "/home/pravich2/R/x86_64-conda-linux-gnu-library/4.1"))
library("SpatialExperiment")
library("spatialLIBD")
library("RColorBrewer")
library(ggplot2)
library(cowplot)
library("scran") ## requires uwot for UMAP
library("scater")
library("BiocParallel")
library("PCAtools")
library(harmony)
library(dplyr)
# Read in the data
datDir <- "/home/pravich2/scratch16-abattle4/prashanthi/dewrinkler/processed-data/"
Br6522_ant <- read10xVisiumWrapper(
    samples = paste0(datDir, "NextSeq/Round3/DLPFC_Br6522_ant_manual_alignment_all/outs/"),
    sample_id = "Br6522_ant",
    type = "sparse",
    data = "filtered",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    reference_gtf = "/data/abattle4/prashanthi/dewrinkler/data/genes.gtf"
)

Br6522_mid <- read10xVisiumWrapper(
    samples = paste0(datDir, "NextSeq/Round3/DLPFC_Br6522_mid_manual_alignment_all/outs/"),
    sample_id = "Br6522_mid",
    type = "sparse",
    data = "filtered",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    reference_gtf = "/data/abattle4/prashanthi/dewrinkler/data/genes.gtf"
)

Br8667_post <- read10xVisiumWrapper(
    samples = paste0(datDir, "NextSeq/Round3/DLPFC_Br8667_post_manual_alignment_all/outs/"),
    sample_id = "Br8667_post",
    type = "sparse",
    data = "filtered",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    reference_gtf = "/data/abattle4/prashanthi/dewrinkler/data/genes.gtf"
)

Br6522_ant$subject <- "Br6522"
Br6522_ant$diagnosis <- "control"
Br6522_ant$sex <- "M"
Br6522_ant$age <- 33.39
Br6522_ant$region <- "anterior"

Br6522_mid$subject <- "Br6522"
Br6522_mid$diagnosis <- "control"
Br6522_mid$sex <- "M"
Br6522_mid$age <- 33.39
Br6522_mid$region <- "middle"

Br8667_post$subject <- "Br8667"
Br8667_post$diagnosis <- "control"
Br8667_post$sex <- "F"
Br8667_post$age <- 37.33
Br8667_post$region <- "posterior"

## Add information used by spatialLIBD
is_mito <- which(seqnames(Br6522_ant) == "chrM")
Br6522_ant$expr_chrM <- colSums(counts(Br6522_ant)[is_mito, , drop = FALSE])
Br6522_ant$expr_chrM_ratio <- Br6522_ant$expr_chrM / Br6522_ant$sum_umi

is_mito <- which(seqnames(Br6522_mid) == "chrM")
Br6522_mid$expr_chrM <- colSums(counts(Br6522_mid)[is_mito, , drop = FALSE])
Br6522_mid$expr_chrM_ratio <- Br6522_mid$expr_chrM / Br6522_mid$sum_umi

is_mito <- which(seqnames(Br8667_post) == "chrM")
Br8667_post$expr_chrM <- colSums(counts(Br8667_post)[is_mito, , drop = FALSE])
Br8667_post$expr_chrM_ratio <- Br8667_post$expr_chrM / Br8667_post$sum_umi


Br6522_ant_layers <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br6522_ant_layers.csv")
Br6522_mid_layers <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br6522_mid_layers.csv")
Br8667_post_layers <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br8667_post_layers.csv")

Br6522_ant_wrinkles <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br6522_ant_wrinkle.csv")
Br6522_mid_wrinkles <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br6522_mid_wrinkle.csv")
Br8667_post_wrinkles <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br8667_post_wrinkle.csv")

scDir <- "/home/pravich2/scratch16-abattle4/prashanthi/dewrinkler/processed-data/NextSeq/Round3/"
Br6522_ant_spot_counts <- read.csv(paste0(scDir, "DLPFC_Br6522_ant_manual_alignment_all/outs/spatial/tissue_spot_counts_new.csv"))
Br6522_mid_spot_counts <- read.csv(paste0(scDir, "DLPFC_Br6522_mid_manual_alignment_all/outs/spatial/tissue_spot_counts_new.csv"))
Br8667_post_spot_counts <- read.csv(paste0(scDir, "DLPFC_Br8667_post_manual_alignment_all/outs/spatial/tissue_spot_counts_new.csv"))

Br6522_ant_spot_counts <- Br6522_ant_spot_counts[match(colnames(Br6522_ant), Br6522_ant_spot_counts$barcode), ]
Br6522_mid_spot_counts <- Br6522_mid_spot_counts[match(colnames(Br6522_mid), Br6522_mid_spot_counts$barcode), ]
Br8667_post_spot_counts <- Br8667_post_spot_counts[match(colnames(Br8667_post), Br8667_post_spot_counts$barcode), ]

Br6522_ant_spot_counts <- Br6522_ant_spot_counts[, colnames(Br6522_ant_spot_counts) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol")]
Br6522_mid_spot_counts <- Br6522_mid_spot_counts[, colnames(Br6522_mid_spot_counts) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol")]
Br8667_post_spot_counts <- Br8667_post_spot_counts[, colnames(Br8667_post_spot_counts) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol")]

colData(Br6522_ant) <- cbind(colData(Br6522_ant), Br6522_ant_spot_counts)
colData(Br6522_mid) <- cbind(colData(Br6522_mid), Br6522_mid_spot_counts)
colData(Br8667_post) <- cbind(colData(Br8667_post), Br8667_post_spot_counts)

manually_annotate <- function(sobj, layers_df, wrinkles_df, sc_df) {
    spots <- colnames(sobj)
    layers <- c()
    wrinkle <- c()
    for (i in c(1:length(spots))) {
        if (spots[i] %in% layers_df$spot_name) {
            layers[i] <- layers_df$ManualAnnotation[layers_df$spot_name == spots[i]]
        } else {
            layers[i] <- "Unknown"
        }
        if (spots[i] %in% wrinkles_df$spot_name) {
            wrinkle[i] <- wrinkles_df$ManualAnnotation[wrinkles_df$spot_name == spots[i]]
        } else {
            wrinkle[i] <- "None"
        }
    }
    sobj$Layers <- layers
    sobj$Wrinkles <- wrinkle
    sobj
}

Br6522_ant <- manually_annotate(Br6522_ant, Br6522_ant_layers, Br6522_ant_wrinkles, Br6522_ant_spot_counts)
Br6522_mid <- manually_annotate(Br6522_mid, Br6522_mid_layers, Br6522_mid_wrinkles, Br6522_mid_spot_counts)
Br8667_post <- manually_annotate(Br8667_post, Br8667_post_layers, Br8667_post_wrinkles, Br8667_post_spot_counts)

Br6522_ant$is_wrinkle <- !Br6522_ant$Wrinkles == "None"
Br6522_mid$is_wrinkle <- !Br6522_mid$Wrinkles == "None"
Br8667_post$is_wrinkle <- !Br8667_post$Wrinkles == "None"

rownames(Br6522_ant) <- Br6522_ant@rowRanges$gene_name
rownames(Br6522_mid) <- Br6522_mid@rowRanges$gene_name
rownames(Br8667_post) <- Br8667_post@rowRanges$gene_name

library(ggspavis)
plotSpots(Br6522_ant, x_coord = "row", y_coord = "col", annotate = "Layers", y_reverse = FALSE)
plotSpots(Br6522_ant, x_coord = "row", y_coord = "col", annotate = "is_wrinkle", y_reverse = FALSE)

plotSpots(Br6522_mid, x_coord = "row", y_coord = "col", annotate = "Layers", y_reverse = FALSE)
plotSpots(Br6522_mid, x_coord = "row", y_coord = "col", annotate = "is_wrinkle", y_reverse = FALSE)

plotSpots(Br8667_post, x_coord = "row", y_coord = "col", annotate = "Layers", y_reverse = FALSE)
plotSpots(Br8667_post, x_coord = "row", y_coord = "col", annotate = "is_wrinkle", y_reverse = FALSE)

# Remove genes with no expression
no_expr <- which(rowSums(counts(Br6522_ant)) == 0)
length(no_expr)
length(no_expr) / nrow(Br6522_ant) * 100
Br6522_ant <- Br6522_ant[-no_expr, ]
Br6522_ant <- Br6522_ant[, !Br6522_ant@colData$sum_umi == 0]

no_expr <- which(rowSums(counts(Br6522_mid)) == 0)
length(no_expr)
length(no_expr) / nrow(Br6522_mid) * 100
Br6522_mid <- Br6522_mid[-no_expr, ]
Br6522_mid <- Br6522_mid[, !Br6522_mid@colData$sum_umi == 0]

no_expr <- which(rowSums(counts(Br8667_post)) == 0)
length(no_expr)
length(no_expr) / nrow(Br8667_post) * 100
Br8667_post <- Br8667_post[-no_expr, ]
Br8667_post <- Br8667_post[, !Br8667_post@colData$sum_umi == 0]

Br6522_ant_qcstats <- perCellQCMetrics(Br6522_ant, subsets = list(
    Mito = which(seqnames(Br6522_ant) == "chrM")
))
Br6522_ant_qcfilter <- quickPerCellQC(Br6522_ant_qcstats, sub.fields = "subsets_Mito_percent")
colSums(as.matrix(Br6522_ant_qcfilter))
Br6522_ant$scran_discard <-
    factor(Br6522_ant_qcfilter$discard, levels = c("TRUE", "FALSE"))
Br6522_ant$scran_low_lib_size <-
    factor(
        isOutlier(
            Br6522_ant$sum_umi,
            type = "lower",
            log = TRUE
        ),
        levels = c("TRUE", "FALSE")
    )
Br6522_ant$scran_low_n_features <-
    factor(Br6522_ant_qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
Br6522_ant$scran_high_subsets_Mito_percent <-
    factor(Br6522_ant_qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

plotSpots(Br6522_ant, x_coord = "row", y_coord = "col", annotate = "scran_low_lib_size", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Anterior")
plotSpots(Br6522_ant, x_coord = "row", y_coord = "col", annotate = "scran_low_n_features", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Anterior")
plotSpots(Br6522_ant, x_coord = "row", y_coord = "col", annotate = "scran_high_subsets_Mito_percent", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Anterior")
plotSpots(Br6522_ant, x_coord = "row", y_coord = "col", annotate = "scran_discard", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Anterior")

Br6522_mid_qcstats <- perCellQCMetrics(Br6522_mid, subsets = list(
    Mito = which(seqnames(Br6522_mid) == "chrM")
))
Br6522_mid_qcfilter <- quickPerCellQC(Br6522_mid_qcstats, sub.fields = "subsets_Mito_percent")
colSums(as.matrix(Br6522_mid_qcfilter))
Br6522_mid$scran_discard <-
    factor(Br6522_mid_qcfilter$discard, levels = c("TRUE", "FALSE"))
Br6522_mid$scran_low_lib_size <-
    factor(
        isOutlier(
            Br6522_mid$sum_umi,
            type = "lower",
            log = TRUE
        ),
        levels = c("TRUE", "FALSE")
    )
Br6522_mid$scran_low_n_features <-
    factor(Br6522_mid_qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
Br6522_mid$scran_high_subsets_Mito_percent <-
    factor(Br6522_mid_qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

plotSpots(Br6522_mid, x_coord = "row", y_coord = "col", annotate = "scran_low_lib_size", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Middle")
plotSpots(Br6522_mid, x_coord = "row", y_coord = "col", annotate = "scran_low_n_features", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Middle")
plotSpots(Br6522_mid, x_coord = "row", y_coord = "col", annotate = "scran_high_subsets_Mito_percent", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Middle")
plotSpots(Br6522_mid, x_coord = "row", y_coord = "col", annotate = "scran_discard", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Middle")

Br8667_post_qcstats <- perCellQCMetrics(Br8667_post, subsets = list(
    Mito = which(seqnames(Br8667_post) == "chrM")
))
Br8667_post_qcfilter <- quickPerCellQC(Br8667_post_qcstats, sub.fields = "subsets_Mito_percent")
colSums(as.matrix(Br8667_post_qcfilter))
Br8667_post$scran_discard <-
    factor(Br8667_post_qcfilter$discard, levels = c("TRUE", "FALSE"))
Br8667_post$scran_low_lib_size <-
    factor(
        isOutlier(
            Br8667_post$sum_umi,
            type = "lower",
            log = TRUE
        ),
        levels = c("TRUE", "FALSE")
    )
Br8667_post$scran_low_n_features <-
    factor(Br8667_post_qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
Br8667_post$scran_high_subsets_Mito_percent <-
    factor(Br8667_post_qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

plotSpots(Br8667_post, x_coord = "row", y_coord = "col", annotate = "scran_low_lib_size", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Middle")
plotSpots(Br8667_post, x_coord = "row", y_coord = "col", annotate = "scran_low_n_features", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Middle")
plotSpots(Br8667_post, x_coord = "row", y_coord = "col", annotate = "scran_high_subsets_Mito_percent", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Middle")
plotSpots(Br8667_post, x_coord = "row", y_coord = "col", annotate = "scran_discard", y_reverse = FALSE, palette = c("red", "gray")) + ggtitle("Br6522 Middle")

spe <- readRDS("/data/abattle4/prashanthi/dewrinkler/data/spe_filtered_final_with_clusters_and_deconvolution_results.rds")

Br6522_ant_samples <- spe@colData$key[spe@colData$sample_id == "Br6522_ant"]
Br6522_mid_samples <- spe@colData$key[spe@colData$sample_id == "Br6522_mid"]
Br8667_post_samples <- spe@colData$key[spe@colData$sample_id == "Br8667_post"]

Br6522_ant_samples <- gsub("_Br6522_ant", "", Br6522_ant_samples)
Br6522_mid_samples <- gsub("_Br6522_mid", "", Br6522_mid_samples)
Br8667_post_samples <- gsub("_Br8667_post", "", Br8667_post_samples)

Br6522_ant <- Br6522_ant[, colnames(Br6522_ant) %in% Br6522_ant_samples]
Br6522_mid <- Br6522_mid[, colnames(Br6522_mid) %in% Br6522_mid_samples]
Br8667_post <- Br8667_post[, colnames(Br8667_post) %in% Br8667_post_samples]

plotSpots(Br6522_ant, x_coord = "row", y_coord = "col", annotate = "Layers", y_reverse = FALSE)
plotSpots(Br6522_mid, x_coord = "row", y_coord = "col", annotate = "Layers", y_reverse = FALSE)
plotSpots(Br8667_post, x_coord = "row", y_coord = "col", annotate = "Layers", y_reverse = FALSE)


# Filter MALAT1
Br6522_ant <- Br6522_ant[!grepl("MALAT1", rownames(Br6522_ant)), ]
# Filter Mitocondrial
Br6522_ant <- Br6522_ant[!grepl("^MT-", rownames(Br6522_ant)), ]
# Filter Ribosomal
Br6522_ant <- Br6522_ant[!grepl("^RP[SL]", rownames(Br6522_ant)), ]


# Filter MALAT1
Br6522_mid <- Br6522_mid[!grepl("MALAT1", rownames(Br6522_mid)), ]
# Filter Mitocondrial
Br6522_mid <- Br6522_mid[!grepl("^MT-", rownames(Br6522_mid)), ]
# Filter Ribosomal
Br6522_mid <- Br6522_mid[!grepl("^RP[SL]", rownames(Br6522_mid)), ]

# Filter MALAT1
Br8667_post <- Br8667_post[!grepl("MALAT1", rownames(Br8667_post)), ]
# Filter Mitocondrial
Br8667_post <- Br8667_post[!grepl("^MT-", rownames(Br8667_post)), ]
# Filter Ribosomal
Br8667_post <- Br8667_post[!grepl("^RP[SL]", rownames(Br8667_post)), ]

Br6522_ant_normal <- Br6522_ant[, Br6522_ant$Wrinkles == "None"]
Br6522_mid_normal <- Br6522_mid[, Br6522_mid$Wrinkles == "None"]
Br8667_post_normal <- Br8667_post[, Br8667_post$Wrinkles == "None"]

# Normalize the count data
set.seed(030122)
Br6522_ant$scran_quick_cluster <- quickCluster(Br6522_ant)
Br6522_mid$scran_quick_cluster <- quickCluster(Br6522_mid)
Br8667_post$scran_quick_cluster <- quickCluster(Br8667_post)

Br6522_ant <- computeSumFactors(Br6522_ant, clusters = Br6522_ant$scran_quick_cluster)
Br6522_mid <- computeSumFactors(Br6522_mid, clusters = Br6522_mid$scran_quick_cluster)
Br8667_post <- computeSumFactors(Br8667_post, clusters = Br8667_post$scran_quick_cluster)

Br6522_ant <- logNormCounts(Br6522_ant)
Br6522_mid <- logNormCounts(Br6522_mid)
Br8667_post <- logNormCounts(Br8667_post)

Br6522_ant_normal$scran_quick_cluster <- quickCluster(Br6522_ant_normal)
Br6522_mid_normal$scran_quick_cluster <- quickCluster(Br6522_mid_normal)
Br8667_post_normal$scran_quick_cluster <- quickCluster(Br8667_post_normal)

Br6522_ant_normal <- computeSumFactors(Br6522_ant_normal, clusters = Br6522_ant_normal$scran_quick_cluster)
Br6522_mid_normal <- computeSumFactors(Br6522_mid_normal, clusters = Br6522_mid_normal$scran_quick_cluster)
Br8667_post_normal <- computeSumFactors(Br8667_post_normal, clusters = Br8667_post_normal$scran_quick_cluster)

Br6522_ant_normal <- logNormCounts(Br6522_ant_normal)
Br6522_mid_normal <- logNormCounts(Br6522_mid_normal)
Br8667_post_normal <- logNormCounts(Br8667_post_normal)

compare_gene_prop <- function(normal, all, title) {
    plot_df <- data.frame(
        "Mean_expr" = rowMeans(all@assays@data$logcounts),
        "Var_expr" = rowVars(all@assays@data$logcounts),
        "Mean_expr_normal" = rowMeans(normal@assays@data$logcounts),
        "Var_expr_normal" = rowVars(normal@assays@data$logcounts)
    )
    p1 <- ggplot(plot_df, aes(x = Mean_expr_normal, y = Mean_expr)) +
        geom_point(colour = "#00BFC4", alpha = 0.5) +
        theme_classic() +
        xlab("Mean (excl. artifacts)") +
        ylab("Mean (all spots)") +
        geom_abline(intercept = 0, slope = 1, lty = 2) +
        theme(axis.title = element_text(size = 10), plot.title = element_text(face = "bold", size = 10))
    p2 <- ggplot(plot_df, aes(x = Var_expr_normal, y = Var_expr)) +
        geom_point(colour = "#00BFC4", alpha = 0.5) +
        theme_classic() +
        xlab("Variance (excl. artifacts)") +
        ylab("Variance (all spots)") +
        geom_abline(intercept = 0, slope = 1, lty = 2) +
        theme(axis.title = element_text(size = 10)) +
        ggtitle("")
    title <- ggdraw() + draw_label(title, fontface = "bold")
    plot_grid(title, plot_grid(p1, p2, rel_widths = c(1, 1), nrow = 1, align = "h", labels = c("I", "II")), nrow = 2, rel_heights = c(0.1, 1))
}
p <- compare_gene_prop(Br6522_ant_normal, Br6522_ant, "Br6522 Anterior")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br6522_ant_gene_mean_gene_var.pdf", width = 4.5, height = 3)
p
dev.off()
p <- compare_gene_prop(Br6522_mid_normal, Br6522_mid, "Br6522 Middle")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br6522_mid_gene_mean_gene_var.pdf", width = 4.5, height = 3)
p
dev.off()
p <- compare_gene_prop(Br8667_post_normal, Br8667_post, "Br8667 Posterior")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br8667_post_gene_mean_gene_var.pdf", width = 4.5, height = 3)
p
dev.off()

Br6522_ant_dec <- modelGeneVar(Br6522_ant)
Br6522_mid_dec <- modelGeneVar(Br6522_mid)
Br8667_post_dec <- modelGeneVar(Br8667_post)

Br6522_ant_top.hvgs <- getTopHVGs(Br6522_ant_dec, prop = 0.1)
Br6522_mid_top.hvgs <- getTopHVGs(Br6522_mid_dec, prop = 0.1)
Br8667_post_top.hvgs <- getTopHVGs(Br8667_post_dec, prop = 0.1)

Br6522_ant <- runPCA(Br6522_ant, subset_row = Br6522_ant_top.hvgs, ncomponents = 50)
Br6522_mid <- runPCA(Br6522_mid, subset_row = Br6522_mid_top.hvgs, ncomponents = 50)
Br8667_post <- runPCA(Br8667_post, subset_row = Br8667_post_top.hvgs, ncomponents = 50)

Br6522_ant_normal_dec <- modelGeneVar(Br6522_ant_normal)
Br6522_mid_normal_dec <- modelGeneVar(Br6522_mid_normal)
Br8667_post_normal_dec <- modelGeneVar(Br8667_post_normal)

Br6522_ant_normal_top.hvgs <- getTopHVGs(Br6522_ant_normal_dec, prop = 0.1)
Br6522_mid_normal_top.hvgs <- getTopHVGs(Br6522_mid_normal_dec, prop = 0.1)
Br8667_post_normal_top.hvgs <- getTopHVGs(Br8667_post_normal_dec, prop = 0.1)

Br6522_ant_normal <- runPCA(Br6522_ant, subset_row = Br6522_ant_normal_top.hvgs, ncomponents = 50)
Br6522_mid_normal <- runPCA(Br6522_mid, subset_row = Br6522_mid_normal_top.hvgs, ncomponents = 50)
Br8667_post_normal <- runPCA(Br8667_post, subset_row = Br8667_post_normal_top.hvgs, ncomponents = 50)

# Find clusters
findClusters_sce <- function(sce) {
    output <- getClusteredPCs(reducedDim(sce))
    npcs <- metadata(output)$chosen
    cat(npcs)
    reducedDim(sce, "PCAsub") <- reducedDim(sce, "PCA")[, 1:npcs, drop = FALSE]
    g <- buildSNNGraph(sce, use.dimred = "PCAsub")
    cluster <- igraph::cluster_walktrap(g)$membership
    # Assigning to the 'colLabels' of the 'sce'.
    colLabels(sce) <- factor(cluster)
    table(colLabels(sce))
    sce
}

Br6522_ant <- findClusters_sce(Br6522_ant)
Br6522_mid <- findClusters_sce(Br6522_mid)
Br8667_post <- findClusters_sce(Br8667_post)

get_cluster_comp <- function(sobj, title) {
    Layers <- unique(sobj@colData$Layers)
    Layers <- Layers[order(Layers)]
    Labels <- levels(sobj@colData$label)
    layer_mat <- matrix(NA, nrow = length(Labels), ncol = length(Layers))
    for (i in c(1:dim(layer_mat)[1])) {
        for (j in c(1:dim(layer_mat)[2])) {
            layer_mat[i, j] <- sum(sobj@colData$Layers[sobj@colData$label == Labels[i]] == Layers[j])
        }
    }
    artifact_mat <- matrix(NA, nrow = length(Labels), ncol = 2)
    for (i in c(1:dim(layer_mat)[1])) {
        artifact_mat[i, 1] <- sum(!sobj@colData$is_wrinkle[sobj@colData$label == Labels[i]])
        artifact_mat[i, 2] <- sum(sobj@colData$is_wrinkle[sobj@colData$label == Labels[i]])
    }
    colnames(artifact_mat) <- c("None", "Artifact")
    colnames(layer_mat) <- Layers
    rownames(layer_mat) <- paste("Cluster", Labels)
    rownames(artifact_mat) <- paste("Cluster", Labels)
    colnames(layer_mat)[colnames(layer_mat) == "Unknown"] <- "NA"
    layer_mat <- reshape2::melt(layer_mat)
    artifact_mat <- reshape2::melt(artifact_mat)
    p1 <- ggplot(layer_mat, aes(x = Var1, y = value, fill = Var2)) +
        geom_bar(position = "dodge", stat = "identity") +
        theme_classic() +
        scale_fill_manual(name = "Layers", values = c(
            "Layer 1" = "#F0027F", "Layer 2" = "#377EB8", "Layer 3" = "#4DAF4A",
            "Layer 4" = "#984EA3", "Layer 5" = "#FFD700", "Layer 6" = "#FF7F00",
            "WM" = "#1A1A1A", "NA" = "transparent"
        )) +
        xlab("") +
        ylab("Number of spots") +
        ggtitle(title) +
        theme(plot.title = element_text(face = "bold"))

    p2 <- ggplot(artifact_mat, aes(x = Var1, y = value, fill = Var2)) +
        geom_bar(position = "fill", stat = "identity") +
        theme_classic() +
        scale_fill_manual(name = "", values = c("None" = "#619CFF", "Artifact" = "#F8766D")) +
        xlab("") +
        ylab("Fraction of spots") +
        theme(plot.title = element_text(face = "bold"))
    list(p1, p2)
}

p_list_1 <- get_cluster_comp(Br6522_ant, "Br6522 Anterior")
p_list_2 <- get_cluster_comp(Br6522_mid, "Br6522 Middle")
p_list_3 <- get_cluster_comp(Br8667_post, "Br8667 Posterior")

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_ant_cluster_layers.pdf", width = 8, height = 2.2)
p_list_1[[1]]
dev.off()

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_ant_cluster_artifact.pdf", width = 8, height = 2.2)
p_list_1[[2]]
dev.off()

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_mid_cluster_layers.pdf", width = 8, height = 2.2)
p_list_2[[1]]
dev.off()

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_mid_cluster_artifact.pdf", width = 8, height = 2.2)
p_list_2[[2]]
dev.off()

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br8667_post_cluster_layers.pdf", width = 8, height = 2.2)
p_list_3[[1]]
dev.off()

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br8667_post_cluster_artifact.pdf", width = 8, height = 2.2)
p_list_3[[2]]
dev.off()




Br6522_ant_normal <- findClusters_sce(Br6522_ant_normal)
Br6522_mid_normal <- findClusters_sce(Br6522_mid_normal)
Br8667_post_normal <- findClusters_sce(Br8667_post_normal)

dev.off()
percent.var <- attr(reducedDim(Br6522_ant, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)", pch = 20)
abline(v = chosen.elbow, col = "red")

percent.var <- attr(reducedDim(Br6522_mid, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)", pch = 20)
abline(v = chosen.elbow, col = "red")

percent.var <- attr(reducedDim(Br8667_post, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)", pch = 20)
abline(v = chosen.elbow, col = "red")

Br6522_ant <- runUMAP(Br6522_ant, dimred = "PCA")
colnames(reducedDim(Br6522_ant, "UMAP")) <- c("UMAP1", "UMAP2")
Br6522_mid <- runUMAP(Br6522_mid, dimred = "PCA")
colnames(reducedDim(Br6522_mid, "UMAP")) <- c("UMAP1", "UMAP2")
Br8667_post <- runUMAP(Br8667_post, dimred = "PCA")
colnames(reducedDim(Br8667_post, "UMAP")) <- c("UMAP1", "UMAP2")

Br6522_ant_normal <- runUMAP(Br6522_ant_normal, dimred = "PCA")
colnames(reducedDim(Br6522_ant_normal, "UMAP")) <- c("UMAP1", "UMAP2")
Br6522_mid_normal <- runUMAP(Br6522_mid_normal, dimred = "PCA")
colnames(reducedDim(Br6522_mid_normal, "UMAP")) <- c("UMAP1", "UMAP2")
Br8667_post_normal <- runUMAP(Br8667_post_normal, dimred = "PCA")
colnames(reducedDim(Br8667_post_normal, "UMAP")) <- c("UMAP1", "UMAP2")


viz_umap <- function(sobj, title) {
    df <- data.frame(reducedDim(sobj, "UMAP"))
    df$Layers <- sobj$Layers
    df$is_wrinkle <- sobj$is_wrinkle
    df$Artifact <- df$is_wrinkle
    df$Artifact[df$is_wrinkle] <- "Artifact"
    df$Artifact[!df$is_wrinkle] <- "Excluding Artifacts"
    df$Artifact <- factor(df$Artifact, levels = c("Excluding Artifacts", "Artifact"))
    df$Layers[df$Layers == "Unknown"] <- NA
    layer_palette <- c(
        "Layer 1" = "#F0027F", "Layer 2" = "#377EB8", "Layer 3" = "#4DAF4A",
        "Layer 4" = "#984EA3", "Layer 5" = "#FFD700", "Layer 6" = "#FF7F00",
        "WM" = "#1A1A1A", "NA" = "transparent"
    )
    ggplot(df, aes(x = UMAP1, y = UMAP2, colour = Layers)) +
        geom_point(size = 0.8, alpha = 0.7) +
        theme_bw() +
        scale_color_manual(name = "Layers", values = layer_palette) +
        facet_wrap(~Artifact) +
        ggtitle(title) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
}

p1 <- viz_umap(Br6522_ant, "Br6522 Anterior")
p2 <- viz_umap(Br6522_mid, "Br6522 Middle")
p3 <- viz_umap(Br8667_post, "Br8667 Posterior")

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_ant_UMAP.pdf", width = 6, height = 3)
p1
dev.off()

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_mid_UMAP.pdf", width = 6, height = 3)
p2
dev.off()

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br8667_post_UMAP.pdf", width = 6, height = 3)
p3
dev.off()

get_PC_var_explained <- function(sobj) {
    df <- data.frame(reducedDim(sobj, "PCA"))
    df <- df[, 1:10]
    df$Layers <- sobj$Layers
    df$is_wrinkle <- sobj$is_wrinkle
    df$libSize <- sobj$sum_gene
    df$percent_mito <- sobj$expr_chrM_ratio
    res <- anova(lm(PC1 ~ Layers + is_wrinkle + libSize + percent_mito, df))
    all_res <- res$`Sum Sq` / sum(res$`Sum Sq`)
    res <- anova(lm(PC2 ~ Layers + is_wrinkle + libSize + percent_mito, df))
    all_res <- rbind(all_res, res$`Sum Sq` / sum(res$`Sum Sq`))
    res <- anova(lm(PC3 ~ Layers + is_wrinkle + libSize + percent_mito, df))
    all_res <- rbind(all_res, res$`Sum Sq` / sum(res$`Sum Sq`))
    res <- anova(lm(PC4 ~ Layers + is_wrinkle + libSize + percent_mito, df))
    all_res <- rbind(all_res, res$`Sum Sq` / sum(res$`Sum Sq`))
    res <- anova(lm(PC5 ~ Layers + is_wrinkle + libSize + percent_mito, df))
    all_res <- rbind(all_res, res$`Sum Sq` / sum(res$`Sum Sq`))
    res <- anova(lm(PC6 ~ Layers + is_wrinkle + libSize + percent_mito, df))
    all_res <- rbind(all_res, res$`Sum Sq` / sum(res$`Sum Sq`))
    res <- anova(lm(PC7 ~ Layers + is_wrinkle + libSize + percent_mito, df))
    all_res <- rbind(all_res, res$`Sum Sq` / sum(res$`Sum Sq`))
    res <- anova(lm(PC8 ~ Layers + is_wrinkle + libSize + percent_mito, df))
    all_res <- rbind(all_res, res$`Sum Sq` / sum(res$`Sum Sq`))
    res <- anova(lm(PC9 ~ Layers + is_wrinkle + libSize + percent_mito, df))
    all_res <- rbind(all_res, res$`Sum Sq` / sum(res$`Sum Sq`))
    res <- anova(lm(PC10 ~ Layers + is_wrinkle + libSize + percent_mito, df))
    all_res <- rbind(all_res, res$`Sum Sq` / sum(res$`Sum Sq`))
    all_res <- all_res[, c(1, 2, 3, 4)]
    rownames(all_res) <- c(1:10)
    colnames(all_res) <- c("Layers", "Is Wrinkle", "Library Size", "Percent mito")
    all_res <- data.frame(all_res)
    all_res$PC <- paste0("PC", c(1:10))
    all_res <- reshape2::melt(all_res, id.vars = "PC")
    all_res$variable <- as.character(all_res$variable)
    all_res$variable[all_res$variable == "Is.Wrinkle"] <- "Is wrinkle"
    all_res$variable[all_res$variable == "Library.Size"] <- "Library size"
    all_res$variable[all_res$variable == "Percent.mito"] <- "Percent mito"
    all_res$variable <- as.factor(all_res$variable)
    all_res
}

Br6522_ant_anova <- get_PC_var_explained(Br6522_ant)
Br6522_mid_anova <- get_PC_var_explained(Br6522_mid)
Br8667_post_anova <- get_PC_var_explained(Br8667_post)

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_ant_PC_anova.pdf", width = 5, height = 3.5)
ggplot(aes(x = variable, y = PC, fill = value), data = Br6522_ant_anova) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#416A1D") +
    labs(y = NULL, x = NULL, fill = "R-squared") +
    geom_text(aes(label = round(value, 2)), size = 3) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks.length = unit(0, "lines")) +
    scale_y_discrete(limits = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")) +
    ggtitle("Br6522 Anterior") +
    theme(plot.title = element_text(face = "bold"))
dev.off()

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_mid_PC_anova.pdf", width = 5, height = 3.5)
ggplot(aes(x = variable, y = PC, fill = value), data = Br6522_mid_anova) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#416A1D") +
    labs(y = NULL, x = NULL, fill = "R-squared") +
    geom_text(aes(label = round(value, 2)), size = 3) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks.length = unit(0, "lines")) +
    scale_y_discrete(limits = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")) +
    ggtitle("Br6522 Middle") +
    theme(plot.title = element_text(face = "bold"))
dev.off()

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br8667_post_PC_anova.pdf", width = 5, height = 3.5)
ggplot(aes(x = variable, y = PC, fill = value), data = Br8667_post_anova) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#416A1D") +
    labs(y = NULL, x = NULL, fill = "R-squared") +
    geom_text(aes(label = round(value, 2)), size = 3) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks.length = unit(0, "lines")) +
    scale_y_discrete(limits = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")) +
    ggtitle("Br8667 Posterior") +
    theme(plot.title = element_text(face = "bold"))
dev.off()

Br6522_ant$Layer_wrinkle <- paste(Br6522_ant$Layers, as.numeric(Br6522_ant$is_wrinkle), sep = "_")
Br6522_mid$Layer_wrinkle <- paste(Br6522_mid$Layers, as.numeric(Br6522_mid$is_wrinkle), sep = "_")
Br8667_post$Layer_wrinkle <- paste(Br8667_post$Layers, as.numeric(Br8667_post$is_wrinkle), sep = "_")

plot_DE_genes <- function(df_list, opt = "all", quantile_thresh = 0.95) {
    df_list <- lapply(df_list, function(df) {
        if (is.null(df)) {
            NULL
        } else {
            df <- data.frame(df)
            df$gene <- rownames(df)
            rownames(df) <- c(1:dim(df)[1])
            df
        }
    })
    if (opt == "all") {
        layer_1 <- df_list[[1]]
        layer_1$logFC.Layer.1_1 <- NULL
        layer_1$Layer <- "Layer 1"
        layer_1 <- layer_1[layer_1$FDR < 0.01, ]
        layer_1 <- layer_1[abs(layer_1$summary.logFC) >= quantile(abs(layer_1$summary.logFC), quantile_thresh), ]
        # layer_1 <- layer_1[abs(layer_1$summary.logFC) >= 0.2, ]
        layer_1 <- layer_1[order(abs(layer_1$summary.logFC)), ]

        layer_WM <- df_list[[7]]
        layer_WM$Layer <- "WM"
        layer_WM$logFC.WM_1 <- NULL
        layer_WM <- layer_WM[layer_WM$FDR < 0.01, ]
        layer_WM <- layer_WM[abs(layer_WM$summary.logFC) >= quantile(abs(layer_WM$summary.logFC), quantile_thresh), ]
        # layer_WM <- layer_WM[abs(layer_WM$summary.logFC) >= 0.2, ]
        layer_WM <- layer_WM[order(abs(layer_WM$summary.logFC)), ]
    }

    layer_2 <- df_list[[2]]
    layer_2$logFC.Layer.2_1 <- NULL
    layer_2$Layer <- "Layer 2"
    layer_2 <- layer_2[layer_2$FDR < 0.01, ]
    layer_2 <- layer_2[abs(layer_2$summary.logFC) >= quantile(abs(layer_2$summary.logFC), quantile_thresh), ]
    # layer_2 <- layer_2[abs(layer_2$summary.logFC) >= 0.2, ]
    layer_2 <- layer_2[order(abs(layer_2$summary.logFC)), ]

    layer_3 <- df_list[[3]]
    layer_3$logFC.Layer.3_1 <- NULL
    layer_3$Layer <- "Layer 3"
    layer_3 <- layer_3[layer_3$FDR < 0.01, ]
    layer_3 <- layer_3[abs(layer_3$summary.logFC) >= quantile(abs(layer_3$summary.logFC), quantile_thresh), ]
    # layer_3 <- layer_3[abs(layer_3$summary.logFC) >= 0.2, ]
    layer_3 <- layer_3[order(abs(layer_3$summary.logFC)), ]

    layer_4 <- df_list[[4]]
    layer_4$logFC.Layer.4_1 <- NULL
    layer_4$Layer <- "Layer 4"
    layer_4 <- layer_4[layer_4$FDR < 0.01, ]
    # layer_4 <- layer_4[abs(layer_4$summary.logFC) >= 0.2, ]
    layer_4 <- layer_4[abs(layer_4$summary.logFC) >= quantile(abs(layer_4$summary.logFC), quantile_thresh), ]
    layer_4 <- layer_4[order(abs(layer_4$summary.logFC)), ]

    layer_5 <- df_list[[5]]
    layer_5$logFC.Layer.5_1 <- NULL
    layer_5$Layer <- "Layer 5"
    layer_5 <- layer_5[layer_5$FDR < 0.01, ]
    layer_5 <- layer_5[abs(layer_5$summary.logFC) >= quantile(abs(layer_5$summary.logFC), quantile_thresh), ]
    # layer_5 <- layer_5[abs(layer_5$summary.logFC) >= 0.2, ]
    layer_5 <- layer_5[order(abs(layer_5$summary.logFC)), ]

    layer_6 <- df_list[[6]]
    layer_6$logFC.Layer.6_1 <- NULL
    layer_6$Layer <- "Layer 6"
    layer_6 <- layer_6[layer_6$FDR < 0.01, ]
    # layer_6 <- layer_6[abs(layer_6$summary.logFC) >= 0.2, ]
    layer_6 <- layer_6[abs(layer_6$summary.logFC) >= quantile(abs(layer_6$summary.logFC), quantile_thresh), ]
    layer_6 <- layer_6[order(abs(layer_6$summary.logFC)), ]

    if (opt == "all") {
        plot_df <- rbind(layer_1, layer_2, layer_3, layer_4, layer_5, layer_6, layer_WM)
    } else {
        plot_df <- rbind(layer_2, layer_3, layer_4, layer_5, layer_6)
    }
    plot_df$Absolute_logFC <- abs(plot_df$summary.logFC)
    plot_df$genes_ordered <- factor(plot_df$gene, levels = unique(plot_df$gene))
    ggplot(plot_df, aes(x = Layer, y = genes_ordered, size = Absolute_logFC, color = FDR)) +
        geom_point(alpha = 0.8) +
        theme_bw() +
        xlab("") +
        ylab("") +
        scale_colour_gradient(low = "#3D07F4", high = "#E2DBF8") +
        scale_size_continuous(name = "|log FC|") +
        theme(panel.grid.major = element_line(colour = "black", linewidth = 0.2, linetype = "dashed")) +
        theme(
            panel.grid.major.x = element_blank()
        )
}

plot_DE_violin_plot <- function(df_list, opt = "all") {
    df_list <- lapply(df_list, function(df) {
        data.frame(df)
    })
    if (opt == "all") {
        layer_1 <- df_list[[1]]
        layer_1$logFC.Layer.1_1 <- NULL
        layer_1$Layer <- "Layer 1"
        # layer_1 <- layer_1[layer_1$FDR < 0.01, ]
    }
    layer_2 <- df_list[[2]]
    layer_2$logFC.Layer.2_1 <- NULL
    layer_2$Layer <- "Layer 2"
    # layer_2 <- layer_2[layer_2$FDR < 0.01, ]

    layer_3 <- df_list[[3]]
    layer_3$Layer <- "Layer 3"
    layer_3$logFC.Layer.3_1 <- NULL
    # layer_3 <- layer_3[layer_3$FDR < 0.01, ]

    layer_4 <- df_list[[4]]
    layer_4$Layer <- "Layer 4"
    layer_4$logFC.Layer.4_1 <- NULL
    # layer_4 <- layer_4[layer_4$FDR < 0.01, ]

    layer_5 <- df_list[[5]]
    layer_5$Layer <- "Layer 5"
    layer_5$logFC.Layer.5_1 <- NULL
    # layer_5 <- layer_5[layer_5$FDR < 0.01, ]

    layer_6 <- df_list[[6]]
    layer_6$Layer <- "Layer 6"
    layer_6$logFC.Layer.6_1 <- NULL
    # layer_6 <- layer_6[layer_6$FDR < 0.01, ]

    if (opt == "all") {
        layer_WM <- df_list[[7]]
        layer_WM$Layer <- "WM"
        layer_WM$logFC.WM_1 <- NULL
        # layer_WM <- layer_WM[layer_WM$FDR < 0.01, ]
    }
    if (opt == "all") {
        plot_df <- rbind(layer_1, layer_2, layer_3, layer_4, layer_5, layer_6, layer_WM)
    } else {
        plot_df <- rbind(layer_2, layer_3, layer_4, layer_5, layer_6)
    }
    plot_df$pvalue_cat <- plot_df$FDR
    plot_df$pvalue_cat[plot_df$FDR > 0.1] <- "FDR > 0.1"
    plot_df$pvalue_cat[plot_df$FDR < 0.1 & plot_df$FDR > 0.05] <- "0.05 < FDR < 0.1"
    plot_df$pvalue_cat[plot_df$FDR < 0.05 & plot_df$FDR > 0.01] <- "0.01 < FDR < 0.05"
    plot_df$pvalue_cat[plot_df$FDR < 0.01 & plot_df$FDR > 0.001] <- "0.001 < FDR < 0.01"
    plot_df$pvalue_cat[plot_df$FDR < 0.001] <- "FDR < 0.001"
    plot_df$pvalue_cat <- factor(plot_df$pvalue_cat, levels = c(
        "FDR > 0.1", "0.05 < FDR < 0.1", "0.01 < FDR < 0.05", "0.001 < FDR < 0.01",
        "FDR < 0.001"
    ))
    plot_df.summary <- plot_df %>%
        group_by(Layer, pvalue_cat) %>%
        summarise(
            mad = mad(summary.logFC),
            median = median(summary.logFC)
        )
    plot_df.summary$lowerLim <- plot_df.summary$median - 3 * plot_df.summary$mad
    plot_df.summary$upperLim <- plot_df.summary$median + 3 * plot_df.summary$mad
    ggplot(plot_df.summary, aes(x = Layer, y = median)) +
        xlab("") +
        ylab("Log Fold-change") +
        geom_hline(yintercept = -0.1, lty = 2, linewidth = 0.3) +
        geom_hline(yintercept = 0.1, lty = 2, linewidth = 0.3) +
        theme_classic() +
        geom_point(aes(color = pvalue_cat), position = position_dodge(0.5)) +
        geom_errorbar(aes(ymin = lowerLim, ymax = upperLim, color = pvalue_cat), position = position_dodge(0.5), width = 0) +
        ylim(c(-0.5, 0.5)) +
        scale_color_manual(name = "", values = c(
            "FDR > 0.1" = "#C3B2F9",
            "0.05 < FDR < 0.1" = "#BBA8F8",
            "0.01 < FDR < 0.05" = "#9B7DF8",
            "0.001 < FDR < 0.01" = "#784FF8",
            "FDR < 0.001" = "#4A12FA"
        ))
}

Br6522_ant_out <- pairwiseTTests(Br6522_ant, groups = Br6522_ant$Layer_wrinkle)
Br6522_ant_out_all <- scran::combineMarkers(
    de.lists = Br6522_ant_out$statistics[paste(Br6522_ant_out$pairs$first, Br6522_ant_out$pairs$second) %in% c(
        "Layer 1_0 Layer 1_1",
        "Layer 2_0 Layer 2_1", "Layer 3_0 Layer 3_1",
        "Layer 4_0 Layer 4_1", "Layer 5_0 Layer 5_1",
        "Layer 6_0 Layer 6_1", "WM_0 WM_1"
    )],
    pairs = Br6522_ant_out$pairs[paste(Br6522_ant_out$pairs$first, Br6522_ant_out$pairs$second) %in% c(
        "Layer 1_0 Layer 1_1",
        "Layer 2_0 Layer 2_1", "Layer 3_0 Layer 3_1",
        "Layer 4_0 Layer 4_1", "Layer 5_0 Layer 5_1",
        "Layer 6_0 Layer 6_1", "WM_0 WM_1"
    ), ],
    pval.type = "all"
)
Br6522_ant_DE_plots <- plot_DE_genes(Br6522_ant_out_all, "all", 0.98) + ggtitle("Br6522 Anterior") + theme(plot.title = element_text(face = "bold"))
Br6522_ant_vln_plots <- plot_DE_violin_plot(Br6522_ant_out_all, "all") + ggtitle("Br6522 Anterior") + theme(plot.title = element_text(face = "bold"))

pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_ant_topDEgenes.pdf", width = 5, height = 5)
Br6522_ant_DE_plots
dev.off()
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_ant_logFC_dotPlot.pdf", width = 6, height = 3.5)
Br6522_ant_vln_plots + theme(legend.position = "none")
dev.off()

Br6522_mid_out <- pairwiseTTests(Br6522_mid, groups = Br6522_mid$Layer_wrinkle)
Br6522_mid_out_all <- scran::combineMarkers(
    de.lists = Br6522_mid_out$statistics[paste(Br6522_mid_out$pairs$first, Br6522_mid_out$pairs$second) %in% c(
        "Layer 1_0 Layer 1_1",
        "Layer 2_0 Layer 2_1", "Layer 3_0 Layer 3_1",
        "Layer 4_0 Layer 4_1", "Layer 5_0 Layer 5_1",
        "Layer 6_0 Layer 6_1", "WM_0 WM_1"
    )],
    pairs = Br6522_mid_out$pairs[paste(Br6522_mid_out$pairs$first, Br6522_mid_out$pairs$second) %in% c(
        "Layer 1_0 Layer 1_1",
        "Layer 2_0 Layer 2_1", "Layer 3_0 Layer 3_1",
        "Layer 4_0 Layer 4_1", "Layer 5_0 Layer 5_1",
        "Layer 6_0 Layer 6_1", "WM_0 WM_1"
    ), ],
    pval.type = "all"
)
Br6522_mid_DE_plots <- plot_DE_genes(Br6522_mid_out_all, "all", 0.05) + ggtitle("Br6522 Middle") + theme(plot.title = element_text(face = "bold"))
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_mid_topDEgenes.pdf", width = 5, height = 5)
Br6522_mid_DE_plots
dev.off()
Br6522_mid_vln_plots <- plot_DE_violin_plot(Br6522_mid_out_all, "all") + ggtitle("Br6522 Middle") + theme(plot.title = element_text(face = "bold"))
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br6522_mid_logFC_dotPlot.pdf", width = 6, height = 3.5)
Br6522_mid_vln_plots + theme(legend.position = "none")
dev.off()

Br8667_post_out <- pairwiseTTests(Br8667_post, groups = Br8667_post$Layer_wrinkle)
Br8667_post_out_some <- scran::combineMarkers(
    de.lists = Br8667_post_out$statistics[paste(Br8667_post_out$pairs$first, Br8667_post_out$pairs$second) %in% c(
        "Layer 1_0 Layer 1_1",
        "Layer 2_0 Layer 2_1", "Layer 3_0 Layer 3_1",
        "Layer 4_0 Layer 4_1", "Layer 5_0 Layer 5_1",
        "Layer 6_0 Layer 6_1", "WM_0 WM_1"
    )],
    pairs = Br8667_post_out$pairs[paste(Br8667_post_out$pairs$first, Br8667_post_out$pairs$second) %in% c(
        "Layer 1_0 Layer 1_1",
        "Layer 2_0 Layer 2_1", "Layer 3_0 Layer 3_1",
        "Layer 4_0 Layer 4_1", "Layer 5_0 Layer 5_1",
        "Layer 6_0 Layer 6_1", "WM_0 WM_1"
    ), ],
    pval.type = "all"
)
Br8667_post_out_all <- list()
Br8667_post_out_all[[1]] <- NULL
Br8667_post_out_all[[2]] <- Br8667_post_out_some[[1]]
Br8667_post_out_all[[3]] <- Br8667_post_out_some[[2]]
Br8667_post_out_all[[4]] <- Br8667_post_out_some[[3]]
Br8667_post_out_all[[5]] <- Br8667_post_out_some[[4]]
Br8667_post_out_all[[6]] <- Br8667_post_out_some[[5]]
Br8667_post_out_all[[7]] <- NULL
Br8667_post_DE_plots <- plot_DE_genes(Br8667_post_out_all, "some", 0.8) + ggtitle("Br8667 Posterior") + theme(plot.title = element_text(face = "bold"))
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br8667_post_topDEgenes.pdf", width = 5, height = 5)
Br8667_post_DE_plots
dev.off()
Br8667_post_vln_plots <- plot_DE_violin_plot(Br8667_post_out_all, "some") + ggtitle("Br8667 Posterior") + theme(plot.title = element_text(face = "bold"))
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/Br8667_post_logFC_dotPlot.pdf", width = 6, height = 3.5)
Br8667_post_vln_plots + theme(legend.position = "none")
dev.off()

get_legend <- function(myggplot) {
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

Br8667_post_vln_plots <- Br8667_post_vln_plots + theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 2, byrow = TRUE))
FDR_legend_forDotPlot <- get_legend(Br8667_post_vln_plots)
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S5/legend_logFC_dotPlot.pdf", width = 6, height = 1)
plot_grid(FDR_legend_forDotPlot)
dev.off()
