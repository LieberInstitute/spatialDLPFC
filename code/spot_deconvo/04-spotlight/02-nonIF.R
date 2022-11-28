library("ggplot2")
library("ggcorrplot")
library("scatterpie")
library("SPOTlight")
library("SingleCellExperiment")
library("SpatialExperiment")
library("scater")
library("scran")
library("here")
library("NMF")
library("sessioninfo")
library("tidyverse")

################################################################################
#   Variable definitions
################################################################################

cell_group <- "layer" # "broad" or "layer"

sce_in <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("sce_", cell_group, ".rds")
)
spe_in <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
)

marker_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("markers_", cell_group, ".txt")
)
marker_stats_path <- here(
    "processed-data", "spot_deconvo", "05-shared_utilities",
    paste0("marker_stats_", cell_group, ".rds")
)

plot_dir <- here(
    "plots", "spot_deconvo", "04-spotlight", "nonIF", cell_group
)
processed_dir <- here(
    "processed-data", "spot_deconvo", "04-spotlight", "nonIF", cell_group
)

#   Column names in colData(sce)
if (cell_group == "broad") {
    cell_type_var <- "cellType_broad_hc"
} else {
    cell_type_var <- "layer_level"
}

#   Column names in rowData(sce)
symbol_var <- "gene_name"

#   Column names in colData(spe)
sample_var <- "sample_id"
cell_count_var <- "count"

#   Used for downsampling single-cell object prior to training
n_cells_per_type <- 100

################################################################################
#   Preprocessing
################################################################################

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

#   Load objects
sce <- readRDS(sce_in)
load(spe_in, verbose = TRUE)
gc()

dec <- modelGeneVar(sce)

#   Plot gene-expression variance
pdf(file.path(plot_dir, "expr_variance.pdf"))
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
dev.off()

#   Get the top 3000 highly variable genes.
hvg <- getTopHVGs(dec, n = 3000)

colLabels(sce) <- colData(sce)[[cell_type_var]]

#   Read in markers found by Louise's method
markers <- readLines(marker_path)
marker_stats <- readRDS(marker_stats_path)
stopifnot(all(markers %in% rownames(sce)))

#   Make sure there is only one row per gene, and that row corresponds to the
#   target cell type for which the gene could potentially be a marker
marker_stats <- marker_stats %>%
    group_by(gene) %>%
    filter(ratio == max(ratio))


#   Filter out any mitochondrial/ribosomal genes
genes <- !grepl(
    pattern = "^(RP[LS]|MT-)",
    x = rowData(sce)[[symbol_var]]
) & (rownames(sce) %in% markers)
print(
    paste(
        length(which(genes)),
        "markers used after filtering out mitochondrial/ribosomal genes."
    )
)

#   Score markers (determine weights) and form a data frame SPOTlight can use as
#   input
mgs <- scoreMarkers(sce, subset.row = genes)

mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i

    #   Only take genes that have already been determined to be markers for this
    #   cell type
    x <- x[
        marker_stats$cellType.target[
            match(rownames(x), marker_stats$gene)
        ] == i,
    ]

    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

#   It's technically possible for the above 'lapply' to generate data frames
#   with non-unique rownames, leading to R automatically changing some rownames
#   and thus illegitimate gene names. Verify this most likely didn't happen
stopifnot(all(rownames(mgs_df) %in% rownames(sce)))

#-------------------------------------------------------------------------------

# split cell indices by identity
idx <- split(seq(ncol(sce)), sce[[cell_type_var]])

#   This was slightly changed from the tutorial for simplicity
cs_keep <- lapply(idx, function(i) sample(i, min(length(i), n_cells_per_type)))
sce <- sce[, unlist(cs_keep)]

################################################################################
#   Train model and deconvolve cell types
################################################################################

res <- SPOTlight(
    x = sce,
    y = spe,
    groups = as.character(sce[[cell_type_var]]),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene"
)

saveRDS(res, file.path(processed_dir, "results.rds"))

################################################################################
#   Visualization
################################################################################

#   Extract deconvolution matrix and NMF model fit
mat <- res$mat
mod <- res$NMF

pdf(file.path(plot_dir, "topic_profiles.pdf"))
plotTopicProfiles(
    x = mod,
    y = sce[[cell_type_var]],
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1
) +
    theme(aspect.ratio = 1)

plotTopicProfiles(
    x = mod,
    y = sce[[cell_type_var]],
    facet = TRUE,
    min_prop = 0.01,
    ncol = 6
)
dev.off()

sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))

pdf(file.path(plot_dir, "visualizations.pdf"))
plotCorrelationMatrix(mat)
plotInteractions(mat, which = "heatmap", metric = "prop")
plotInteractions(mat, which = "heatmap", metric = "jaccard")
plotInteractions(mat, which = "network")
dev.off()

ct <- colnames(mat)
mat[mat < 0.1] <- 0

# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
    "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db",
    "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
    "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d"
)

pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

spe$res_ss <- res[[2]][colnames(spe)]
xy <- spatialCoords(spe)
spe$x <- xy[, 1]
spe$y <- xy[, 2]

#   Split spatial-related plots by sample
for (sample_id in unique(spe[[sample_var]])) {
    #   Subset objects to this sample
    this_sample_indices <- which(spe[[sample_var]] == sample_id)
    temp_spe <- spe[, this_sample_indices]
    temp_mat <- mat[this_sample_indices, ]

    #   Scatterpie
    pdf(file.path(plot_dir, paste0("scatterpie_", sample_id, ".pdf")))
    print(
        plotSpatialScatterpie(
            x = temp_spe,
            y = temp_mat,
            cell_types = colnames(temp_mat),
            img = FALSE,
            scatterpie_alpha = 1,
            pie_scale = 0.4
        ) +
            scale_fill_manual(
                values = pal,
                breaks = names(pal)
            )
    )
    dev.off()

    #   Residuals
    pdf(file.path(plot_dir, paste0("residuals_", sample_id, ".pdf")))
    print(
        ggcells(temp_spe, aes(x, y, color = res_ss)) +
            geom_point() +
            scale_color_viridis_c() +
            coord_fixed() +
            theme_bw()
    )
    dev.off()
}

#   Save final objects
saveRDS(sce, file.path(processed_dir, "sce.rds"))
saveRDS(spe, file.path(processed_dir, "spe.rds"))

################################################################################
#   Export 'clusters.csv' file of cell counts
################################################################################

#   '/' chracters in cell types are problematic and should be replaced with '_'
cell_types <- gsub("/", "_", colnames(res$mat))

#   Create a data frame with cell counts for all samples, and add the 'key'
#   column. Note here we scale cell-type proportions by total cells per spot,
#   the latter of which is computed prior to running SPOTlight
clusters <- data.frame(res$mat * spe[[cell_count_var]])
colnames(clusters) <- cell_types

clusters$key <- spe$key
clusters <- clusters[, c("key", cell_types)]

#   Write individual 'clusters.csv' files for each sample
for (sample_id in unique(spe[[sample_var]])) {
    clusters_small <- clusters[spe[[sample_var]] == sample_id, ]

    #   Make sure processed directory exists for this sample
    dir.create(
        file.path(processed_dir, sample_id),
        recursive = TRUE,
        showWarnings = FALSE
    )

    write.csv(
        clusters_small,
        file.path(processed_dir, sample_id, "clusters.csv"),
        row.names = FALSE, quote = FALSE
    )
}

session_info()
