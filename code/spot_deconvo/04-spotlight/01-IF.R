library('ggplot2')
library('ggcorrplot')
library('scatterpie')
library('SPOTlight')
library('SingleCellExperiment')
library('SpatialExperiment')
library('scater')
library('scran')
library('here')
library('NMF')
library('sessioninfo')

################################################################################
#   Variable definitions
################################################################################

sce_in <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"
spe_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)

marker_path <- here(
    "processed-data", "spot_deconvo", "01-tangram", "markers.txt"
)

plot_dir = here(
    "plots", "spot_deconvo", "04-spotlight", "IF"
)
processed_dir = here(
    "processed-data", "spot_deconvo", "04-spotlight", "IF"
)

#   Column names in colData(sce)
cell_type_var = 'cellType_broad_hc'

#   Column names in rowData(sce)
ensembl_id_var = 'gene_id'
symbol_var = 'gene_name'

#   Used for downsampling single-cell object prior to training
n_cells_per_type <- 100

################################################################################
#   Preprocessing
################################################################################

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

#   Load objects
load(sce_in, verbose = TRUE)
spe <- readRDS(spe_in)
gc()

rownames(sce) = rowData(sce)[[ensembl_id_var]]

#   Drop rare cell types for single-cell data
print("Distribution of cells to keep (FALSE) vs. drop (TRUE):")
table(sce[[cell_type_var]] == "EndoMural")
sce <- sce[, sce[[cell_type_var]] != "EndoMural"]

dec <- modelGeneVar(sce)

#   Plot gene-expression variance
pdf(file.path(plot_dir, 'expr_variance.pdf'))
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
dev.off()

# Get the top 3000 genes.
hvg <- getTopHVGs(dec, n = 3000)

#-------------------------------------------------------------------------------
#   TODO: import Louise's markers instead. Well, can we assign each of our
#   marker genes a score in a way compatible with SPOTlight?
#-------------------------------------------------------------------------------

colLabels(sce) <- colData(sce)[[cell_type_var]]

# Get vector indicating which genes are neither ribosomal or mitochondrial
#   The regular expression in the tutorial appears to be bad, and was fixed here
genes <- !grepl(pattern = "^(Rp[ls]|MT-)", x = rowData(sce)[[symbol_var]]) # "^Rp[l|s]|Mt"

# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)

mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

#-------------------------------------------------------------------------------

# split cell indices by identity
idx <- split(seq(ncol(sce)), sce[[cell_type_var]])

#   This was slightly changed from the tutorial for simplicity
cs_keep <- lapply(idx, function(i) sample(i, min(length(i), n_cells_per_type)))
sce <- sce[, unlist(cs_keep)]

################################################################################
#   Train model and deconvolve cell types
################################################################################

#   Train model
mod_ls <- trainNMF(
    x = sce,
    y = spe,
    groups = as.character(sce[[cell_type_var]]),
    mgs = mgs_df,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene"
)

# Run deconvolution
res <- runDeconvolution(
    x = spe,
    mod = mod_ls[["mod"]],
    ref = mod_ls[["topic"]]
)

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

################################################################################
#   Visualization
################################################################################

#   Extract deconvolution matrix and NMF model fit
mat <- res$mat
mod <- res$NMF

pdf(file.path(plot_dir, 'topic_profiles.pdf'))
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

pdf(file.path(plot_dir, 'visualizations.pdf'))
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
for (sample_id in unique(spe$sample_id)) {
    this_sample_indices = which(spe$sample_id == sample_id)
    temp_spe = spe[, this_sample_indices]
    temp_mat = mat[this_sample_indices,]
    
    #   Scatterpie
    pdf(file.path(plot_dir, paste0('scatterpie_', sample_id, '.pdf')))
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
    pdf(file.path(plot_dir, paste0('residuals_', sample_id, '.pdf')))
    print(
        ggcells(temp_spe, aes(x, y, color = res_ss)) +
            geom_point() +
            scale_color_viridis_c() +
            coord_fixed() +
            theme_bw()
    )
    dev.off()
}

session_info()
