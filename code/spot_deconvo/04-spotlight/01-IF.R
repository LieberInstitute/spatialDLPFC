library('ggplot2')
library('SPOTlight')
library('SingleCellExperiment')
library('SpatialExperiment')
library('scater')
library('scran')
library('here')

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

cell_type_var = 'cellType_broad_hc'

#   Used for downsampling single-cell object prior to training
n_cells_per_type <- 100

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

#   Load objects
load(sce_in, verbose = TRUE)
spe_IF <- readRDS(spe_in)
gc()

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
genes <- !grepl(pattern = "^(Rp[ls]|MT-)", x = rownames(sce)) # "^Rp[l|s]|Mt"

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
