# https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/asd_snRNAseq_recast.R
# https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/ad_snRNAseq_recast.R
# http://research.libd.org/spatialLIBD/reference/layer_stat_cor.html
# http://research.libd.org/spatialLIBD/reference/layer_stat_cor_plot.html

library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)
library(edgeR)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# load pseudo bulked spe object. made in /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/analysis/09_region_differential_expression/09_region_differential_expression.R
load(file = here::here("processed-data", "rdata", "spe", "pseudo_bulked_spe", paste0("sce_pseudobulk_bayesSpace_normalized_filtered_k", k, ".Rdata")), verbose = TRUE)
sce_pseudobulk_bayesSpace <- spe_pseudo


###############################
##### get mean expression  ####
mat <- assays(sce_pseudobulk_bayesSpace)$logcounts # make matrix of just the log normalized counts
save(mat, file = here::here("processed-data", "rdata", "spe", "07_spatial_registration", paste0("dlpfc_pseudobulked_mat_k", k, ".Rdata")))


#####################
## Build a group model

# convert variables to factors
colData(sce_pseudobulk_bayesSpace)$spatial.cluster <- as.factor(colData(sce_pseudobulk_bayesSpace)[[paste0("bayesSpace_harmony_", k)]])
colData(sce_pseudobulk_bayesSpace)$region <- as.factor(colData(sce_pseudobulk_bayesSpace)$region)
colData(sce_pseudobulk_bayesSpace)$age <- as.integer(colData(sce_pseudobulk_bayesSpace)$age)
colData(sce_pseudobulk_bayesSpace)$sex <- as.factor(colData(sce_pseudobulk_bayesSpace)$sex)
colData(sce_pseudobulk_bayesSpace)$diagnosis <- as.factor(colData(sce_pseudobulk_bayesSpace)$diagnosis)
colData(sce_pseudobulk_bayesSpace)$subject <- as.factor(colData(sce_pseudobulk_bayesSpace)$subject)

# create matrix where the rownames are the sample:clusters and the columns are the other variales (spatial.cluster + region + age + sex)
mod <- with(
    colData(sce_pseudobulk_bayesSpace),
    model.matrix(~ 0 + spatial.cluster + region + age + sex)
) # removed diagnosis cuz it only has 1 level
colnames(mod) <- gsub("cluster", "", colnames(mod))
# why is regtionanterior missing in the colnames(mod)?

## get duplicate correlation #http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/dupcor.html
corfit <- duplicateCorrelation(mat, mod,
    block = sce_pseudobulk_bayesSpace$sample_id
)
save(corfit, file = here::here("processed-data", "rdata", "spe", "07_spatial_registration", paste0("dlpfc_pseudobulked_bayesSpace_dupCor_k", k, ".Rdata")))

## Next for each layer test that layer vs the rest
cell_idx <- splitit(sce_pseudobulk_bayesSpace$spatial.cluster)

eb0_list_cell <- lapply(cell_idx, function(x) {
    res <- rep(0, ncol(sce_pseudobulk_bayesSpace))
    res[x] <- 1
    m <- with(
        colData(sce_pseudobulk_bayesSpace),
        model.matrix(~ res +
            region + age + sex)
    )
    eBayes(
        lmFit(
            mat,
            design = m,
            block = sce_pseudobulk_bayesSpace$sample_id,
            correlation = corfit$consensus.correlation
        )
    )
})
save(eb0_list_cell, file = here::here("processed-data", "rdata", "spe", "07_spatial_registration", paste0("dlpfc_pseudobulked_bayesSpace_specific_Ts_k", k, ".Rdata")))


##########
## Extract the p-values

pvals0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cell) <- rownames(mat)

t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) <- rownames(mat)
fdrs0_contrasts_cell <- apply(pvals0_contrasts_cell, 2, p.adjust, "fdr")

data.frame(
    "FDRsig" = colSums(fdrs0_contrasts_cell < 0.05 &
        t0_contrasts_cell > 0),
    "Pval10-6sig" = colSums(pvals0_contrasts_cell < 1e-6 &
        t0_contrasts_cell > 0),
    "Pval10-8sig" = colSums(pvals0_contrasts_cell < 1e-8 &
        t0_contrasts_cell > 0)
)

# FDRsig Pval10.6sig Pval10.8sig
# 1  12608         178          71
# 2    498         171          68
# 3     52           6           2
# 4  18978       13168        6562
# 5  16379         182          67
# 6      5           0           0
# 7     61           9           3

###################
# load modeling outputs from manual annotations???
# load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb_contrasts.Rdata")
# load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb0_list.Rdata")

ground_truth <- spatialLIBD::fetch_data("modeling_results")

cor_stats_layer <- layer_stat_cor(
    t0_contrasts_cell,
    modeling_results = ground_truth,
    model_type = "enrichment",
    top_n = 100
)

### heatmap ### here can also use layer_stat_cor_plot() from spatialLIBD

## plot output directory
dir_plots <-
    here::here("plots", "07_spatial_registration")
dir.create(dir_plots, showWarnings = FALSE)

# http://research.libd.org/spatialLIBD/reference/layer_stat_cor_plot.html newer function for plotting



pdf(
    file = here::here(
        "plots",
        "07_spatial_registration",
        paste0(
            "dlpfc_pseudobulked_bayesSpace_vs_mannual_annotations_k",
            k,
            ".pdf"
        )
    ),
    width = 8
)
layer_stat_cor_plot_AS(
    cor_stats_layer,
    max = 1,
    cex = 2.5
)

dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
