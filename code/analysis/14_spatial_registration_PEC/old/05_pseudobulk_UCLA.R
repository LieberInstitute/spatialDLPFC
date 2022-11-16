library("zellkonverter")
library(SingleCellExperiment)
library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)
library(edgeR)

#######################
### load UCLA dataset
########################
sce <- readH5AD(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version2/UCLA-ASD/UCLA-ASD-snRNAseq_annotated.h5ad")

# identify annotation/cluster labels
rowData(sce)$gene_name <- rownames(rowData(sce)) # save gene name as column of rowData
rownames(sce) <- rowData(sce)$featureid # have to make row names of object the ensembl id instead of gene names
colData(sce)$anno <- as.factor(colData(sce)$anno)
colnames(colData(sce))[7] <- "sample_id"
names(assays(sce)) <- "counts"

# save sce here because it takes so long to convert the h5AD file
saveRDS(sce, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/UCLA_sce.rds")
sce$anno <- droplevels(sce$anno)

message("Make psuedobulk object")
spe_pseudo <- aggregateAcrossCells(
    sce,
    DataFrame(
        Anno = sce$anno,
        sample_id = sce$sample_id
    )
)
colnames(spe_pseudo) <- paste0(spe_pseudo$sample_id, "_", spe_pseudo$anno)

message("Filter lowly expressed genes")
rowData(spe_pseudo)$filter_expr <- filterByExpr(spe_pseudo)
summary(rowData(spe_pseudo)$filter_expr)

spe_pseudo <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr), ]

message("Normalize expression")
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
## We don't need this 'x' object anymore
rm(x)

# save pseudobulked object
saveRDS(spe_pseudo, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/UCLA_spe_pseudo.rds")

mat <- assays(spe_pseudo)$logcounts

mod <- with(
    colData(spe_pseudo),
    model.matrix(~ 0 + anno)
)
colnames(mod) <- gsub("anno", "", colnames(mod))

corfit <- duplicateCorrelation(mat, mod,
    block = spe_pseudo$sample_id
)

cluster_idx <- splitit(spe_pseudo$anno)

message("Run enrichment statistics")
eb0_list_cluster <- lapply(cluster_idx, function(x) {
    res <- rep(0, ncol(spe_pseudo))
    res[x] <- 1
    m <- with(
        colData(spe_pseudo),
        model.matrix(~res)
    )
    eBayes(
        lmFit(
            mat,
            design = m,
            block = spe_pseudo$sample_id,
            correlation = corfit$consensus.correlation
        )
    )
})

message("extract and reformat enrichment results")
##########
## Extract the p-values

pvals0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cluster) <- rownames(mat)

t0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cluster) <- rownames(mat)
fdrs0_contrasts_cluster <- apply(pvals0_contrasts_cluster, 2, p.adjust, "fdr")

data.frame(
    "FDRsig" = colSums(fdrs0_contrasts_cluster < 0.05 &
        t0_contrasts_cluster > 0),
    "Pval10-6sig" = colSums(pvals0_contrasts_cluster < 1e-6 &
        t0_contrasts_cluster > 0),
    "Pval10-8sig" = colSums(pvals0_contrasts_cluster < 1e-8 &
        t0_contrasts_cluster > 0)
)

save(pvals0_contrasts_cluster,
    t0_contrasts_cluster,
    fdrs0_contrasts_cluster,
    file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/UCLA_modeling_results.rda"
)

# vs manual annotations
modeling_results <- fetch_data(type = "modeling_results")
cor <- layer_stat_cor(
    t0_contrasts_cluster,
    modeling_results,
    model_type = names(modeling_results)[2],
    reverse = FALSE,
    top_n = NULL
)


pdf("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/14_spatial_registration_PEC/UCLA/spatial_registration_plot_UCLA_v_manual.pdf")
layer_stat_cor_plot(cor)
dev.off()

# load my k = 9 modeling results
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/07_layer_differential_expression/parsed_modeling_results_k9.Rdata")
cor <- layer_stat_cor(
    t0_contrasts_cluster,
    modeling_results,
    model_type = names(modeling_results)[2],
    reverse = FALSE,
    top_n = NULL
)
pdf("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/14_spatial_registration_PEC/UCLA/spatial_registration_plot_UCLA_v_k9.pdf")
layer_stat_cor_plot(cor)
dev.off()

# load my k = 16 modeling results
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/07_layer_differential_expression/parsed_modeling_results_k16.Rdata")
cor <- layer_stat_cor(
    t0_contrasts_cluster,
    modeling_results,
    model_type = names(modeling_results)[2],
    reverse = FALSE,
    top_n = NULL
)
pdf("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/14_spatial_registration_PEC/UCLA/spatial_registration_plot_UCLA_v_k16.pdf")
layer_stat_cor_plot(cor)
dev.off()
