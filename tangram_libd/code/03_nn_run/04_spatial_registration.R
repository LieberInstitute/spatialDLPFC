## Code adapted from
## https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/dlpfc_snRNAseq_annotation.R

## Code styled with:
# styler::style_file(here::here("tangram_libd", "code", "03_nn_run", "04_spatial_registration.R"), transformers = biocthis::bioc_style())
library("SingleCellExperiment")
library("rafalib")
library("scater")
library("limma")
library("spatialLIBD")
library("here")
library("sessioninfo")

## Define output directories and create them
dir_plots <- here("tangram_libd", "plots", "03_nn_run", "spatial_registration")
dir_rdata <- here("tangram_libd", "processed-data", "03_nn_run", "spatial_registration")
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

## Load Matt's snRNA-seq DLPFC data
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/forAmazonS3/SCE_DLPFC-n3_tran-etal.rda", verbose = TRUE)

## Explore the data a bit
table(sce.dlpfc.tran$donor)
table(sce.dlpfc.tran$processBatch)
sort(table(sce.dlpfc.tran$cellType))
# Inhib_E    Inhib_F      Tcell Macrophage      Mural    Excit_D    Excit_E
#       7          8          9         10         18        132        187
# Excit_F    Inhib_A    Inhib_C      Micro    Inhib_D    Inhib_B    Excit_C
#     243        333        365        388        413        454        524
# Excit_A        OPC    Excit_B      Astro      Oligo
#     529        572        773        782       5455

sort(table(sce.dlpfc.tran$prelimCluster))
#   39   99   28   81   73   40  101   33   76   90  103  105   26   36   79    1
#    4    4    7    7    8    9    9   10   10   10   10   10   11   11   11   12

#   25    3   62   71   74   97   68   80   95   57   82  104   16   43   58   94
#   14   16   16   16   16   17   18   18   19   20   20   20   21   21   21   21

#   84   86   59    2   52   47  102   70   19   41   67   48   60   50   92   98
#   22   22   23   24   24   25   25   26   27   27   27   28   29   30   30   30

#   91   32   42   96   22   46  107   20   30   88   75   87   51   34   54   12
#   31   32   32   33   36   37   38   39   46   47   48   48   49   52   52   53

#   38   44   11   89   24   63   72   85   37   45   10   77   64    9   18   17
#   56   56   57   57   58   58   58   58   59   62   64   64   65   66   66   74

#  106   15   93   31   13  100    6   29   56   49    4    5   83   14    8   55
#   74   78   79   82   87   88   96   97  105  109  110  111  111  117  122  133

#   21   27   61   65   78   23   69    7   53   35   66
#  143  167  191  205  230  358  388  484  654 1610 2666

table(sce.dlpfc.tran$prelimCluster)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
#  12   24   16  110  111   96  484  122   66   64   57   53   87  117   78   21

#  17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32
#  74   66   27   39  143   36  358   58   14   11  167    7   97   46   82   32

#  33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  10   52 1610   11   59   56    4    9   27   32   21   56   62   37   25   28

#  49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64
# 109   30   49   24  654   52  133  105   20   21   23   29  191   16   58   65

#  65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80
# 205 2666   27   18  388   26   16   58    8   16   48   10   64  230   11   18

#  81   82   83   84   85   86   87   88   89   90   91   92   93   94   95   96
#   7   20  111   22   58   22   48   47   57   10   31   30   79   21   19   33

#  97   98   99  100  101  102  103  104  105  106  107
#  17   30    4   88    9   25   10   20   10   74   38
length(table(sce.dlpfc.tran$prelimCluster))

sce.dlpfc.tran$PseudoSample <-
    with(colData(sce.dlpfc.tran), paste0(donor, ":", cellType))
length(table(sce.dlpfc.tran$PseudoSample))
length(unique(sce.dlpfc.tran$cellType))
19 * 3


## sum counts
cIndexes <- splitit(sce.dlpfc.tran$PseudoSample)
umiComb <- sapply(cIndexes, function(ii) {
      rowSums(assays(sce.dlpfc.tran)$counts[, ii, drop = FALSE])
  })

## filter pheno
phenoComb <- colData(sce.dlpfc.tran)[
    !duplicated(sce.dlpfc.tran$PseudoSample),
    c("cellType", "PseudoSample", "donor")
]
rownames(phenoComb) <- phenoComb$PseudoSample
phenoComb <- phenoComb[colnames(umiComb), ]
phenoComb <- DataFrame(phenoComb)
# phenoComb$cellType <- droplevels(phenoComb$cellType)

sce_pseudobulk <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = phenoComb,
        rowData = rowData(sce.dlpfc.tran)
    ))

save(sce_pseudobulk, file = file.path(dir_rdata, "sce_pseudobulk.Rdata"))

## extract expression
mat <- assays(sce_pseudobulk)$logcounts

## Build a group model
mod <- with(
    colData(sce_pseudobulk),
    model.matrix(~ 0 + cellType)
)
colnames(mod) <- gsub("cellType", "", colnames(mod))

## get duplicate correlation
corfit <- duplicateCorrelation(mat, mod,
    block = sce_pseudobulk$donor
)

## Next for each cell type test that specific cell type vs the rest
cell_idx <- splitit(sce_pseudobulk$cellType)

eb0_list_cell <- lapply(cell_idx, function(x) {
    res <- rep(0, ncol(sce_pseudobulk))
    res[x] <- 1
    m <- with(
        colData(sce_pseudobulk),
        model.matrix(~res)
    )
    eBayes(
        lmFit(
            mat,
            design = m,
            block = sce_pseudobulk$donor,
            correlation = corfit$consensus.correlation
        )
    )
})
save(eb0_list_cell, file = file.path(dir_rdata, "dlpfc_snRNAseq_pseudobulked_specific_Ts.Rdata"))

## Extract enrichment t-statistics
t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) <- rowData(sce_pseudobulk)$gene_id
t0_contrasts_cell[1:3, 1:3]
t0_contrasts_cell <- as.data.frame(t0_contrasts_cell, check.names = FALSE)

## Download Nat Neuro 2021 modeling results
modeling_results <- fetch_data(type = "modeling_results")

## Compute the correlations
cor_stats_layer <- layer_stat_cor(
    t0_contrasts_cell,
    modeling_results = modeling_results,
    model_type = "enrichment"
)

## Visualize correlations
pdf(file.path(dir_plots, "spatial_registration_Tran2021.pdf"), width = 11)
layer_stat_cor_plot(cor_stats_layer, max = 0.3)
dev.off()

## Attempt to keep only the top 100 genes per layer
top100 <- unlist(lapply(grep("^t_stat_", colnames(modeling_results$enrichment)), function(i) {
    order(modeling_results$enrichment[, i], decreasing = TRUE)[1:100]
}))

modeling_results_top100 <- modeling_results
modeling_results_top100$enrichment <- modeling_results_top100$enrichment[top100, ]

cor_stats_layer_top100 <- layer_stat_cor(
    t0_contrasts_cell,
    modeling_results = modeling_results_top100,
    model_type = "enrichment"
)

## Visualize correlations
pdf(file.path(dir_plots, "spatial_registration_Tran2021_top100LayerGenes.pdf"), width = 11)
## The max is typically oligo/WM, so let's use the 2nd max
x <- as.numeric(cor_stats_layer_top100)
layer_stat_cor_plot(cor_stats_layer_top100, max = max(x[x != max(x)]))
dev.off()
## We can see good/decent layer-specificity for the excitatory neurons but
## not for the inhibitory neurons. Combined with Astro and Oligo, we get
## all layers and WM, which could be all we need for spatially registering
## this snRNA-seq data, which then when applied via Tangram to the new
## spatialDLPFC data could help us identify the layers in this new spatial
## DLPFC dataset.



### Do it again but at the prelimCluster level in case they are more
### layer-specific than the cellType's from Matt.

sce.dlpfc.tran$PseudoSample <-
    with(colData(sce.dlpfc.tran), paste0(donor, ":", prelimCluster))


## sum counts
cIndexes <- splitit(sce.dlpfc.tran$PseudoSample)
umiComb <- sapply(cIndexes, function(ii) {
      rowSums(assays(sce.dlpfc.tran)$counts[, ii, drop = FALSE])
  })

## filter pheno
phenoComb <- colData(sce.dlpfc.tran)[
    !duplicated(sce.dlpfc.tran$PseudoSample),
    c("prelimCluster", "PseudoSample", "donor")
]
rownames(phenoComb) <- phenoComb$PseudoSample
phenoComb <- phenoComb[colnames(umiComb), ]
phenoComb <- DataFrame(phenoComb)

sce_pseudobulk <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = phenoComb,
        rowData = rowData(sce.dlpfc.tran)
    ))

most_common <- sapply(cIndexes, function(i) {
    names(sort(table(sce.dlpfc.tran$cellType[i]), decreasing = TRUE))[1]
})
stopifnot(identical(names(most_common), sce_pseudobulk$PseudoSample))
sce_pseudobulk$most_common <- most_common

save(sce_pseudobulk, file = file.path(dir_rdata, "sce_pseudobulk_prelimCluster.Rdata"))

## extract expression
mat <- assays(sce_pseudobulk)$logcounts

## Build a group model
mod <- with(
    colData(sce_pseudobulk),
    model.matrix(~ 0 + prelimCluster)
)
colnames(mod) <- gsub("prelimCluster", "", colnames(mod))

## get duplicate correlation
corfit <- duplicateCorrelation(mat, mod,
    block = sce_pseudobulk$donor
)

## Next for each cell type test that specific cell type vs the rest
cell_idx <- splitit(sce_pseudobulk$prelimCluster)

eb0_list_cell <- lapply(cell_idx, function(x) {
    res <- rep(0, ncol(sce_pseudobulk))
    res[x] <- 1
    m <- with(
        colData(sce_pseudobulk),
        model.matrix(~res)
    )
    eBayes(
        lmFit(
            mat,
            design = m,
            block = sce_pseudobulk$donor,
            correlation = corfit$consensus.correlation
        )
    )
})
save(eb0_list_cell, file = file.path(dir_rdata, "dlpfc_snRNAseq_pseudobulked_specific_Ts_prelimCluster.Rdata"))

## Extract enrichment t-statistics
t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) <- rowData(sce_pseudobulk)$gene_id
t0_contrasts_cell <- as.data.frame(t0_contrasts_cell, check.names = FALSE)

cor_stats_layer_top100 <- layer_stat_cor(
    t0_contrasts_cell,
    modeling_results = modeling_results_top100,
    model_type = "enrichment"
)

m <- match(rownames(cor_stats_layer_top100), as.character(sce_pseudobulk$prelimCluster))
stopifnot(all(!is.na(m)))
rownames(cor_stats_layer_top100) <- paste0(rownames(cor_stats_layer_top100), "_", sce_pseudobulk$most_common[m])

## Visualize correlations
pdf(file.path(dir_plots, "spatial_registration_Tran2021_top100LayerGenes_prelimCluster.pdf"), width = 20)
## The max is typically oligo/WM, so let's use the 2nd max
x <- as.numeric(cor_stats_layer_top100)
layer_stat_cor_plot(cor_stats_layer_top100, max = max(x[x != max(x)]))
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

