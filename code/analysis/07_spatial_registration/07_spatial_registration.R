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

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load spe object
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)

#load clusters
spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"), 
  prefix = ""
)


#pseudobulk my data based on bayesSpace clustering and sample
spe$PseudoSample = paste0(spe$sample_id, ":",colData(spe)[[paste0("bayesSpace_harmony_",k)]])
cIndexes = splitit(spe$PseudoSample) # gives you the index for each cluster:sample

#sum umis for each pseudobulked group (cluster). produces a data frame where the rows are genes and the columns are the pseudobulked samples:clusters
# and the values are the total number of counts for each gene in each sample:cluster
umiComb <- sapply(cIndexes, function(ii)
  rowSums(assays(spe)$counts[, ii, drop = FALSE])) #makes data frame of count sums where gene is row and sample:cluster is column

phenoComb = colData(spe)[!duplicated(spe$PseudoSample),] #creates new colData dataframe with pseudobulked colData
rownames(phenoComb) = phenoComb$PseudoSample #renames rows of new colData frame to be the cluster:samples
phenoComb = phenoComb[colnames(umiComb), ]
phenoComb = DataFrame(phenoComb)

sce_pseudobulk_bayesSpace <-
  logNormCounts(SingleCellExperiment(
    list(counts = umiComb),
    colData = phenoComb,
    rowData = rowData(spe)
  ))

save(sce_pseudobulk_bayesSpace, file = here::here("processed-data","rdata","spe","07_spatial_registration",paste0("sce_pseudobulk_bayesSpace_k",k,".Rdata")))

###############################
##### get mean expression  ####
mat <- assays(sce_pseudobulk_bayesSpace)$logcounts #make matrix of just the log normalized counts

## filter
gIndex = rowMeans(mat) > 0.2 # find the genes for which the mean expression is greater than 0.2
mat_filter = mat[gIndex, ] #subset matrix on just those genes.  want to remove lowly expressed genes. 

#####################
## Build a group model

#convert variables to factors 
colData(sce_pseudobulk_bayesSpace)$spatial.cluster <- as.factor(colData(sce_pseudobulk_bayesSpace)[[paste0("bayesSpace_harmony_",k)]])
colData(sce_pseudobulk_bayesSpace)$region <- as.factor(colData(sce_pseudobulk_bayesSpace)$region)
colData(sce_pseudobulk_bayesSpace)$age <- as.integer(colData(sce_pseudobulk_bayesSpace)$age)
colData(sce_pseudobulk_bayesSpace)$sex <- as.factor(colData(sce_pseudobulk_bayesSpace)$sex)
colData(sce_pseudobulk_bayesSpace)$diagnosis <- as.factor(colData(sce_pseudobulk_bayesSpace)$diagnosis)
colData(sce_pseudobulk_bayesSpace)$subject <- as.factor(colData(sce_pseudobulk_bayesSpace)$subject)

#create matrix where the rownames are the sample:clusters and the columns are the other variales (spatial.cluster + region + age + sex)
mod <- with(colData(sce_pseudobulk_bayesSpace),
            model.matrix(~ 0 + spatial.cluster + region + age + sex)) #removed diagnosis cuz it only has 1 level 
colnames(mod) <- gsub('cluster', '', colnames(mod)) 
#why is regtionanterior missing in the colnames(mod)?

## get duplicate correlation #http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/dupcor.html
corfit <- duplicateCorrelation(mat_filter, mod,
                               block = sce_pseudobulk_bayesSpace$subject)
save(corfit, file = here::here("processed-data","rdata","spe","07_spatial_registration",paste0("dlpfc_pseudobulked_bayesSpace_dupCor_k",k,".Rdata")))

## Next for each layer test that layer vs the rest
cell_idx <- splitit(sce_pseudobulk_bayesSpace$spatial.cluster)

eb0_list_cell <- lapply(cell_idx, function(x) {
  res <- rep(0, ncol(sce_pseudobulk_bayesSpace))
  res[x] <- 1
  m <- with(colData(sce_pseudobulk_bayesSpace),
            model.matrix(~ res +
                           region + age + sex))
  eBayes(
    lmFit(
      mat_filter,
      design = m,
      block = sce_pseudobulk_bayesSpace$subject,
      correlation = corfit$consensus.correlation
    )
  )
})
save(eb0_list_cell, file = here::here("processed-data","rdata","spe","07_spatial_registration",paste0("dlpfc_pseudobulked_bayesSpace_specific_Ts_k",k,".Rdata")))


##########
## Extract the p-values

pvals0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cell) = rownames(mat_filter)

t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) = rownames(mat_filter)
fdrs0_contrasts_cell = apply(pvals0_contrasts_cell, 2, p.adjust, 'fdr')

data.frame(
  'FDRsig' = colSums(fdrs0_contrasts_cell < 0.05 &
                       t0_contrasts_cell > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts_cell < 1e-6 &
                            t0_contrasts_cell > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts_cell < 1e-8 &
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
#load modeling outputs from manual annotations???
load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb_contrasts.Rdata")
load("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/rda/eb0_list.Rdata")

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) = rownames(eb_contrasts)
fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) = rownames(eb_contrasts)

############
# line up ## ##from here to line 175 I supposed to use the function leo created for spatialLIBD called layer_stat_cor()

mm = match(rownames(pvals0_contrasts), rownames(pvals0_contrasts_cell))

pvals0_contrasts = pvals0_contrasts[!is.na(mm), ]
t0_contrasts = t0_contrasts[!is.na(mm), ]
fdrs0_contrasts = fdrs0_contrasts[!is.na(mm), ]

pvals0_contrasts_cell = pvals0_contrasts_cell[mm[!is.na(mm)], ]
t0_contrasts_cell = t0_contrasts_cell[mm[!is.na(mm)], ]
fdrs0_contrasts_cell = fdrs0_contrasts_cell[mm[!is.na(mm)], ]

cor_t = cor(t0_contrasts_cell, t0_contrasts)
signif(cor_t, 2)

# WM Layer1  Layer2 Layer3 Layer4 Layer5 Layer6
# 1 -0.420 -0.150  0.5000  0.420   0.21  0.046 -0.071
# 2  0.630 -0.031 -0.3600 -0.530  -0.32 -0.230  0.096
# 3  0.099  0.370 -0.1400  0.059  -0.18 -0.190 -0.200
# 4 -0.340 -0.350  0.0760  0.180   0.43  0.390  0.047
# 5 -0.150 -0.380  0.0360 -0.081   0.14  0.290  0.390
# 6 -0.095 -0.200  0.0480  0.056   0.12  0.130  0.110
# 7  0.140  0.520 -0.0019 -0.093  -0.26 -0.310 -0.240

### just layer specific genes from ones left
layer_specific_indices = mapply(function(t, p) {
  oo = order(t, decreasing = TRUE)[1:100]
},
as.data.frame(t0_contrasts),
as.data.frame(pvals0_contrasts))
layer_ind = unique(as.numeric(layer_specific_indices))

cor_t_layer = cor(t0_contrasts_cell[layer_ind, ],
                  t0_contrasts[layer_ind, ])
signif(cor_t_layer, 3)

### heatmap ### here can also use layer_stat_cor_plot() from spatialLIBD
theSeq = seq(-.85, .85, by = 0.01)
my.col <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

dd = dist(1-cor_t_layer)
hc = hclust(dd)
cor_t_layer_toPlot = cor_t_layer[hc$order, c(1, 7:2)]
colnames(cor_t_layer_toPlot) = gsub("ayer", "", colnames(cor_t_layer_toPlot))
rownames(cor_t_layer_toPlot)[rownames(cor_t_layer_toPlot) == "Oligodendrocytes"] = "OLIGO" # does thismatter? 

pdf(file = here::here("plots","07_spatial_registration",paste0("dlpfc_pseudobulked_bayesSpace_vs_mannual_annotations_k",k,".pdf")), width = 8)
print(
  levelplot(
    cor_t_layer_toPlot,
    aspect = "fill",
    at = theSeq,
    col.regions = my.col,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 1.5), y = list(cex = 1.5))
  )
)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()