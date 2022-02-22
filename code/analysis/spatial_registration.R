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

#load spe object
load(file = here::here("processed-data","rdata","spe","spe_final.Rdata"),verbose = TRUE)

#pseudobulk my data
spe$PseudoSample = paste0(spe$sample_id, ":", spe$spatial.cluster)
cIndexes = splitit(spe$PseudoSample) # gives you the index for each cluster 

#sum umis for each pseudobulked group (cluster). produces a data frame where the rows are genes and the columns are the pseudobulked samples (clusters)
# and the values are the total number of counts for each gene in each cluster
umiComb <- sapply(cIndexes, function(ii)
  rowSums(assays(spe)$counts[, ii, drop = FALSE])) #

phenoComb = colData(spe)[!duplicated(spe$PseudoSample),] #creates new colData dataframe with pseudobulked colData
rownames(phenoComb) = phenoComb$PseudoSample #renames rows of new colData frame to be the clusters
phenoComb = phenoComb[colnames(umiComb), ]
phenoComb = DataFrame(phenoComb)

sce_pseudobulk_bayesSpace_k7 <-
  logNormCounts(SingleCellExperiment(
    list(counts = umiComb),
    colData = phenoComb,
    rowData = rowData(spe)
  ))

save(sce_pseudobulk_bayesSpace_k7, file = here::here("processed-data","rdata","spe","sce_pseudobulk_bayesSpace_k7.Rdata"))

###############################
##### get mean expression  ####
mat <- assays(sce_pseudobulk_bayesSpace_k7)$logcounts

## filter
gIndex = rowMeans(mat) > 0.2 # find the genes for which the mean expression is greater than 0.2
mat_filter = mat[gIndex, ] #subset matrix on just those genes.  want to remove lowly expressed genes. 

#####################
## Build a group model

#convert variables to factors 
colData(sce_pseudobulk_bayesSpace_k7)$spatial.cluster <- as.factor(colData(sce_pseudobulk_bayesSpace_k7)$spatial.cluster)
colData(sce_pseudobulk_bayesSpace_k7)$region <- as.factor(colData(sce_pseudobulk_bayesSpace_k7)$region)
colData(sce_pseudobulk_bayesSpace_k7)$age <- as.integer(colData(sce_pseudobulk_bayesSpace_k7)$age)
colData(sce_pseudobulk_bayesSpace_k7)$sex <- as.factor(colData(sce_pseudobulk_bayesSpace_k7)$sex)
colData(sce_pseudobulk_bayesSpace_k7)$diagnosis <- as.factor(colData(sce_pseudobulk_bayesSpace_k7)$diagnosis)
colData(sce_pseudobulk_bayesSpace_k7)$subject <- as.factor(colData(sce_pseudobulk_bayesSpace_k7)$subject)


mod <- with(colData(sce_pseudobulk_bayesSpace_k7),
            model.matrix(~ 0 + spatial.cluster + region + age + sex)) #removed diagnosis cuz it only has 1 level 
colnames(mod) <- gsub('cluster', '', colnames(mod)) #not neccesary 

## get duplicate correlation
corfit <- duplicateCorrelation(mat_filter, mod,
                               block = sce_pseudobulk_bayesSpace_k7$subject)
save(corfit, file = here::here("processed-data","rdata","spe","dlpfc_pseudobulked_bayesSpace_k7_dupCor.Rdata"))

## Next for each layer test that layer vs the rest
cell_idx <- splitit(sce_pseudobulk_bayesSpace_k7$spatial.cluster)

eb0_list_cell <- lapply(cell_idx, function(x) {
  res <- rep(0, ncol(sce_pseudobulk_bayesSpace_k7))
  res[x] <- 1
  m <- with(colData(sce_pseudobulk_bayesSpace_k7),
            model.matrix(~ res +
                           region + age + sex))
  eBayes(
    lmFit(
      mat_filter,
      design = m,
      block = sce_pseudobulk_bayesSpace_k7$subject,
      correlation = corfit$consensus.correlation
    )
  )
})
save(eb0_list_cell, file = here::here("processed-data","rdata","spe","dlpfc_pseudobulked_bayesSpace_k7_specific_Ts.Rdata"))







