# https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/asd_snRNAseq_recast.R
# https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/ad_snRNAseq_recast.R
# http://research.libd.org/spatialLIBD/reference/layer_stat_cor.html
# http://research.libd.org/spatialLIBD/reference/layer_stat_cor_plot.html

library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)

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

