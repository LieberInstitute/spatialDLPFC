library("SpatialExperiment")
library("mclust")

load(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/pilot_dlpfc_data/spe_pilot_102121.Rdata")

table_percent <- function(input.table){
  list(
    "counts" = addmargins(input.table),
    
    "percent_all" = addmargins(round(input.table/sum(input.table)*100,2)), #divide my all spots
    
    "percent_col" = addmargins(round(sweep(input.table,2,colSums(input.table),"/")*100,2)), #divide by spots in bayesSpace cluster
    
    "percent_row" = addmargins(round(sweep(input.table,1,rowSums(input.table),"/")*100,2)) #divide by spots in Kristen's annotation
    
  )
  
}

options(width=120)

table.layer <- with(colData(spe), table(
  layer_guess_reordered_short, spatial.cluster
))

table_percent(table.layer)

adjustedRandIndex(spe$spatial.cluster,spe$layer_guess_reordered_short)

table.subject <- with(colData(spe), table(
  subject, spatial.cluster
))
table_percent(table.subject)

adjustedRandIndex(spe$spatial.cluster,spe$subject)


table.layer <- with(colData(spe), table(
  layer_guess_reordered_short, batch_corr_SNN_k10_k7
))

table_percent(table.layer)

adjustedRandIndex(spe$batch_corr_SNN_k10_k7,spe$layer_guess_reordered_short)

table.layer <- with(colData(spe), table(
  layer_guess_reordered_short, SNN_k50_k7
))

table_percent(table.layer)

adjustedRandIndex(spe$SNN_k50_k7,spe$layer_guess_reordered_short)

#plot ARI for pilot data comparing different clustering algorithms to Kristen's manual annnotations
library(ggplot2)
Method <- c("BayesSpace", "BayesSpace Batch Corrected", "Graph-Based", "Graph-Based Batch Corrected")
ARI<- c(0.21, 0.34, 0.13, 0.16)

df <- data.frame(Method, ARI)

print (df)

pdf(here::here("plots","pilot_ARI.pdf"))
  ggplot(df, aes(x=Method, y=ARI)) +
    geom_point(size=2) +
    theme_bw()
dev.off()

#plot ARI for my data comparing different clustering methods against BayesSpace batch corrected
cluster_colNames <- paste0("SNN_k10_k",4:28)
for (i in seq_along(cluster_colNames)){
  colData(spe) <- cbind(colData(spe),clust_k5_list[i])
}
colnames(colData(spe))[37:61] <- cluster_colNamesload("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/graph_based_harmony/clust_k5_list_harmony.Rdata")

adjustedRandIndex(spe$spatial.cluster,spe$SNN_k10_k7)
# [1] 0.1958784

adjustedRandIndex(spe$spatial.cluster,spe$pseudobulk_PCA.y)
#[1] 0.1958784

adjustedRandIndex(spe$spatial.cluster,spe$cluster.init) 
#[1] 0.3465709
adjustedRandIndex(spe$spatial.cluster,spe$pseudobulk_PCA.y)  
#[1] 0.1958784
adjustedRandIndex(spe$spatial.cluster,spe$SNN_k10_k24)
#[1] 0.2009123
adjustedRandIndex(spe$spatial.cluster,spe$SNN_k10_k8)
#[1] 0.1908048

