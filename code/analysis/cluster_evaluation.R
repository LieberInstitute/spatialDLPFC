library("SpatialExperiment")
library("mclust")
library("spatialLIBD")

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

spe<- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","graph_based_harmony"),
  prefix = "harmony_graph_based_"
)

spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","semi_supervised_harmony_across_samples"),
  prefix = "harmony_semi_supervised_"
)

spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","graph_based_pcs"),
  prefix = "pcs_graph_based_"
)


adjustedRandIndex(spe$spatial.cluster,spe$harmony_semi_supervised_SNN_k10_k7)
# 0.1958784

adjustedRandIndex(spe$spatial.cluster,spe$harmony_graph_based_SNN_k10_k7)
# 0.1958784

adjustedRandIndex(spe$spatial.cluster,spe$pcs_graph_based_SNN_k10_k7)
# 0.213768

Method <- c("Graph-Based PCs", "Graph-Based Batch Corrected", "Semi-supervised Batch Corrected")
ARI<- c(0.213768, 0.1958784, 0.1958784)

df <- data.frame(Method, ARI)

print (df)

pdf(here::here("plots","my_data_ARI.pdf"))
ggplot(df, aes(x=Method, y=ARI)) +
  geom_point(size=2) +
  theme_bw()
dev.off()


#### spaGCN ARI for pilot data
spe.temp <- cluster_import(
  spe,
  cluster_dir = here::here("..","spython","spagcn","processed-data","03-our_data_analysis"),
  prefix = "spaGCN_"
)


with(colData(spe),addmargins(table(spatial.cluster,pseudobulk_PCA.y,sample_id)))


