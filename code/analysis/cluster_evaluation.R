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


