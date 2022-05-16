load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")

for (i in 4:28){
  print(i)
  print(length(unique(colData(spe)[[paste0("bayesSpace_harmony_",i)]])))
}

#k = 24 is missing a cluster
#k = 25 is missing a cluster
#k = 26 is missing a cluster
#k = 27 is missing a cluster
#k = 28 is missing two clusters