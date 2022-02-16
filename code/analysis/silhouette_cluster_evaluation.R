# qrsh -l bluejay,mem_free=75G,h_vmem=75G -now n

library(SpatialExperiment)
library(here)

# Clustering
library(cluster) 
library(factoextra)

# adapted from: https://medium.com/codesmart/r-series-k-means-clustering-silhouette-794774b46586
# silhouette_score <- function(k){
#   km <- kmeans(df, centers = k, nstart=25)
#   ss <- silhouette(km$cluster, dist(df))
#   mean(ss[, 3])
# }
# k <- 2:10
# avg_sil <- sapply(k, silhouette_score)
# plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)


# load SPE object 
load(file = here::here("processed-data","rdata","spe","spe_final.Rdata"))

#import all BayesSpace cluster
cluster_export(
  spe,
  bayesSpace_name,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results","bayesSpace_k" )
)

clusterColnames <- paste0("bayesSpace_harmony_",4:15)
for(i in seq_along(clusterColnames)){
  ss <- silhouette(spe$clusterColnames[i], dist(reducedDims(spe)$HARMONY))
}

#1:10pm
ss <- silhouette(spe$spatial.cluster, dist(reducedDims(spe)$HARMONY))

