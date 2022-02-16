###########################
#fasthplus, need to use hpb()
##########################

#install_github(repo="ntdyjack/fasthplus", ref = "main")
library(fasthplus)
library(SpatialExperiment)
library(here)

#load spe object
load(file = here::here("processed-data","rdata","spe","spe_final.Rdata"))

#hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations
hpb(D= reducedDims(spe)$HARMONY,L=as.vector(spe$spatial.cluster),t=10,r=10)

## error is line 30 of bsf which is called by hpb() in line 49