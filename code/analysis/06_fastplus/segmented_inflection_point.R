library(SpatialExperiment)
library(here)
library("sessioninfo")
library(spatialLIBD)
library(segmented)

#load spe object
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)

fasthplus <- read.csv(file = here::here("processed-data","rdata","spe","06_fasthplus","fasthplus_results_no_WM.csv"))

f2 <- lm(X1.h ~ k, data = fasthplus)

seg2 <- segmented(f2,
                  seg.Z = ~k,
                  npsi=2
)

# Call: segmented.lm(obj = f2, seg.Z = ~k, npsi = 2)
# 
# Meaningful coefficients of the linear terms:
#   (Intercept)            k         U1.k         U2.k  
# 0.863357    -0.052563     0.055860    -0.002778  
# 
# Estimated Break-Point(s):
#   psi1.k  psi2.k  
# 4.497  16.296  