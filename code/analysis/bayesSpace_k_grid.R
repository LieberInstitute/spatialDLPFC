library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library(BayesSpace)
library("RColorBrewer")

load(file = here::here("processed-data","rdata","spe","spe_final.Rdata"),verbose = TRUE)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

set.seed(20220127)

spe = spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 10000) 

spe$bayesSpace_temp<-spe$spatial.cluster
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
  spe,
  bayesSpace_name,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results" )
)

sample_ids <- unique(colData(spe)$sample_id)
mycolors <- brewer.pal(7, "Dark2")

pdf(file = here::here("plots",paste0("vis_clus_bayesSpace_harmony_",k,".pdf")))
for (i in seq_along(sample_ids)){
  my_plot <- vis_clus(
    spe = spe,
    clustervar = bayesSpace_name,
    sampleid = sample_ids[i],
    colors =  mycolors
  )
  print(my_plot)
}
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

#sgejobs::job_single("bayesSpace_k_grid",create_shell = TRUE, queue = "bluejay", memory = "80G", task_num = 15)

