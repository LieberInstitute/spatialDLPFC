library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)


#sample_names <- c("DLPFC_Br2743_ant_manual_alignment", "DLPFC_Br2743_mid_manual_alignment","DLPFC_Br2743_post_manual_alignment","DLPFC_Br3942_ant_manual_alignment","DLPFC_Br3942_mid_manual_alignment","DLPFC_Br3942_post_manual_alignment","DLPFC_Br6423_ant_manual_alignment","DLPFC_Br6423_mid_manual_alignment","DLPFC_Br6423_post_manual_alignment","DLPFC_Br8492_ant_manual_alignment","DLPFC_Br8492_mid_manual_alignment","DLPFC_Br8492_post_manual_alignment")

#for (i in 7:8) {
  sample_name <-"DLPFC_Br2743_mid_manual_alignment"
  #sample_name <-"DLPFC_Br3942_post_manual_alignment"
  #sample_name <- sample_names[i]
  visium_dir <- paste0("/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/outputs/NextSeq/",sample_name,"/outs/")
  #read data into sce object
  sce <- readVisium(visium_dir)
  
  #pre-processing
  set.seed(100)
  dlpfc <- spatialPreprocess(sce, platform = "Visium")
  
  #selecting the number of clusters
  dlpfc <- qTune(dlpfc, qs=seq(2, 10), platform="Visium", d=7)
  
  pdf(file = paste0("/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/plots/BayesSpace/",sample_name, "BayesSpace_likelihood.pdf"))
  #pdf(file = paste0("/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/plots/BayesSpace/DLPFC_Br2743_ant_manual_alignment/BayesSpace_likelihood.pdf"))
  qPlot(dlpfc)
  dev.off()
  
  # using smaller number of MCMC iterations for faster runtime
  set.seed(100)
  dlpfc <- spatialCluster(dlpfc, q = 7, platform = "Visium", 
                          nrep = 10000, burn.in = 500, save.chain = TRUE)
  
  head(colData(dlpfc))
  
  pdf(file = paste0("/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/plots/BayesSpace/",sample_name, "BayesSpace_clusterPlot.pdf"))
  clusterPlot(dlpfc)
  dev.off()
  
  #dlpfc <- readRDS("/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/plots/BayesSpace/dlpfc_Br2743_ant_dlpfc.RDS")
  #enhanced resolution
  dlpfc_enhanced <- spatialEnhance(dlpfc, q = 7, platform = "Visium", 
                                   nrep = 100000, burn.in = 500, save.chain = TRUE)
  #clusterPlot(dlpfc_enhanced)
  
  head(colData(dlpfc_enhanced))
  
  #featurePlot(dlpfc_enhanced)
  
  
  pdf(file = paste0("/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/plots/BayesSpace/", sample_name, "BayesSpace_enhanced_clusterPlot.pdf"))
  clusterPlot(dlpfc_enhanced)
  dev.off()
#}
