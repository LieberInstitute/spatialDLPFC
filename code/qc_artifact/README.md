Scripts used to evaluate the effect of artifacts on downstream analysis and produce Supplementary Figures 5, 6 and 7. The scripts listed below should be run in the order listed

1. Gene_based_analysis_artifact.R: 
  - Compute and compare QC metrics - Library size, Number of detected genes, Percent mitochondrial gene expression, and Number of nuclei/ spot between artifact and non-artifact spots. 
  - Fit the observed gene expression using negative binomial error with either only library size as the independent effect or with library size and artifact. Generates Supplementary Fig 5 A-C
  
2. Spot_based_analysis_artifact.R: 
  - Normalize and compare properties of gene expression with and without artifacts. 
  - Learn lower dimensional embedding of all spots using UMAP to visualize clustering of artifact spots. 
  - Using ANOVA determine the drivers of variance in the top 10 PCs. 
  - Cluster samples using SNN and examine the composition of clusters by layers and artifacts. 
  - Perform differential gene expression testing with T-tests within each sample and layer. (Generates Supplementary Fig 5 D, Fig 6)

3. Simulate_heterotypic_artifacts.R:
 - Simulate heterotypic artifacts and normalize pooled observed and simulated data
 - Perform PCA and UMAP embedding of combined data
 - Find the 10 nearest neighbors for each spot using SNN. Compute the mean and variance of the identity of the neighbors. (Generates Supplementary Fig 7 A, B)
 
4. Compare_cellType_correlation_diff.R: This script compares the correlation between pairs of 7 unique cell types for all spots and spots excluding artifacts and compares them (Supplementary Fig 7C)