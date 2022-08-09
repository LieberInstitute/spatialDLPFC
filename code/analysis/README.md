# Code


### 01_build_spe/
`01_build_spe.R`
Read 10x data and create SpatialExperiment object. Identify mito genes. Add cell segementation. Filter genes with filterByExpr() and spots with low library size with scran. Find highly variable genes. Do dimensionality reduction. Batch correct with Harmony. 

`01_build_spe.sh`

`sfigu_harmony.R`
Create supplemental figure showing effects of batch correction. 

### 01a_marker_genes/
`01_marker_genes.R`
Make spot plots of well-known marker genes and new marker genes. 

`01_marker_genes.sh`

### 02_graph_based_clustering/
`01_graph_based_clustering.R`
Perform graph-based clustering within and across samples. 

`01_graph_based_clustering.sh`

### 03_BayesSpace
`01_BayesSpace.R`
Use BayesSpace to cluster across all samples for k = 2 to k = 28. 

`01_BayesSpace.sh`
Array job for k = 2 to k = 15

`01a_BayesSpace_k_16_28.sh`
Array job for k = 16 to k = 28

`02_cluster_check.R`
Noticed at some values of k there were missing cluster labels.

`03_BayesSpace_big_plots.R`
Vis_clus plots of BayesSpace clusters for k = 2 to k = 28. Made custom version of vis_clus() function. 

`03_BayesSpace_big_plots.sh`
Array job for k = 2 to k = 28

`04_BayesSpace_grid_plots.R`
custom version of vis_grid_clus() plots. 

`04_BayesSpace_grid_plots.sh`

`05_BayesSpace_sfigu_plots.R`
Create k9, k16, and k28 with custom vis_grid_clus() plots for supplemental figure of manuscript. 

### 04_semi_supervised_clustering
`01_semi_supervised_clustering.R`
Run semi-supervised clustering with genes from pilot dataset.

`01_semi_supervised_clustering.sh`

### 05_ARI
`01_ARI.R`
Do ARI calculations comparing BayesSpace k=7 clustering of pilot data to manual annotations of pilot data

### 06_fastplus
`01_fasthplus.R`
Run fasthplus discordance internal validity metric for BayesSpace clustering k = 2 to k = 28. 

`01_fasthplus.sh`
Array job k = 3 to k = 15

`01a_fasthplus_16_23.sh`
Array job k = 16 to k = 23

`01b_fasthplus_24_28.sh`
Array job k = 24 to k = 28

`02_segmented_inflection_point.R`
Find inflection points of fasthplus graph

`03_plot_fasthplus.R`
Make plot of fasthplus measure for manuscript

### 07_spatial_registration
`01_spatial_registration.R`
Perform spatial registration on BayesSpace clusters for k = 2 to k = 28 against manual annotations.

`01_spatial_registration.sh`

`01a_spatial_registration_16_23.sh`

`01b_spatial_registration_24_28.sh`

`02_bad_clusters.R`
Originally created script to investigate BayesSpace clusters that didn't look like any of the manual annotations according to spatial registration. Then made heatmaps of the number of spots from each BayesSpace at k = 9 fall into which of the BayesSpace clusters at k = 16. Also did this with k = 9 vs k = 28. This is a supplementary figure called spots_in_cluster.

`03_custom_spatial_registration_plots.R`
Make custom versions of spatial registration plots for manuscript. 

`03_custom_spatial_registration_plots.sh`

`04_color_bar.R`
Make color bar to go under spatial registration plots of figure 2.

`layer_matrix_plot_AS.R`
My version of layer_matrix_plot() function.

`layer_stat_cor_plot_AS.R`
My version of layer_stat_cor_plot() function.

`05_t_stat_cor_plot.R`
Make t-stat correlation plot for figure 2. 








