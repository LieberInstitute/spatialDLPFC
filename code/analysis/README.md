# Code

`.sh` files run corresponding Rscripts, some have additional details noted below

### 01_build_spe/
`01_build_spe.R`
Read 10x data and create SpatialExperiment object. Identify mito genes. Add cell segementation. Filter genes with filterByExpr() and spots with low library size with scran. Find highly variable genes. Do dimensionality reduction. Batch correct with Harmony. 

`02_add_clusters.R`
`03_add_deconvolution.R`

`check_pxl_functions.R`
`check_pxl_in_spe_objects.R`

`sfigu_harmony.R`
Create supplemental figure showing effects of batch correction. 

### 01a_marker_genes/
`01_marker_genes.R`
Make spot plots of well-known marker genes and new marker genes. 

### 02_graph_based_clustering/
`01_graph_based_clustering.R`
Perform graph-based clustering within and across samples. 

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

`05_BayesSpace_subset_k2.R`

### 04_semi_supervised_clustering
`01_semi_supervised_clustering.R`
Run semi-supervised clustering with genes from pilot dataset.

### 05_ARI
`01_ARI.R`
Do ARI calculations comparing BayesSpace k=7 clustering of pilot data to manual annotations of pilot data

`02_plot_ARI.R`

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

### 07_layer_differential_expression
`01_create_pseudobulk_data.R`
`02_explore_expr_variability.R`
`03_model_BayesSpace.R`

`04_select_layer_DE_plots.R` Create violin plots for select examples of layer DEGs

`custom_plotExpression.R` Function for creating violin plots over pseudo bulked spe data

### 08_spatial_registration
`01_layer_correlation_annotation.R`
Perform spatial registration on BayesSpace clusters for k = 2 to k = 28 against manual annotations.

`02_t_stat_cor_plot.R`
Make t-stat correlation plot for figure 1. 

`03_jaccard.R`
Computer and plot Jaccard Index aka correspondence between clusters at different ks 

`libd_intermediate_layer_colors.R`
Create color pallet to expand on  `spatialLIBD::libd_layer_colors` add an intermediate color for hybrid layers like L2/3

### 09_position_differential_expression
`01_model_position.R` 

### 10_Clinical_Gene_Set_Enrichment
`01_gene_sets_HumanPilot.R`
Extract and format DE gene lists for bulk studies 

`02_enrichment_HumanPilot_sets.R`
Perform clinical geneset enrichment on several published datasets and k = 9 pseudobulked data. 

`03_gene_sets_singleNuc.R`
Extract and format DE gene lists for snRNA-seq studies 

`04_enrichment_singleNuc.R`
Run geneset enrichemnt on snRNA-seq genes lists

`gene_set_enrichment_plot_complex.R`
Adapt `spatialLIBD::gene_set_enrichment_plot.R` to use `ComplexHeatmap` tools to add barplot annotations for number of genes

### 11_gene_ontology
deprecated

### 12_spatial_registration_sc
`01_spatial_registration_sc.R`
Compute registration stats for HC clusters of snRNA-seq data

`02_cellType_correlation_annotation.R`
 Correlate enrichment stats against manual annotations and visium k = 9, k = 16 data using top 100 genes.
 Annotate by layer and add layer-level annotation for snRNA-seq data. Create heatmap for layer, k9 & k16 correlation values & annotations.

`04_spatial_registration_sn_velm.R` 
Compute registration stats for cell type populations in Velmeshev et al. snRNA-seq data set

`05_velm_correlation_annotation.R`
 Correlate Velmeshev et al. enrichment stats against manual, k=9 ,and k=16. Annotate and create corelation heatmap.
 
### 13_nnSVG
Exploratory part of this analysis, may be refined later

### 14_spatial_registration_PEC
`01_pseudobulk_data.R` Data was psuedo-bulked by `subclass` and `IndividualID` and saved so it could be used
for both the all data and Dx (by primary diagnosis, i.e. seperate for case and control samples). Done in parallel by data set.

`02_compute_registration_stats.R` Then internal steps of `spatialLIBD::registration_wrapper` (`registration_block_cor`, `registration_mode`, and `registration_stats_enrichment`) was then run on the pseudo bulk data.  Done in parallel by data set.

`03_PEC_correlation_annotation.R` The enrichment statistics are used for spatial registration with layer, k09, and k16 data. Correlation values are used for layer annotation. Plot Summary dot plot of annotation across all 8 data sets.

`04_PEC_check_expression.R` brief check of expression of some genes of interest across the data sets.

`05_compute_registration_Dx.R` Compute registration sets, but separate the data by Dx.

`06_PEC_correlation_annotation_Dx.R` Spatial registration and annotation of the Dx enrichment statistics

`07_PEC_annotation_Dx_plots.R` Create dotplot of Dx annotations

### 15_cell_composition
Some exploratory down-stream analyses using spot deconvolution results to understand the cell composition in ant/mid/post controlling for spatial domains

`01-EDA.R`, `011-EDA_Visium.R`, `012-EDA_SpotDeconv.R` explores the spot deconvoluted data, including section specific analysis, layer specific analysis of the spot deconvolued data.

### 16_position_differential_expression_noWM
Repeat `09_position_differential_expression/` process but exclude all WM associated bayesSpace clusers

`01_create_pseudobulk_data.R`
`02_model_position.R`

