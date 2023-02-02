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

### 07_layer_differential_expression
`01_preliminary_analysis.R`
This is where I created the pseudobulked objects for all values of k, filter and normalize. Then run pca. Save pseudobulked object for each k. 

`01_preliminary_analysis.sh`

`02_explanatory_variables_cluster.R`
Plot explanatory variables using getVarianceExplained().

`02_explanatory_variables_cluster.sh`

`03_layer_DE.R`
Differential expression among layer performed in 3 ways: enrichment, pairwise, and anova. Done for all k's. Modeling_results saved. 

`03_layer_DE.sh`

`03a_manual_layer_DE.R`
Attempt to do layer differential expression a different way. Didn't end up using. 

`03a_manual_layer_DE.sh`

`04_layer_parse_modeling_results.R`
Parse and reformat all modeling results for futher downstream analysis and shiny app.

`04_layer_parse_modeling_results.sh`

`05_plotting_layer_DE.R`
Expression plots for top differentially expressed genes. Plots are in the trash folder inside the corresponding plots folder for this set of scripts. I guess we decided not to use them. 

`05_plotting_layer_DE.sh`

`06_top_500_enriched_genes.R`
Create csv's of top 500 enriched genes for Kristen

`07_enriched_genes_heatmap.R`
Create heatmap of top 5 enriched genes for each layer.

`08_other_enrichment_plots.R`
Plots of the numbers of enriched and downregulated genes in each layer. 

### 08_spatial_registration
`01_layer_correlation_annotation.R`
Perform spatial registration on BayesSpace clusters for k = 2 to k = 28 against manual annotations.

`02_t_stat_cor_plot.R`
Make t-stat correlation plot for figure 1. 

`03_jaccard.R`
Computer and plot Jaccard Index aka correspondence between clusters at different ks 

`libd_intermediate_layer_colors.R`
Create color pallet to expand on  `spatialLIBD::libd_layer_colors` add an intermediate color for hybrid layers like L2/3

### 09_region_differential_expression
`01_preliminary_analysis_region.R`
Create pseudobulked objects and filter genes using region as the group variable. Do this for all values of k. Run pca as well. 

`01_preliminary_analysis_region.sh`

`02_explanatory_variables_region.R`

`02_explanatory_variables_region.sh`

`03_region_DE.R`
Differential expression analysis among regions performed in 3 ways: enrichment, pairwise and anova. Done for all k's modeling results saved. 

`03_region_DE.sh`

`04_region_parse_modeling_results.R`
Parse and reformat modeling results for downstream analysis and shiny app.

`04_region_parse_modeling_results.sh`

`05_enrichment_heatmap.R`
Make heatmap of top 5 enriched genes in each region:cluster for k = 9. Save data object of top 25 enriched genes per region here: `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/09_region_differential_expression/top_25_sig_genes.rds` Save csv of enriched genes in same location. 

### 10_Clinical_Gene_Set_Enrichment
`01_Clinical_Gene_Set_Enrichment.R`
Perform clinical geneset enrichment on several published datasets and k = 9 pseudobulked data. 

### 11_gene_ontology
`gene_ontology.R`
Performed gene ontology analysis for k = 9 pseudobulked data.

`gene_ontology.sh`

### 12_spatial_registration_sc
`01_spatial_registration_sc.R`
Perform spatial registartion of snRNA-seq data against manual annotations and visium k = 9, k = 16 data, and k = 28. Did this using top 100 gene as well as all genes. Also did this for just the excitatory clusters of the snRNA-seq data. 

### 13_nnSVG
Exploratory part of this analysis, may be refined later
`00_mod_check.R`
`01_nnSVG_pairwise_modTest.R`
`01_nnSVG_pairwise.R`
`02_compile_nnSVG_output.R`
`03_summarize_nnSVG_output.R`

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





