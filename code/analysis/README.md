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

### 08_layer_differential_expression
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
Perform clinical geneset enrichment on several published datasets and my k = 9 pseudobulked data. 

### 11_gene_ontology
`gene_ontology.R`
Performed gene ontology analysis for k = 9 pseudobulked data.

`gene_ontology.sh`

### 12_spatial_registration_sc
`01_spatial_registration_sc.R`
Perform spatial registartion of snRNA-seq data against manual annotations and visium k = 9, k = 16 data, and k = 28. Did this using top 100 gene as well as all genes. Also did this for just the excitatory clusters of the snRNA-seq data. 

### 13_nnSVG
`01_nnSVG.R`
Ran nnSVG just regularly across the whole tissue. Results `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/13_nnSVG/whole_tissue`

`01_nnSVG.sh`

`02_nnSVG_pairwise_subset.R`
Run nnSVG by subsetting the data for just the two clusters of interest. Results objects contain "subset" in the name. Did this for 6 pairs of clusters that Kristen was interested in. Results are here `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/13_nnSVG/pairwise` I submitted two jobs to each the shared queue, bluejay, and caracol. Bash scripts below. 

`02_nnSVG_pairwise_subset_caracol.sh`

`02_nnSVG_pairwise_subset_bluejay.sh`

`02_nnSVG_pairwise_subset_shared.sh`

`03_nnSVG_pairwise_combine_clusters_2.R`
The second way we ran nnSVG was to merge/combine the clusters that we're trying to compare and run nnSVG with cluster as a covariate. We were able to do this all at once for pairs that didn't have any over-lapping clusters. The scripts for this contain the word "combine_clusters", the results objects contain "merge_pairs." Scripts labeled with _2 was before making certain changed to the script that Lukas suggested. Results are here: `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/13_nnSVG/pairwise` This script combines the following pairs of clusters: 5 and 9, 4 and 16, 7 and 13. 

`03_nnSVG_pairwise_combine_clusters_2.sh`

`04_nnSVG_pairwise_combine_clusters_3.R`
_3 contains modifications that Lukas suggested which were to add the filter_genes_counts and filter_genes_pcspots parameters to the filter_genes() function. Results are here: `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/13_nnSVG/pairwise` This script is for pairs: 5 and 9, 4 and 16, 7 and 13. 

`05_nnSVG_pairwise_combine_clusters_12_13.R`
Run nnSVG with modifications from Lukas after combining clusters 12 and 13.

`05_nnSVG_pairwise_combine_clusters_12_13.sh`

`06_nnSVG_pairwise_combine_clusters_12_7.R`
Run nnSVG with modifications from Lukas after combining clusters 12 and 7.

`06_nnSVG_pairwise_combine_clusters_12_7.sh`

`07_nnSVG_pairwise_combine_clusters_12_16.R`
Run nnSVG with modifications from Lukas after combining clusters 12 and 6.

`07_nnSVG_pairwise_combine_clusters_12_16.sh`

`08_nnSVG_pairwise_subset_results.R`
Make plots of top 20 svgs. Only did this for the clusters pairs: 13 and 7, 16 and 4 because those were the only ones that ran all the way through. The jobs got stuck and didn't finish. Plots here: `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/13_nnSVG/pairwise/subset_13_7` and `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/13_nnSVG/pairwise/subset_16_4`

`09_nnSVG_pairwise_combine_clusters_results.R`
Make plots for top 20 svg. Plots here: `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/13_nnSVG/pairwise/combine_2`

### 14_spatial_registration_PEC
`01_pseudobulk_DevBrain.R`
Spatial registration of DevBrain dataset against manual annotations, my k = 9, and k = 16 data. Plots are here: `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/14_spatial_registration_PEC/DevBrain`

`02_pseudobulk_SZBD.R`
This study contained 2 annotated datasets with different numbers of cells. I only performed spatial registration on the first one. This dataset took a really long time to pseudobulk. 

`03_pseudobulk_CMC.R`

`04_pseudobulk_IsoHUB.R`

`05_pseudobulk_UCLA.R`




