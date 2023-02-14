Spatial domain-level documentation
==================================

This document describes the layer-level portion of the shiny web application made by the  [`spatialLIBD`](https://bioconductor.org/packages/spatialLIBD) Bioconductor package. You can either find the documentation about this package through [Bioconductor](https://bioconductor.org/packages/spatialLIBD) or at the [`spatialLIBD` documentation website](http://lieberinstitute.github.io/spatialLIBD). Below we explain the options common across tabs and each of the tabs at the layer-level data. As explained in the documentation, the layer-level data is the result of pseudo-bulking the spot-level data to compress it, reduce sparsity and power more analyses.

## Slides and videos

You might find the following slides useful for understanding the features from this part of the web application. Particularly slides 10-12 and 15-22.

<script async class="speakerdeck-embed" data-id="dde92cd6dfc04f9589770e074915658f" data-ratio="1.33333333333333" src="//speakerdeck.com/assets/embed.js"></script>

These slides were part of our 2021-04-27 webinar for BioTuring that you can watch from [their website](https://bioturing.com/sources/webinar/60752954a433e26dd8affcbd) or YouTube:

<iframe width="560" height="315" src="https://www.youtube.com/embed/S8884Kde-1U" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

A recording of an earlier version of this talk is also available on YouTube.

<iframe width="560" height="315" src="https://www.youtube.com/embed/aD2JU-vUv54" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

You might also be interested in this video demonstration of `spatialLIBD` for the [LIBD rstats club](http://research.libd.org/rstatsclub/). Particularly starting at minute 26 with 25 seconds.

<iframe width="560" height="315" src="https://www.youtube.com/embed/LZ2kvCiRVdM?start=1584" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Raw summary

Before the documentation, this tab displays the [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment) object that contains the spatial domain-level data (layer-level in other `spatialLIBD` apps). It's basically useful to know that the data has been loaded and that you can start navigating the app. If you wish to download this data, use the following command.

```{r}
## Check that you have a recent version of spatialLIBD installed
stopifnot(packageVersion("spatialLIBD") >= "1.11.6")

## Download sce data
sce_pseudo <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium_pseudobulk")
```

Throughout the rest of this document, we'll refer to this object by the name `sce_pseudo`.

This tab also shows the statistical modeling results, described below, that you can access locally and re-shape using the following code.

```{r}
## Reproduce locally with
modeling_results <- fetch_data("spatialDLPFC_Visium_modeling_results"")
sig_genes <-
        spatialLIBD::sig_genes_extract_all(
            n = nrow(sce_pseudo),
            modeling_results = modeling_results,
            sce_layer = sce_pseudo
        )
```

## Common options

* `Model results`: the statistical modeling results to use. We computed three different types of models:
  1. `enrichment`: one spatial domain against all the the other spatial domains. Results in t-statistics.
  2. `pairwise`: one spatial domain against another one. Results in t-statistics with two-sided p-values.
  3. `anova`: changes among the spatial domains (adjusting for the mean expression) using the data from all spatial domains (`all`). Note that great changes between the white matter and grey matter will greatly influence these results.

## Reduced dim

In this panel you can visualize the spatial domain-level data (`sce_pseudo`) across reduced dimensionality representations derived from the gene expression data from the spatial domain-level pseudo-bulked data. Select which dimensionality reduction method to use with `Reduced Dimension` (PCA, MDS or with the `scater::runPCA` function which provides more stable results and shows the percent of variable explained). Then use `Color by` to choose which variable to color data by, which can be useful to identify groups of pseudo-bulked samples. The options are:

* `BayesSpace`: the main spatial domain resolution used in this website
* `age`: age of the n = 10 donors
* `diagnosis`: all donors are neurotypical controls in this dataset
* `ncells`: number of spots that were combined when pseudo-bulking
* `position`: whether the sample is from anterior, middle or posterior
* `sample_id`: sample identifier
* `sex`: sex of the n = 10 donors
* `subject`: donor identified

```{r}
## Reproduce locally with
scater::plotReducedDim(sce_pseudo)
```

## Model boxplots

This tab allows you to make a boxplot of the `logcounts` gene expression from the spatial domain-level data (`sce_pseudo`) for a given `gene`; you can search your gene by typing either the symbol or the Ensembl gene ID. The model result information displayed in the title of the plot is based on which `model results` you selected and whether you are using the short title version or not (controlled by a checkbox). We provide two different color scales you can use: the color blind friendly `viridis` as well as a custom one we used for the `paper`. Through the `Model test` selector, you can choose which particular comparison to display. For example, `Sp09D01` for the `enrichment` model means that you would display the results of comparing Sp09D01 against the rest of the spatial domains. `Sp09D01-Sp09D02` for the `pairwise` model means that you would display the results of Sp09D01 being greater than Sp09D02, while `Sp09D02-Sp09D01` is the reverse scenario. Under `pairwise`, the spatial domains not used are display in gray.

Below the plot you can find the subset of the table of results  (`sig_genes` from earlier), sort the table by the different columns, and download it as a CSV if you want. For more details about what each of these columns mean, check the [`spatialLIBD` vignette documentation](http://LieberInstitute.github.io/spatialLIBD/articles/spatialLIBD.html#extract-significant-genes).

```{r}
## Reproduce locally with
spatialLIBD::layer_boxplot()
```

## Gene Set Enrichment

This tab allows you to upload a CSV file that has a particular format as illustrated [in this example file](https://github.com/LieberInstitute/spatialLIBD/blob/master/data-raw/asd_sfari_geneList.csv). This CSV file should contain:

* one column per gene set of interest labeled as column names on the first row,
* no row names, 
* and human Ensembl gene IDs as values in the cells. 

Once you have uploaded a CSV file following this specific format, you can then check if the genes on each of your gene sets are enriched among the statistics from `model results` (`enrichment`, etc) that have a false discovery rate (FDR) adjusted p-value less than `FDR cutoff` (0.1 by default).

Similar to the `Model boxplots` tab, you can interact with the results table or download it.

```{r}
## Reproduce locally with
spatialLIBD::gene_set_enrichment()
spatialLIBD::gene_set_enrichment_plot()
```

## Spatial registration

If you have a single nucleus or single cell RNA-sequencing (snRNA-seq)  (scRNA-seq) dataset, you might group your cells into clusters. Once you do, you could compress the data by pseudo-bulking (like we did to go from `spe` to `sce_pseudo`). You could then compute `enrichment` (`pairwise`, `anova`) statistics for your cell clusters. If you do so, you can then upload a specially formatted CSV file just like the one in [this example file](https://github.com/LieberInstitute/spatialLIBD/blob/master/data-raw/tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer.csv). This file has:

* column names,
* human Ensembl gene IDs as the row names (first column, no name for the column),
* statistics (numeric values) for the cells.

Once you have uploaded a CSV file following this specific format, you can then assess whether the correlation between your statistics and the ones from our spatial domains for the subset of genes (Ensembl ids) present in both. The resulting heatmap and interactive correlation matrix (which again you can interact with and download) can be useful if you are in the process of labeling your sn/scRNA-seq clusters or simply want to compare them against the spatial domain-specific data we have provided. This can also be used for new spatially-resolved transcriptomics datasets.

Finally, you can change the `Maximum correlation` for visualization purposes on the heatmap as it will change the dynamic range for the colors.

```{r}
## Reproduce locally with
spatialLIBD::layer_stat_cor()
spatialLIBD::layer_stat_cor_plot()
```
