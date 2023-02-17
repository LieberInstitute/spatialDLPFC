This is a [Samui](https://github.com/chaichontat/samui) app, which displays data the Visium samples from the [LIBS spatial DLPFC project](http://research.libd.org/spatialDLPFC/). Here, you can interactively visualize individual spots, along with many associated features, such as spot deconvolution results and gene expression.

# Features

## Deconvolution

Spot deconvolution was performed with two software tools (`Tangram` and `Cell2location`) and at two cell-type resolutions (`broad` and `layer`). The total number of cells of a given cell type can be displayed by selecting the feature with this naming format:

`[cell-type resolution]_[software tool]_[cell type]`

For example, the feature `broad_cell2location_excit` shows the number of excitatory (a broad cell type) cells predicted by `Cell2location` in each spot.


## Genes

View the log counts of all measured genes in the experiment at the spot level.
