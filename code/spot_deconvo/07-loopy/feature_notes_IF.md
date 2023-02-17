This is a [Samui](https://github.com/chaichontat/samui) app, which displays data the Visium-SPG samples from the [LIBS spatial DLPFC project](http://research.libd.org/spatialDLPFC/). Here, you can interactively visualize individual spots, along with many associated features, such as spot deconvolution results and gene expression.

# Features

## Deconvolution

Spot deconvolution was performed with two software tools (`Tangram` and `Cell2location`) and at two cell-type resolutions (`broad` and `layer`). The total number of cells of a given cell type can be displayed by selecting the feature with this naming format:

`[cell-type resolution]_[software tool]_[cell type]`

For example, the feature `broad_cell2location_excit` shows the number of excitatory (a broad cell type) cells predicted by `Cell2location` in each spot.

## CellsFiltered

These features show information about individual cells segmented by `Cellpose`. The `gfap`. `neun`, `olig2`, and `tmem119` features show the mean fluorescence of the respective immunofluorescent image channels within the segmented cell, using 8-bit pixel intensities. `area` and `dist` show the pixel area and pixel distance to the nearest spot's center, respectively, for the cell.

## Kmeans

Gene expression for each spot was using to perform k-means clustering at different values of k.

## Genes

View the log counts of all measured genes in the experiment at the spot level.
