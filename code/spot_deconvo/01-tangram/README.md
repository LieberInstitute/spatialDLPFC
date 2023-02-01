Scripts to run Tangram, one of the 3 deconvolution methods benchmarked and run for spot deconvolution.

- `*-prepare_anndatas_*`: subset the spatial `AnnData` object to one sample, filter to only use the marker genes as training genes, and otherwise prepare for running `Tangram`.
- `*-align_tangram_*`: run `Tangram` in "deconvolution mode", outputting CSVs with each cell types as columns, spots as rows, and counts as values.
- `08-sample_table.*`: generates a table mapping between two sets of sample-ID conventions (helper script for all spot deconvolution methods)
- `09-clean_raw_results.R`: given the initial spot-deconvolution results, create a new copy of CSV results in a format that can be trivially used with `spatialLIBD::cluster_import` to add this info to the corresponding `SpatialExperiment` object.
- `custom_tg_code.py`: define a modified version of `tangram.count_cell_annotations` that doesn't require cell centroids, a trick to produce cell-type counts without needing to perform segmentation with `squidpy`, which was prohibitively slow on our data.
