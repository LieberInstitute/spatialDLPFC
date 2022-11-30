Several processing steps for spot deconvolution (e.g. moving R objects to python, finding marker genes, etc) are shared between methods (`tangram`, `cell2location`, `SPOTlight`) or involve comparison of methods. This directory contains code for such processing steps.

We're performing spot deconvolution at two different cell-type resolutions, and every R script in this directory has a variable called `cell_group` which can be set to either `"layer"` or `"broad"` to appropriately run each script at the right resolution.

- `01-r_to_python.*`: starting from the IF and non-IF `SpatialExperiment`s and the `SingleCellExperiment`, filter out some cells, convert R objects to `AnnData` objects for use in python, and rank all genes as markers.
- `02-find_markers.*`: find 25 "mean ratio" markers for each cell type and produce many exploratory plots to assess the quality of these markers
- `0*-deconvo_figure_*.*`: (deprecated) produce a spatial plot showing individual cells present in each spot, as determined by the 3 computational methods
- `05-assess_methods.*`: a large script producing all of the benchmarking plots and IF deconvolution results, including the figures we'll use in the paper
- `06-compare_total_counts.*`: exploratory script to see if there was concordance in total cell counts for spatially adjacent tissues
- `07-add_to_spe_IF.*`: for the Visium-IF data, add cell-type counts to the `SpatialExperiment` object and save a new copy. This includes counts from each deconvolution tool as well as the "ground-truth" counts from `cellpose` + the trained `DecisionTreeClassifier`.
- `08-gather_results.*`: compile all deconvolution results into 8 CSV files total for easy use downstream (e.g. in `05-assess-methods.*`). Creates 2 CSVs per cell group (`broad`/`layer`) and per dataset (`IF` and `nonIF`): one raw version, with the cell types at original resolution, and a collapsed version where cell types are collapsed onto those measured from the IF images/ CART (just `neuron`, `oligo`, `micro`, `astro`, and `other`).
- `09-result_plots_nonIF.*`: an analogous script to `05-assess_methods.*` but for nonIF instead of IF
- `10-combined_plots.*`: the `05-assess_methods.*` and `09-result_plots_nonIF.*` scripts are structured such that analyses can be run at *either* `broad` or `layer`-level resolution. In contrast, this script produces plots that compare *both* `broad` and `layer`-level results (for IF data only)
- `11-export_results.*`: `spatialLIBD::cluster_import` did not seem to support importing the existing results (either directly from the deconvolution tools or from the outputs of `08-gather_results.*`). This script produces `spatialLIBD::cluster_import`-compatible CSVs for each deconvo tool, cell group, and dataset.
