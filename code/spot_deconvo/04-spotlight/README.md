Scripts to run `SPOTlight`, one of the 3 deconvolution methods benchmarked and run for spot deconvolution.

- `01-IF.*`: run `SPOTlight` on the 4 IF samples, outputting CSVs with each cell types as columns, spots as rows, and counts as values. Note here that to obtain counts, we multiply `Cellpose`-determined total per-spot counts by `SPOTlight`-determined cell-type proportions.
- `02-nonIF.*`: run `SPOTlight` on the 30 nonIF samples, similarly as the IF data.
- `03-move_results.*`: originally, `SPOTlight` was run in a way that was not reproducible, and the `SingleCellExperiment` object was subset to only 100 cells per cell type, which was suspected to lead to poor results. This script keeps a copy of the old results (for comparison) and moves the newer reproducible results to the original location for use downstream.
- `04-compare_settings.*`: using the `01-IF.*` and `02-nonIF.*` scripts, `SPOTlight` was run, attempting different choices for the number of cells per type to subset the `SingleCellExperiment` object: 100, 1000, and no subsetting. This script compares deconvolution results between these 3 subsetting choices (via plots).
