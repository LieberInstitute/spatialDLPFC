Scripts to run `Cell2location`, one of the 3 deconvolution methods benchmarked and run for spot deconvolution.

- `*-prepare_anndata_*`: attach the high-resolution images to the `AnnData` object, name object attributes for compatibility with the `Cell2location` tutorial code, and otherwise prepare for running `Cell2location`.
- `*-registration_*`: run `Cell2location`, outputting CSVs with each cell types as columns, spots as rows, and counts as values.
