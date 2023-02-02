
## Anlaysis

1. `01_pseudobulk_data.R` Data was psuedo-bulked by `subclass` and `IndividualID` and saved so it could be used
for both the all data and Dx (by primary diagnosis, i.e. seperate for case and control samples). Done in parallel by data set.

2. `02_compute_registration_stats.R` Then internal steps of `spatialLIBD::registration_wrapper` (`registration_block_cor`, `registration_mode`, and `registration_stats_enrichment`) was then run on the pseudo bulk data.  Done in parallel by data set.

3. `03_PEC_correlation_annotation.R` The enrichment statistics are used for spatial registration with layer, k09, and k16 data. Correlation values are used for layer annotation. Plot Summary dot plot of annotation across all 8 data sets.

4. `04_PEC_check_expression.R` brief check of expression of some genes of interest across the data sets.

5. `05_compute_registration_Dx.R` Compute registration sets, but separate the data by Dx.

6. `06_PEC_correlation_annotation_Dx.R` Spatial registration and annotation of the Dx enrichment statistics

7. `07_PEC_annotation_Dx_plots.R` Create dotplot of Dx annotations


## Dataset Info  

Data downloaded on 1/31/2023 (V5)

**What assay to use?**

From Docs:
"Given the default operations in
Pegasus, we had to explicitly create a new matrix to store the raw count values: “raw_new” is
the name of the matrix containing the raw counts for each study, while the default “X” and
“raw.X” contain the log-normalized counts."

**What annotations to use?**

"The final annotations are contained under the “subclass” column"

Subclass with `make.names()` applied = `cellType`
