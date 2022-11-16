
library("purrr")
library("here")
library("jaffelab")
library("sgejobs")

raw_data_dir <- here("raw-data", "psychENCODE", "version3", "scRNAseq_AllenCTHarmonized")
list.files(raw_data_dir)
# [1] "CMC"                           "DevBrain-snRNAseq"             "Documentation"                
# [4] "IsoHuB"                        "SYNAPSE_METADATA_MANIFEST.tsv" "SZBDMulti-Seq"                
# [7] "UCLA-ASD"                      "Urban-DLPFC"  

datasets <- c("CMC", "DevBrain-snRNAseq", "IsoHuB", "SZBDMulti-Seq", "UCLA-ASD", "Urban-DLPFC")
names(datasets) <- datasets

map(datasets, ~list.files(here(raw_data_dir, .x), pattern = "-snRNAseq_annotated.h5ad"))

map(datasets, ~list.files(here(raw_data_dir, .x), pattern = ".h5ad"))

h5ad_files <- c(CMC = "CMC-CellHashing_annotated.h5ad",
                  `DevBrain-snRNAseq` = "DevBrain-snRNAseq_annotated.h5ad",
                   IsoHuB = "IsoHuB-snRNAseq_annotated.h5ad",
                   `SZBDMulti-Seq` = "SZBDMulti-Seq_annotated.h5ad",
                  `UCLA-ASD` = "UCLA-ASD-snRNAseq_annotated_mismatches_removed.h5ad",# "UCLA-ASD-snRNAseq_annotated.h5ad" which file?
                  `Urban-DLPFC` = "Urban-DLPFC-snRNAseq_annotated.h5ad")

ss(h5ad_files, "-")
# CMC DevBrain-snRNAseq            IsoHuB     SZBDMulti-Seq          UCLA-ASD       Urban-DLPFC 
# "CMC"        "DevBrain"          "IsoHuB"       "SZBDMulti"            "UCLA"           "Urban" 

map2(h5ad_files, names(h5ad_files), ~file.exists(here("raw-data", "psychENCODE", "version3", "scRNAseq_AllenCTHarmonized", .y, .x)))

job_loop(
  loops = list(PE_data = h5ad_files),
  name = '01_pseudobulk_data',
  create_shell = TRUE
)

## What outputs exist?
output_dir <- here("raw-data", "psychENCODE", "version3", "scRNAseq_AllenCTHarmonized")
list.files(raw_data_dir)


walk(datasets, ~message(gsub("DATASET", .x,
                    "file.remove('01_pseudobulk_data_DATASET.sh')\nsgejobs::job_single('01_pseudobulk_data_DATASET', create_shell = TRUE, memory = '25G', command = 'Rscript 01_pseudobulk_data.R DATASET DATASET.h5ad')\n")))
