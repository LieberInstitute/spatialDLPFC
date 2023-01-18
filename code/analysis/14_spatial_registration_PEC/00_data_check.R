
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

map(datasets, ~ list.files(here(raw_data_dir, .x), pattern = "-snRNAseq_annotated.h5ad"))

map(datasets, ~ list.files(here(raw_data_dir, .x), pattern = ".h5ad"))

h5ad_files <- c(
    CMC = "CMC/CMC-CellHashing_annotated.h5ad",
    `DevBrain-snRNASeq` = "DevBrain-snRNAseq/DevBrain-snRNAseq_annotated.h5ad",
    IsoHuB = "IsoHuB/IsoHuB-snRNAseq_annotated.h5ad",
    `SZBDMulti-Seq` = "SZBDMulti-Seq/SZBDMulti-Seq_annotated.h5ad",
    `UCLA-ASD` = "UCLA-ASD/UCLA-ASD-snRNAseq_annotated_mismatches_removed.h5ad", # "UCLA-ASD-snRNAseq_annotated.h5ad" which file?
    `Urban-DLPFC` = "Urban-DLPFC/Urban-DLPFC-snRNAseq_annotated.h5ad"
)

ss(h5ad_files, "/")
# CMC DevBrain-snRNAseq            IsoHuB     SZBDMulti-Seq          UCLA-ASD       Urban-DLPFC
# "CMC"        "DevBrain"          "IsoHuB"       "SZBDMulti"            "UCLA"           "Urban"

map(h5ad_files, ~ file.exists(here("raw-data", "psychENCODE", "version3", "scRNAseq_AllenCTHarmonized", .x)))

## Check input from 01_pseudobulk
map(h5ad_files, function(input_file) {
    dataset <- gsub("-(s|S).*$", "", dirname(input_file))
    message("\n#### Running: ", dataset, " ####")
    # filepath <- here("raw-data", "psychENCODE", "version2", dataset, paste0(dataset, "-snRNAseq_annotated.h5ad"))

    ## for v3 data
    filepath <- here("raw-data", "psychENCODE", "version3", "scRNAseq_AllenCTHarmonized", input_file)
    stopifnot(file.exists(filepath))
})


job_loop(
    loops = list(PE_data = h5ad_files),
    name = "01_pseudobulk_data",
    create_shell = TRUE
)

## What outputs exist?
output_dir <- here("processed-data", "rdata", "spe", "14_spatial_registration_PEC")
list.files(output_dir, pattern = "pseudobulk")
# [1] "pseudobulk_CMC.rds"               "pseudobulk_DevBrain-snRNAseq.rds" "pseudobulk_IsoHuB.rds"           
# [4] "pseudobulk_SZBDMulti.rds"         "pseudobulk_UCLA-ASD.rds"          "pseudobulk_Urban-DLPFC.rds"

list.files(output_dir, pattern = "registration_stats")
# [1] "registration_stats_CMC.rds"               "registration_stats_DevBrain-snRNAseq.rds"
# [3] "registration_stats_IsoHuB.rds"            "registration_stats_UCLA-ASD.rds"         
# [5] "registration_stats_Urban-DLPFC.rds"   

#### V4 files ####
raw_data_dir_v4 <- here("raw-data", "psychENCODE", "version4")
h5ad_files_v4 <- list.files(raw_data_dir_v4, pattern = "h5ad")
names(h5ad_files_v4) <- ss(basename(h5ad_files_v4),"_")

# UCLA-ASD ## use mismatches_removed?
# "UCLA-ASD_annotated_mismatches_removed.h5ad" 

datasets_v4 <- unique(names(h5ad_files_v4))
# [1] "CMC"                 "DevBrain-snRNAseq"   "IsoHuB"              "MultiomeBrain-DLPFC" "SZBDMulti-Seq"      
# [6] "UCLA-ASD"

# job_loop(
#   loops = list(PE_data = datasets_v4),
#   name = "01_pseudobulk_data_check",
#   create_shell = TRUE
# )


#### Check Dx in metadata ####

metadata_fn <- map(datasets_v4, ~here("raw-data", "psychENCODE", "version4", paste0("individual_", .x, "_metadata.csv")))
names(metadata_fn) <- datasets_v4

## no metadata for multiomeBrain
metadata <- map(metadata_fn[map_lgl(metadata_fn, file.exists)], ~read.csv(.x))

map(metadata , ~.x |> dplyr::count(primaryDiagnosis))

# $CMC
# primaryDiagnosis  n
# 1          control 53
# 2    Schizophrenia 48
# 
# $`DevBrain-snRNAseq`
# primaryDiagnosis  n
# 1 Autism Spectrum Disorder 11
# 2           Not Applicable  7
# 3        Williams Syndrome  3
# 
# $IsoHuB
# primaryDiagnosis n
# 1          control 3
# 2   Not Applicable 2
# 
# $`SZBDMulti-Seq`
# primaryDiagnosis  n
# 1 Bipolar Disorder 24
# 2   Not Applicable 24
# 3    Schizophrenia 24
# 
# $`UCLA-ASD`
# primaryDiagnosis  n
# 1 Autism Spectrum Disorder 62
# 2                  control 77
# 3                     <NA>  1


dx_data <- map(metadata , ~.x |> 
                 dplyr::select(individualID, primaryDiagnosis) |>
                 dplyr::mutate(primaryDiagnosis = gsub("Not Applicable|control","Control",primaryDiagnosis)))

## for mulitome (for now)

sce_pseudo <- readRDS(file = here(
  "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
  paste0("pseudobulk_MultiomeBrain-DLPFC.rds")
))

dx_data$`MultiomeBrain-DLPFC` <- unique(data.frame(individualID = sce_pseudo$individualID, 
                       primaryDiagnosis = ss(as.character(sce_pseudo$sampleID), "-", 2)))

map(dx_data , ~.x |> dplyr::count(primaryDiagnosis))

# $CMC
# primaryDiagnosis  n
# 1          Control 53
# 2    Schizophrenia 48
# 
# $`DevBrain-snRNAseq`
# primaryDiagnosis  n
# 1 Autism Spectrum Disorder 11
# 2                  Control  7
# 3        Williams Syndrome  3
# 
# $IsoHuB
# primaryDiagnosis n
# 1          Control 5
# 
# $`SZBDMulti-Seq`
# primaryDiagnosis  n
# 1 Bipolar Disorder 24
# 2          Control 24
# 3    Schizophrenia 24
# 
# $`UCLA-ASD`
# primaryDiagnosis  n
# 1 Autism Spectrum Disorder 62
# 2                  Control 77
# 
# $`MultiomeBrain-DLPFC`
# primaryDiagnosis  n
# 1          Bipolar 10
# 2          Control  5
# 3    Schizophrenia  6

walk2(dx_data, names(dx_data), ~write.csv(.x, file = here("processed-data", "rdata", "spe", "14_spatial_registration_PEC",
                                                          paste0("primaryDiagnosis_",.y,".csv")),
                                          row.names = FALSE))

