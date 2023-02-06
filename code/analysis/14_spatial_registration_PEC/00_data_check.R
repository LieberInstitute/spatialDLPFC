library("SingleCellExperiment")
library("tidyverse")
library("here")
library("jaffelab")
library("sgejobs")

data_dir <-
    here(
        "processed-data",
        "rdata",
        "spe",
        "14_spatial_registration_PEC"
    )


#### Define PEC & store cell type details ####
pec_cell_types <- c(
    # Non-neuronal cells (8)
    "Astro",
    "Endo",
    "Immune",
    "Micro",
    "OPC",
    "Oligo",
    "PC",
    "SMC",
    # Excit (9)
    "L2/3 IT",
    "L4 IT",
    "L5 ET",
    "L5 IT",
    "L5/6 NP",
    "L6 CT",
    "L6 IT",
    "L6 IT Car3",
    "L6b",
    # Inhib (10)
    "Chandelier",
    "Lamp5",
    "Lamp5 Lhx6",
    "Pax6",
    "Pvalb",
    "Sncg",
    "Sst",
    "Sst Chodl",
    "VLMC",
    "Vip"
)

pec_cell_type_tb <- tibble(
    cell_type = factor(pec_cell_types, levels = pec_cell_types),
    cluster = make.names(cell_type),
    ct_cat = unlist(map2(c("Non-neuronal", "Excit", "Inhib"), c(8, 9, 10), rep))
)

save(pec_cell_type_tb, file = here(data_dir, "pec_cell_type_tb.Rdata"))

## color set
pec_dataset_colors <- c(
    CMC = "#00d995",
    `DevBrain-snRNAseq` = "#362600",
    IsoHuB = "#bc17d7",
    LIBD = "#ff6266",
    `MultiomeBrain-DLPFC` = "#02e3fc",
    PTSDBrainomics = "#b3005a",
    `SZBDMulti-Seq` = "#beae00",
    `UCLA-ASD` = "#afacff"
)

save(pec_dataset_colors, file = here(data_dir, "pec_dataset_colors.Rdata"))

#### Chck V3 data ####
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


# job_loop(
#     loops = list(PE_data = h5ad_files),
#     name = "01_pseudobulk_data",
#     create_shell = TRUE
# )


#### V4 files ####
raw_data_dir_v4 <- here("raw-data", "psychENCODE", "version4")
h5ad_files_v4 <- list.files(raw_data_dir_v4, pattern = "h5ad")
names(h5ad_files_v4) <- ss(basename(h5ad_files_v4), "_")

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

#### V5 files ####
raw_data_dir_v5 <- here("raw-data", "psychENCODE", "version5")
h5ad_files_v5 <- list.files(raw_data_dir_v5, pattern = "h5ad")
(names(h5ad_files_v5) <- ss(basename(h5ad_files_v5), "_"))

# CMC                            DevBrain-snRNAseq
# "CMC_annotated.h5ad"           "DevBrain-snRNAseq_annotated.h5ad"
# IsoHuB                                         LIBD
# "IsoHuB_annotated.h5ad"                        "LIBD_annotated.h5ad"
# MultiomeBrain-DLPFC                               PTSDBrainomics
# "MultiomeBrain-DLPFC_annotated.h5ad"              "PTSDBrainomics_annotated.h5ad"
# SZBDMulti-Seq                                     UCLA-ASD
# "SZBDMulti-Seq_annotated.h5ad" "UCLA-ASD_annotated_mismatches_removed.h5ad"
# UCLA-ASD
# "UCLA-ASD_annotated.h5ad"
# UCLA-ASD ## use mismatches_removed?
# "UCLA-ASD_annotated_mismatches_removed.h5ad"

datasets_v5 <- unique(names(h5ad_files_v5))

## add LIBD & PTSD
# job_loop(
#     loops = list(PE_data = c("LIBD", "PTSDBrainomics")),
#     name = "01_pseudobulk_data_check",
#     create_shell = TRUE
# )


#### Check Dx in metadata ####

metadata_fn <- map(datasets_v5, ~ here("raw-data", "psychENCODE", "metadata", paste0("individual_", .x, "_metadata.csv")))
names(metadata_fn) <- datasets_v5


## no metadata for multiomeBrain
metadata <- map(metadata_fn[map_lgl(metadata_fn, file.exists)], ~ read.csv(.x))

map(metadata, ~ .x |> dplyr::count(primaryDiagnosis))

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
# $LIBD
# primaryDiagnosis  n
# 1          control 10
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


dx_data <- map(metadata, ~ .x |>
    dplyr::select(individualID, primaryDiagnosis) |>
    dplyr::mutate(primaryDiagnosis = gsub("Not Applicable|control", "Control", primaryDiagnosis)))

## for mulitome (for now)

sce_pseudo <- readRDS(file = here(
    "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
    paste0("pseudobulk_MultiomeBrain-DLPFC.rds")
))

dx_data$`MultiomeBrain-DLPFC` <- unique(data.frame(
    individualID = sce_pseudo$individualID,
    primaryDiagnosis = ss(as.character(sce_pseudo$sampleID), "-", 2)
))

map(dx_data, ~ .x |> dplyr::count(primaryDiagnosis))

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

walk2(dx_data["PTSDBrainomics"], names(dx_data["PTSDBrainomics"]), ~ write.csv(.x,
    file = here(
        "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
        paste0("primaryDiagnosis_", .y, ".csv")
    ),
    row.names = FALSE
))

## What outputs exist?
output_dir <- here("processed-data", "rdata", "spe", "14_spatial_registration_PEC")
list.files(output_dir, pattern = "pseudobulk")
# [1] "pseudobulk_CMC.rds"                 "pseudobulk_DevBrain-snRNAseq.rds"   "pseudobulk_IsoHuB.rds"
# [4] "pseudobulk_LIBD.rds"                "pseudobulk_MultiomeBrain-DLPFC.rds" "pseudobulk_PTSDBrainomics.rds"
# [7] "pseudobulk_SZBDMulti-Seq.rds"       "pseudobulk_SZBDMulti.rds"           "pseudobulk_UCLA-ASD.rds"

list.files(output_dir, pattern = "registration_stats")
# [1] "registration_stats_CMC.rds"               "registration_stats_DevBrain-snRNAseq.rds"
# [3] "registration_stats_IsoHuB.rds"            "registration_stats_UCLA-ASD.rds"
# [5] "registration_stats_Urban-DLPFC.rds"


## No metadata for LIBD or PTSD 2/1

pb_sce <- readRDS(here(output_dir, "pseudobulk_PTSDBrainomics.rds"))
colData(pb_sce)

samp <- as.character(unique(pb_sce$sampleID))

gsub("-(\\d\\d\\d\\d)", "\\1", samp)


cat(as.character(unique(pb_sce$sampleID)), sep = "\n")
