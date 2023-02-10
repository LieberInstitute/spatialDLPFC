library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(here)

# library(ggspavis)
# library(tidySingleCellExperiment)


spe_dat <- readRDS(here("processed-data/rdata/spe/01_build_spe",
                        "spe_filtered_final_with_clusters_and_deconvolution_results.rds")
)

# Coloc of EFNA5 -> EPHA5
coloc_spe <- calc_coloc(spe_dat, "EFNA5", "EPHA5", sample_id = "Br8667_mid")
vis_coloc(coloc_spe, "EFNA5", "EPHA5", sample_id = "Br8667_mid",
          save.path = here("plots/15_cell_composition", "coloc"))
vis_coloc_spd(coloc_spe, "EFNA5", "EPHA5", sample_id = "Br8667_mid",
          save.path = here("plots/15_cell_composition", "coloc"))

# Coloc of FYN -> EPHA5
calc_coloc(spe_dat, "FYN", "EPHA5", sample_id = "Br8667_mid") |>
vis_coloc( "FYN", "EPHA5", sample_id = "Br8667_mid",
          save.path = here("plots/15_cell_composition", "coloc"))
