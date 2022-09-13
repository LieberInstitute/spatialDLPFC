library('zellkonverter')
library("SingleCellExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")


#### load dataset  ####
args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
message("Running - ", dataset)

sce_pseudo <- readRDS(file = here("processed-data","rdata","spe","14_spatial_registration_PEC",
                 paste0("pseudobulk_", dataset, ".rds")))

message("\nSCE Dimesions:")
dim(sce_pseudo)

min_ncells <- 10 
message(
  Sys.time(),
  " dropping ",
  sum(sce_pseudo$ncells < min_ncells),
  " pseudo-bulked samples that are below 'min_ncells'."
)
sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= min_ncells]

message("Cell Types:")
## must be syntactically valid
table(sce_pseudo$cellType)
## Drop Levels
sce_pseudo$cellType <- droplevels(sce_pseudo$cellType)

#### Run models ####
registration_mod <- registration_model(sce_pseudo)
block_cor <- registration_block_cor(sce_pseudo, registration_model = registration_mod)

gene_name = "gene_name"
gene_ensembl = "featureid"

results_enrichment <-
  registration_stats_enrichment(
    sce_pseudo,
    block_cor = block_cor,
    gene_ensembl = gene_ensembl,
    gene_name = gene_name
  )

results_pairwise <-
  registration_stats_pairwise(
    sce_pseudo,
    registration_model = registration_mod,
    gene_ensembl = gene_ensembl,
    gene_name = gene_name
  )

results_anova <-
  registration_stats_anova(
    sce_pseudo,
    block_cor = block_cor,
    gene_ensembl = gene_ensembl,
    gene_name = gene_name,
    prefix = prefix
  )

modeling_results <- list(
  "anova" = results_anova,
  "enrichment" = results_enrichment,
  "pairwise" = results_pairwise
)

## Save results
saveRDS(modeling_results,
        file = here("processed-data","rdata","spe","14_spatial_registration_PEC",
                    paste0("registration_stats_",dataset,".rds")))

# sgejobs::job_single('02_compute_registration_stats_DevBrain', create_shell = TRUE, memory = '25G', command = "Rscript 02_compute_registration_stats.R DevBrain")
# sgejobs::job_single('02_compute_registration_stats_CMC', create_shell = TRUE, memory = '25G', command = "Rscript 02_compute_registration_stats.R CMC")
# sgejobs::job_single('02_compute_registration_stats_IsoHuB', create_shell = TRUE, memory = '25G', command = "Rscript 02_compute_registration_stats.R IsoHuB")
# sgejobs::job_single('02_compute_registration_stats_UCLA', create_shell = TRUE, memory = '25G', command = "Rscript 02_compute_registration_stats.R UCLA-ASD")

# sgejobs::job_single('02_compute_registration_stats_SZBD', create_shell = TRUE, memory = '25G', command = "Rscript 02_compute_registration_stats.R SZBDMulti")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
