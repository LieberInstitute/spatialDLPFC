
library("SpatialExperiment")
library("spatialLIBD")
library("nnSVG")
library("here")
library("sessioninfo")

cluster_k = "bayesSpace_harmony_16"

## get sample and domains from args
args = commandArgs(trailingOnly=TRUE)
sample = args[1]
domain_x = args[2]
domain_y = args[3]

## Example
# sample = "Br2720_ant"
# domain_x = 5
# domain_y = 9

message("Running sample: ", sample, ", Domains", cluster_k, ": ", domain_x, " vs. ", domain_y)

message(Sys.time(), " - Loading spe")
load(here("processed-data", "rdata","spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

stopifnot(sample %in% spe$sample_id)
stopifnot(domain_x %in% spe[[cluster_k]])
stopifnot(domain_y %in% spe[[cluster_k]])

message(Sys.time(), " - Subset spe")

spe <- spe[, spe$sample_id == sample & spe[[cluster_k]] %in% c(domain_x, domain_y)]
message("ncol: ", ncol(spe))

## stop if < 65 spots 
stopifnot(ncol(spe) > 65)

## filter out mito gnees
spe <- filter_genes(spe)

message("nrow: ", nrow(spe))
# re-calculate logcounts after filtering (?) - check this

message(Sys.time(), " - Run nnSVG")

## create Model
mod <- model.matrix(~ spe$bayesSpace_harmony_9)

# set seed for reproducibility
set.seed(1019)
# Run nnSVG
spe <- nnSVG(spe, X = mod, n_threads = 1)

message(Sys.time(), " - Done! save data")
nnSVG_data <- as.data.frame(rowData(spe))
nnSVG_data$Sample <- sample
nnSVG_data$combo <- paste0(order(c(domain_x, domain_y)), collapse = ",")

# directory to save bayesspace informed results
data_dir <- here("processed-data", "rdata", "spe", "13_nnSVG", "01_nnSVG_pairwise")

# save bayesspace nnSVG results
fn_out <- paste0("nnSVG_k16-", sample,"-",domain_x, "v", domain_y , ".RData")

message("Saving to: ", fn_out)
save(nnSVG_data, file = here(data_dir, fn_out))

## Use sgejobs to write sh files 
# sgejobs::job_single('01_nnSVG_pairwise_single_test', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 01_nnSVG_pairwise.R Br2720_ant 5 9")

# combos <- readLines("nnSVG_input.txt")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
