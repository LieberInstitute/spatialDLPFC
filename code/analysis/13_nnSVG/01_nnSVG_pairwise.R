library("SpatialExperiment")
library("spatialLIBD")
library("nnSVG")
library("here")
library("sessioninfo")

cluster_k <- "bayesSpace_harmony_16"

## get sample and domains from args
args <- commandArgs(trailingOnly = TRUE)

domains <- args[1]
domains_n <- as.integer(strsplit(domains, "v")[[1]])

domain_x <- domains_n[[1]]
domain_y <- domains_n[[2]]

## sample i from array job
sample_i <- as.integer(args[2])

message(Sys.time(), " - Loading spe")
load(here("processed-data", "rdata", "spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

samples <- unique(spe$sample_id)
sample <- samples[[sample_i]]

stopifnot(sample %in% spe$sample_id)
stopifnot(domain_x %in% spe[[cluster_k]])
stopifnot(domain_y %in% spe[[cluster_k]])

message("Running sample: ", sample, ", Domains", cluster_k, ": ", domain_x, " vs. ", domain_y)

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
message("K9 model breakdown:")
table(spe$bayesSpace_harmony_9)

mod <- model.matrix(~ spe$bayesSpace_harmony_9)

# set seed for reproducibility
# set.seed(1019)
seed <- 2822
message("seed = ", seed)
set.seed(seed)
# Run nnSVG
spe <- nnSVG(spe, X = mod, n_threads = 1)

message(Sys.time(), " - Done! save data")
nnSVG_data <- as.data.frame(rowData(spe))
nnSVG_data$Sample <- sample
nnSVG_data$domains <- domains

# directory to save bayesspace informed results
data_dir <- here("processed-data", "rdata", "spe", "13_nnSVG", "01_nnSVG_pairwise")

# save bayesspace nnSVG results
fn_out <- paste0("nnSVG_k16-", domains, "-", sample, ".RData")

message("Saving to: ", fn_out)
save(nnSVG_data, file = here(data_dir, fn_out))

#### Use sgejobs to write sh files ###
# sgejobs::job_single('01_nnSVG_pairwise_single_test', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 01_nnSVG_pairwise.R Br2720_ant 5 9")

## Loop through domains, use array to run samples
# domains <- readLines("nnSVG_domains.txt")
# sgejobs::job_loop(
#   loops = list(domains = domains),
#   task_num = 30,
#   name = '01_nnSVG_pairwise_loop',
#   create_shell = TRUE
# )

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
