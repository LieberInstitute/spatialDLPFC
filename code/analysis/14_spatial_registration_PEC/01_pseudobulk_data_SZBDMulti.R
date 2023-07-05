library("zellkonverter")
library("SingleCellExperiment")
library("spatialLIBD")
library("jaffelab")
library("here")
library("sessioninfo")

#### load dataset  ####
subset_filepaths <- list.files(here('raw-data', 'psychENCODE', 'version6', 'SZBDMulti-Seq_subset'), full.names = TRUE)

sce_pseudo_list <- purrr::map(subset_filepaths, function(filepath){
  stopifnot(file.exists(filepath))
  
  message(Sys.time(), " - Reading data from: ", filepath)
  sce <- readH5AD(file = filepath)
  sce <- readH5AD(file = subset_filepaths[[1]])
  
  message("\nSCE Dimesions:")
  dim(sce)
  
  print(colnames(colData(sce)))
  message("Cell Types:")
  ## must be syntactically valid
  colData(sce)$cellType <- as.factor(make.names(colData(sce)$subclass))
  table(sce$cellType)
  
  # identify annotation/cluster labels
  rowData(sce)$gene_name <- rownames(sce) # save gene name as column of rowData
  rownames(sce) <- rowData(sce)$featureid # have to make row names of object the ensembl id instead of gene names
  
  ## Logcounts
  # default “X” contain the log-normalized counts
  message(Sys.time(), " revert to counts")
  
  ## check for all 0s (just first 100 cols for mem)
  stopifnot(any(assays(sce)$X[, 1:100] != 0))
  
  counts(sce) <- assays(sce)$X # change to normalized counts
  # counts(sce)[counts(sce) != 0] <- (2^counts(sce)[counts(sce) != 0])-1 # Replace just non-zero values
  counts(sce)@x <- 2^(counts(sce)@x) - 1 ## remove log2(counts + 1)
  
  #### Pseudobulk ####
  message(Sys.time(), " Pseudobulk")
  sce_pseudo <- registration_pseudobulk(sce,
                                        var_registration = "cellType",
                                        # var_sample_id = "sampleID",
                                        var_sample_id = "individualID", ## for SZDBMulti-Seq
                                        covars = NULL
  )
  
  message("\nSCE Pseudobulk Dimesions:")
  dim(sce_pseudo)
  return(sce_pseudo)
})

## Save results
saveRDS(sce_pseudo_list,
    file = here(
        "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
        paste0("pseudobulk_list_SZBDMulti.rds")
    )
)

# sgejobs::job_single('01_pseudobulk_data_SZBDMulti', create_shell = TRUE, memory = '200G', command = "Rscript 01_pseudobulk_data_SZBDMulti.R")


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

