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
  
  message("\nSCE Dimesions:")
  print(nrow(sce))
  print(ncol(sce))
  
  print(colnames(colData(sce)))
  print(unique(sce$individualID))
  
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

  sce_pseudo <- scuttle::aggregateAcrossCells(
    sce,
    DataFrame(
      registration_variable = sce[["cellType"]],
      registration_sample_id = sce[["individualID"]]
    )
  )
  colnames(sce_pseudo) <-
    paste0(
      sce_pseudo$registration_sample_id,
      "_",
      sce_pseudo$registration_variable
    )
  
  message("\nSCE Pseudobulk Dimesions:")
  dim(sce_pseudo)
  return(sce_pseudo)
})

# saveRDS(sce_pseudo_list,
#         file = here(
#           "processed-data", "rdata", "spe", "14_spatial_registration_PEC","pseudobulk_list_SZBDMulti-seq.rds"
#         )
# )

sce_pseudo <- do.call("cbind", sce_pseudo_list)     
dim(sce_pseudo)
# [1] 34361  1870

## complete the rest of registration_pseudobulk
min_ncells = 10
message(
  Sys.time(),
  " dropping ",
  sum(sce_pseudo$ncells < min_ncells),
  " pseudo-bulked samples that are below 'min_ncells'."
)

# 2023-07-17 14:56:21.081636 dropping 307 pseudo-bulked samples that are below 'min_ncells'.

sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= min_ncells]
dim(sce_pseudo)
# [1] 34361  1563

if (is.factor(sce_pseudo$registration_variable)) {
  ## Drop unused var_registration levels if we had to drop some due
  ## to min_ncells
  sce_pseudo$registration_variable <- droplevels(sce_pseudo$registration_variable)
}

## Drop lowly-expressed genes
message(Sys.time(), " drop lowly expressed genes")
keep_expr <-
  edgeR::filterByExpr(sce_pseudo, group = sce_pseudo$registration_variable)
sce_pseudo <- sce_pseudo[which(keep_expr), ]

## Compute the logcounts
message(Sys.time(), " normalize expression")
logcounts(sce_pseudo) <-
  edgeR::cpm(edgeR::calcNormFactors(sce_pseudo),
             log = TRUE,
             prior.count = 1
  )


## Save results
saveRDS(sce_pseudo,
        file = here(
          "processed-data", "rdata", "spe", "14_spatial_registration_PEC",
          paste0("pseudobulk_", "SZBDMulti-Seq", ".rds")
        )
)


# sgejobs::job_single('01_pseudobulk_data_SZBDMulti', create_shell = TRUE, memory = '200G', command = "Rscript 01_pseudobulk_data_SZBDMulti.R")


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

