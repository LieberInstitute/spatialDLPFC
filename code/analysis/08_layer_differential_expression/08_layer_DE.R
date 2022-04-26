# library(sgejobs)
# sgejobs::job_loop(
#    loops = list(spetype = c(
#        "wholegenome", "targeted"
#     )),
#     name = "03_model_pathology",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "15G")
# To execute the script builder, use: sh 03_model_pathology.sh

# Required libraries
#library("getopt")

## Specify parameters
# spec <- matrix(c(
#   "spetype", "s", 2, "character", "SPE spetype: wholegenome or targeted",
#   "help", "h", 0, "logical", "Display help"
# ), byrow = TRUE, ncol = 5)
# opt <- getopt(spec = spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
# if (!is.null(opt$help)) {
#   cat(getopt(spec, usage = TRUE))
#   q(status = 1)
# }

## For testing
# if (FALSE) {
#   opt <- list(spetype = "wholegenome")
# }


library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("rafalib")
library("limma")


## output directory
# dir_rdata <- here::here("processed-data","rdata","spe", "08_layer_differential_expression", opt$spetype)
# dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
# stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
# dir_plots <- here::here("plots", "08_layer_differential_expression", opt$spetype)
# dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
# stopifnot(file.exists(dir_plots))
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## load spe data
load(file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("sce_pseudobulk_bayesSpace_normalized_filtered_k",k,".Rdata")))

#boxplots of spots per cluster
pdf(file = here::here("plots","08_layer_differential_expression",paste0("ncells_per_cluster_k",k,".pdf")))
boxplot(ncells ~ colData(spe_pseudo)[[paste0("bayesSpace_harmony_",k)]], data = colData(spe_pseudo))
dev.off()

summary(spe_pseudo$ncells)

#drop samples with too few cells 
dim(spe_pseudo)
spe_pseudo <- spe_pseudo[,spe_pseudo$ncells >= 10]
dim(spe_pseudo)


## Extract the data
mat <- assays(spe_pseudo)$logcounts

#make mat_formula
var_oi = paste0("bayesSpace_harmony_",k)
covars = c("region","age","sex")
mat_formula <- as.formula(paste("~","0","+",var_oi,"+",paste(covars, collapse=" + ")))

## Compute correlation
## Adapted from https://github.com/LieberInstitute/Visium_IF_AD/blob/7973fcebb7c4b17cc3e23be2c31ac324d1cc099b/code/10_spatial_registration/01_spatial_registration.R#L134-L150
mod<- model.matrix(mat_formula,
                   data = colData(spe_pseudo))



message(Sys.time(), " running duplicateCorrelation()")
corfit <- duplicateCorrelation(mat, mod,
                               block = spe_pseudo$sample_id
)
message("Detected correlation: ", corfit$consensus.correlation)

######### ENRICHMENT t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity.R#L1423-L1443
#cluster_idx <- splitit(spe_pseudo$bayesSpace_harmony_9) #split by layers not path_grups
cluster_idx <- splitit(colData(spe_pseudo)[,var_oi]) 

message(Sys.time(), " running the enrichment model")
eb0_list <- lapply(cluster_idx, function(x) {
  res <- rep(0, ncol(spe_pseudo))
  res[x] <- 1
  mres_formula <- as.formula(paste("~","res","+",paste(covars, collapse=" + ")))
  m <- with(colData(spe_pseudo),
            model.matrix(res_formula)) 
  eBayes(
    lmFit(
      mat,
      design = m,
      block = spe_pseudo$sample_id,
      correlation = corfit$consensus.correlation
    )
  )
})

######### PAIRWISE t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity.R#L1355-L1383

## Build a group model
#mod <- with(colData(spe_pseudo), model.matrix(~ 0 + bayesSpace_harmony_9 + region + age + sex))
#colnames(mod) <- gsub("bayesSpace_harmony_9", "", colnames(mod))
#colnames(mod) <- gsub("\\+", "pos", colnames(mod)) #can remove this

message(Sys.time(), " runnign the baseline pairwise model")
fit <-
  lmFit(
    mat,
    design = mod,
    block = spe_pseudo$sample_id,
    correlation = corfit$consensus.correlation
  )
eb <- eBayes(fit)


## Define the contrasts for each pathology group vs another one
message(Sys.time(), " run pairwise models")
cluster_combs <- combn(colnames(mod), 2)
cluster_constrats <- apply(cluster_combs, 2, function(x) {
  z <- paste(x, collapse = "-")
  makeContrasts(contrasts = z, levels = mod)
})
rownames(cluster_constrats) <- colnames(mod)
colnames(cluster_constrats) <-
  apply(cluster_combs, 2, paste, collapse = "-")
eb_contrasts <- eBayes(contrasts.fit(fit, cluster_constrats))


######### ANOVA t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity_fstats.R#L24-L85

## From layer_specificity.R
fit_f_model <- function(sce) {  #will want to do this with and without white matter, look at original code from link above
  message(paste(Sys.time(), "starting the model run"))
  
  ## Extract the data
  mat <- assays(sce)$logcounts
  
  ## For dropping un-used levels
  #sce$bayesSpace_harmony_9 <- factor(sce$bayesSpace_harmony_9)
  colData(sce)[[var_oi]] <- as.factor(colData(sce)[[var_oi]])
  
  ## Build a group model
  #mod <- with(colData(sce), model.matrix(~ 0 + bayesSpace_harmony_9 + region + age + sex)) #remember to adjust for age or sex 
  #colnames(mod) <- gsub("bayesSpace_harmony_9", "", colnames(mod))
  #colnames(mod) <- gsub("\\+", "pos", colnames(mod))
  
  ## Takes like 2 min to run
  corfit <-
    duplicateCorrelation(mat, mod, block = sce$subject)
  message(paste(Sys.time(), "correlation:", corfit$consensus.correlation))
  fit <-
    lmFit(
      mat,
      design = mod,
      block = sce$subject,
      correlation = corfit$consensus.correlation
    )
  eb <- eBayes(fit)
  return(eb)
}

ebF_list <-
  lapply(list("noWM" = spe_pseudo), fit_f_model)

## Extract F-statistics
f_stats <- do.call(cbind, lapply(names(ebF_list), function(i) {
  x <- ebF_list[[i]]
  top <-
    topTable(
      x,
      coef = 2:ncol(x$coefficients), # CAREFUL make sure you pick columns from mod that are for your coefiicients of interest
      sort.by = "none",
      number = length(x$F)
    )
  # identical(p.adjust(top$P.Value, 'fdr'), top$adj.P.Val)
  res <- data.frame(
    "f" = top$F,
    "p_value" = top$P.Value,
    "fdr" = top$adj.P.Val,
    "AveExpr" = top$AveExpr,
    stringsAsFactors = FALSE
  )
  colnames(res) <- paste0(i, "_", colnames(res))
  return(res)
}))
f_stats$ensembl <- rownames(spe_pseudo)
f_stats$gene <- rowData(spe_pseudo)$gene_name
rownames(f_stats) <- NULL

head(f_stats)


save(
  f_stats,
  eb0_list,
  eb_contrasts,
  file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("cluster_modeling_results_k",k,".Rdata"))
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()