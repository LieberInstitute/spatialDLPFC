library("here")
library("sessioninfo")
library("SingleCellExperiment")
library(rafalib)
library("limma")

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## load spe data
spe_pseudo <-
    readRDS(
        file = here::here(
            "processed-data",
            "rdata",
            "spe",
            "pseudo_bulked_spe",
            paste0("spe_pseudobulk_bayesSpace_normalized_filtered_region_k", k, ".RDS")
        )
    )

# boxplots of spots per cluster
# pdf(file = here::here("plots", "08_layer_differential_expression", paste0("ncells_per_cluster_k", k, ".pdf")))
# boxplot(ncells ~ spe_pseudo$BayesSpace, data = colData(spe_pseudo))
# dev.off()


## Extract the data
mat <- assays(spe_pseudo)$logcounts

# make mat_formula
# var_oi = paste0("bayesSpace_harmony_",k)
var_oi <- "region"
covars <- c("BayesSpace", "age", "sex")
mat_formula <- eval(str2expression(paste("~", "0", "+", var_oi, "+", paste(covars, collapse = " + "))))

# make sure everything is  a factor
colData(spe_pseudo)[[var_oi]] <- as.factor(colData(spe_pseudo)[[var_oi]])
colData(spe_pseudo)$region <- as.factor(colData(spe_pseudo)$region)
colData(spe_pseudo)$age <- as.numeric(colData(spe_pseudo)$age)
colData(spe_pseudo)$sex <- as.factor(colData(spe_pseudo)$sex)
colData(spe_pseudo)$diagnosis <- as.factor(colData(spe_pseudo)$diagnosis)
colData(spe_pseudo)$subject <- as.factor(colData(spe_pseudo)$subject)

## Compute correlation
## Adapted from https://github.com/LieberInstitute/Visium_IF_AD/blob/7973fcebb7c4b17cc3e23be2c31ac324d1cc099b/code/10_spatial_registration/01_spatial_registration.R#L134-L150
mod <- model.matrix(mat_formula,
    data = colData(spe_pseudo)
)


message(Sys.time(), " running duplicateCorrelation()")
corfit <- duplicateCorrelation(mat, mod,
    block = spe_pseudo$sample_id
)
message("Detected correlation: ", corfit$consensus.correlation)

######### ENRICHMENT t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity.R#L1423-L1443
cluster_idx <- splitit(colData(spe_pseudo)[, var_oi])

message(Sys.time(), " running the enrichment model")
eb0_list <- lapply(cluster_idx, function(x) {
    res <- rep(0, ncol(spe_pseudo))
    res[x] <- 1
    res_formula <- paste("~", "res", "+", paste(covars, collapse = " + "))
    m <- with(
        colData(spe_pseudo),
        model.matrix(eval(str2expression(res_formula)))
    )
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
mod_p <- model.matrix(~ 0 + region,
    data = colData(spe_pseudo)
)
# colnames(mod) <- gsub("bayesSpace_harmony_9", "", colnames(mod))

message(Sys.time(), " running the baseline pairwise model")
fit <-
    lmFit(
        mat,
        design = mod_p,
        block = spe_pseudo$sample_id,
        correlation = corfit$consensus.correlation
    )
eb <- eBayes(fit)


## Define the contrasts for each pathology group vs another one
message(Sys.time(), " run pairwise models")
# cluster_combs <- combn(colnames(mod)[grep("BayesSpace",colnames(mod))], 2)
cluster_combs <- combn(colnames(mod_p), 2)
cluster_constrats <- apply(cluster_combs, 2, function(x) {
    z <- paste(x, collapse = "-")
    makeContrasts(contrasts = z, levels = mod_p)
})
rownames(cluster_constrats) <- colnames(mod_p)
colnames(cluster_constrats) <-
    apply(cluster_combs, 2, paste, collapse = "-")

# cluster_constrats <- cluster_constrats[grep("BayesSpace",rownames(cluster_constrats)),]
eb_contrasts <- eBayes(contrasts.fit(fit, cluster_constrats))


######### ANOVA t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity_fstats.R#L24-L85

## From layer_specificity.R
fit_f_model <- function(sce) { # will want to do this with and without white matter, look at original code from link above
    message(paste(Sys.time(), "starting the model run"))

    ## Extract the data
    mat <- assays(sce)$logcounts

    ## For dropping un-used levels
    # sce$bayesSpace_harmony_9 <- factor(sce$bayesSpace_harmony_9)
    colData(sce)[[var_oi]] <- as.factor(colData(sce)[[var_oi]])

    ## Build a group model
    # already made in beginning of script #remember to adjust for age or sex


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
    y <- ncol(x$coefficients) - 4
    top <-
        topTable(
            x,
            coef = 2:y,
            # coef = 2:ncol(x$coefficients), # CAREFUL make sure you pick columns from mod that are for your coefiicients of interest. will have 8
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
    file = here::here("processed-data", "rdata", "spe", "09_region_differential_expression", paste0("cluster_modeling_results_k", k, ".Rdata"))
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
