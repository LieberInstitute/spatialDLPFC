
library("here")
library("sessioninfo")
library("SingleCellExperiment")


## output directory
# dir_rdata <- here::here("processed-data", "11_grey_matter_only", opt$spetype)
# dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
# stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
# dir_plots <- here::here("plots", "11_grey_matter_only", opt$spetype)
# dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
# stopifnot(file.exists(dir_plots))

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## Load the modeling results
load(file = here::here("processed-data", "rdata", "spe", "09_region_differential_expression", paste0("cluster_modeling_results_k", k, ".Rdata")), verbose = TRUE)

## load spe data
spe_pseudo <-
    readRDS(
        file = here::here(
            "processed-data",
            "rdata",
            "spe",
            "pseudo_bulked_spe",
            paste0("spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k", k, ".RDS")
        )
    )


## From
## https://github.com/LieberInstitute/HumanPilot/blob/879f11c7c57efd72334332b40feb3ad623e067c8/Analysis/Layer_Guesses/misc_numbers.R#L233-L246
## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) <- rownames(eb_contrasts)
fdrs0_contrasts <- apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) <- rownames(eb_contrasts)
summary(fdrs0_contrasts < 0.05)


## From
## https://github.com/LieberInstitute/HumanPilot/blob/879f11c7c57efd72334332b40feb3ad623e067c8/Analysis/Layer_Guesses/misc_numbers.R#L302-L317
f_merge <- function(p, fdr, t) {
    colnames(p) <- paste0("p_value_", colnames(p))
    colnames(fdr) <- paste0("fdr_", colnames(fdr))
    colnames(t) <- paste0("t_stat_", colnames(t))
    res <- as.data.frame(cbind(t, p, fdr))
    res$ensembl <- rownames(res)
    ## Check it's all in order
    stopifnot(identical(rownames(res), rownames(spe_pseudo)))
    res$gene <- rowData(spe_pseudo)$gene_name
    rownames(res) <- NULL
    return(res)
}

results_specificity <-
    f_merge(p = pvals0_contrasts, fdr = fdrs0_contrasts, t = t0_contrasts)
options(width = 400)
head(results_specificity)

## From
## https://github.com/LieberInstitute/HumanPilot/blob/879f11c7c57efd72334332b40feb3ad623e067c8/Analysis/Layer_Guesses/misc_numbers.R#L259-L264
pvals_contrasts <- eb_contrasts$p.value
fdrs_contrasts <- apply(pvals_contrasts, 2, p.adjust, "fdr")
dim(pvals_contrasts)
# [1] 27853    21
summary(fdrs_contrasts < 0.05)

results_pairwise <-
    f_merge(p = pvals_contrasts, fdr = fdrs_contrasts, t = eb_contrasts$t)
colnames(results_pairwise)
sort(colSums(fdrs_contrasts < 0.05))

## From
## https://github.com/LieberInstitute/HumanPilot/blob/879f11c7c57efd72334332b40feb3ad623e067c8/Analysis/Layer_Guesses/misc_numbers.R#L205-L215
f_sig <- function(type, cut = 0.05) {
    cbind(
        "n" = addmargins(table(f_stats[[type]] < cut)),
        "ratio" = addmargins(table(f_stats[[type]] < cut)) / nrow(f_stats)
    )
}
f_sig("noWM_fdr")

## From
## https://github.com/LieberInstitute/HumanPilot/blob/879f11c7c57efd72334332b40feb3ad623e067c8/Analysis/Layer_Guesses/misc_numbers.R#L383-L396
## Match the colnames to the new style
f_rename <- function(x, old, new = old) {
    old_patt <- paste0("_", old, "$")
    i <- grep(old_patt, colnames(x))
    tmp <- gsub(old_patt, "", colnames(x)[i])
    tmp <- paste0(new, "_", tmp)
    colnames(x)[i] <- tmp
    return(x)
}
results_anova <-
    f_rename(f_rename(f_rename(
        f_rename(f_stats, "f", "f_stat"), "p_value"
    ), "fdr"), "Amean")
head(results_anova)

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_specificity,
    "pairwise" = results_pairwise
)

save(
    modeling_results,
    file = here::here(
        "processed-data",
        "rdata",
        "spe",
        "09_region_differential_expression",
        paste0("parsed_modeling_results_k", k, ".Rdata")
    )
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
