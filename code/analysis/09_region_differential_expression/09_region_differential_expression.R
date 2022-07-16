library(edgeR)
library(scater)
library(scran)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# load pseudobulked object created from preliminary_analysis.R
load(file = here::here("processed-data", "rdata", "spe", "pseudo_bulked_spe", paste0("sce_pseudobulk_bayesSpace_k", k, ".Rdata")))


# http://bioconductor.org/books/3.14/OSCA.multisample/multi-sample-comparisons.html#creating-pseudo-bulk-samples


# loop through with pseudoBulkDGE
dim(spe_pseudo)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 10]
dim(spe_pseudo)

### Pair-wise analysis with psedudoBulkDGE######

# head(model.matrix(~factor(region) + factor(subject),as.data.frame(colData(spe_pseudo)))) #;ppl at rclub session on model matrix
# possibly have to loop through other coefficients other than region such as subject, age,
de.results.mid <- pseudoBulkDGE(spe_pseudo,
    label = colData(spe_pseudo)[[paste0("bayesSpace_harmony_", k)]], # tells it to do it one cluster at a time. to do it globally, don't need label.
    design = ~ factor(region) + factor(sex) + factor(age),
    method = "voom",
    coef = "factor(region)middle" # comes from topTable from limma, specifies the coefficient you want to do the t-test on, this is pairwise comparison. can see how they did it before with limma for the pilot study.
    # in order to run anova have to provide more than one coefficient
)

head(de.results.mid[[1]])
table(de.results.mid[[1]]$FDR < 0.05) # shows how many differentially expressed genes between middle and anterior in cluster one
rowData(spe_pseudo)[which(de.results.mid[[1]]$FDR < 0.05), ]
de.results[[1]][which(de.results.mid[[1]]$FDR < 0.05), ]

de.results.post <- pseudoBulkDGE(spe_pseudo,
    label = spe_pseudo$BayesSpace,
    method = "voom",
    design = ~ factor(region) + factor(sex) + factor(age),
    coef = "factor(region)posterior" # comes from topTable from limma, specifies the coefficient you want to do the t-test on
    # in order to run anova have to provide more than one coefficient Run 3 times for all comparisons (ant v mid), mid v post, post v ant.  this is pairwise
    # give it two coefficients to run it as an ANOVA, more than one coefficient is doing fstatistics
)

de.results.ant <- pseudoBulkDGE(spe_pseudo,
    label = spe_pseudo$BayesSpace,
    method = "voom",
    design = ~ factor(region) + factor(sex) + factor(age),
    coef = "factor(region)anterior" # comes from topTable from limma, specifies the coefficient you want to do the t-test on
    # in order to run anova have to provide more than one coefficient Run 3 times for all comparisons (ant v mid), mid v post, post v ant.  this is pairwise
    # give it two coefficients to run it as an ANOVA, more than one coefficient is doing fstatistics
)

save(
    de.results.mid,
    de.results.ant,
    de.results.post,
    file = here::here("processed-data", "rdata", "spe", "09_region_differential_expression", paste0("pairwise_de_results_region_k", k, ".Rdata"))
)


### ANOVA #####
