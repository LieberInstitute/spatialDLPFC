library(edgeR)
library(scater)
library(scran)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# load pseudobulked object created from preliminary_analysis.R
load(file = here::here("processed-data", "rdata", "spe", "pseudo_bulked_spe", paste0("sce_pseudobulk_bayesSpace_k", k, ".Rdata")))

# current = spe_pseudo
# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(spe_pseudo), samples = colData(spe_pseudo))
y

discarded <- spe_pseudo$ncells < 10
y <- y[, !discarded]
summary(discarded)

y <- calcNormFactors(y)
y$samples

keep <- filterByExpr(y, group = spe_pseudo[[paste0("bayesSpace_harmony_", k)]])
y <- y[keep, ]
summary(keep)


# pdf(file = here::here("plots","08_layer_differential_expression",paste0("MDS_k",k,".pdf")))
# plotMDS(cpm(y, log=TRUE),
#         col=ifelse(y$samples[[paste0("bayesSpace_harmony_",k)]]))
# dev.off()

# design <- model.matrix(~factor(subject) + factor(bayesSpace_harmony_9), y$samples) ## need to change this to use k
design <- model.matrix(~ factor(subject) + factor(y$samples[[paste0("bayesSpace_harmony_", k)]]), y$samples)

design

y <- estimateDisp(y, design)
summary(y$trended.dispersion)


pdf(file = here::here("plots", "08_layer_differential_expression", paste0("BCV_k", k, ".pdf")))
plotBCV(y)
dev.off()

fit <- glmQLFit(y, design, robust = TRUE)
summary(fit$var.prior)

summary(fit$df.prior)

pdf(file = here::here("plots", "08_layer_differential_expression", paste0("QLDisp_k", k, ".pdf")))
plotQLDisp(fit)
dev.off()

res <- glmQLFTest(fit, coef = ncol(design))
summary(decideTests(res))

topTags(res)

save(res, file = here::here("processed-data", "rdata", "spe", "08_layer_differential_expression", paste0("res_layer_k", k, ".Rdata")))
