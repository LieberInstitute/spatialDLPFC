library(edgeR)
library(scater)
library(scran)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load pseudobulked object created from preliminary_analysis.R
load(file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("sce_pseudobulk_bayesSpace_k",k,".Rdata")))

#current = spe_pseudo
# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(spe_pseudo), samples=colData(spe_pseudo))
y

discarded <- spe_pseudo$ncells < 10
y <- y[,!discarded]
summary(discarded)

keep <- filterByExpr(y, group=spe_pseudo[[paste0("bayesSpace_harmony_",k)]])
y <- y[keep,]
summary(keep)

y <- calcNormFactors(y)
y$samples

# pdf(file = here::here("plots","08_layer_differential_expression",paste0("MDS_k",k,".pdf")))
# plotMDS(cpm(y, log=TRUE), 
#         col=ifelse(y$samples[[paste0("bayesSpace_harmony_",k)]]))
# dev.off()

#design <- model.matrix(~factor(subject) + factor(bayesSpace_harmony_9), y$samples) ## need to change this to use k
design <- model.matrix(~factor(subject) + factor(y$samples[[paste0("bayesSpace_harmony_",k)]]), y$samples) ## need to change this to use k

design

y <- estimateDisp(y, design)
summary(y$trended.dispersion)


pdf(file = here::here("plots","08_layer_differential_expression",paste0("BCV_k",k,".pdf")))
plotBCV(y)
dev.off()

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)

summary(fit$df.prior)

pdf(file = here::here("plots","08_layer_differential_expression",paste0("QLDisp_k",k,".pdf")))
plotQLDisp(fit)
dev.off()

res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))

topTags(res)

save(res, file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("res_layer_k",k,".Rdata")))

###create violin plots for differentially expressed genes
load(file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("sce_pseudobulk_bayesSpace_normalized_filtered_k",k,".Rdata")))

## Palette taken from `scater`
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D")

pdf(file = here::here("plots","08_layer_differential_expression",paste0("topTags_violin_k",k,".pdf")),height = 8, width = 8)
  print(
    plotExpression(spe_pseudo, exprs_values = "logcounts", features=rownames(topTags(res))[c(1:2,5:7,9)],
                   x=paste0("bayesSpace_harmony_",k), colour_by=paste0("bayesSpace_harmony_",k), point_alpha=1, point_size=1, ncol=2,
                   add_legend=F, theme_size=16)
    +  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
             axis.title.y = element_text(angle = 90, size = 14),
             panel.grid.major=element_line(colour="grey95", size=0.8),
             panel.grid.minor=element_line(colour="grey95", size=0.4)) 
  )

dev.off()



