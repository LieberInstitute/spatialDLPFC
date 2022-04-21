library(edgeR)
library(scater)
library(scran)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load pseudobulked object created from preliminary_analysis.R
load(file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("sce_pseudobulk_bayesSpace_k",k,".Rdata")))
load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("res_layer_k",k,".Rdata")))

##create violin plots for differentially expressed genes
#load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("res_layer_k",k,".Rdata")))


pdf(file = here::here("plots","08_layer_differential_expression",paste0("topTags_violin_k",k,".pdf")),height = 8, width = 8)
print(
  plotExpression(logNormCounts(spe_pseudo, size.factors=NULL), features=rownames(topTags(res)),
                 x=paste0("bayesSpace_harmony_",k), colour_by=paste0("bayesSpace_harmony_",k), point_alpha=1, point_size=1, ncol=2,
                 add_legend=F, theme_size=16)
  +  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
           axis.title.y = element_text(angle = 90, size = 14),
           panel.grid.major=element_line(colour="grey95", size=0.8),
           panel.grid.minor=element_line(colour="grey95", size=0.4)) 
)

dev.off()