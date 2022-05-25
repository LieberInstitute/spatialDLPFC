# https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/asd_snRNAseq_recast.R
# https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/ad_snRNAseq_recast.R
# http://research.libd.org/spatialLIBD/reference/layer_stat_cor.html
# http://research.libd.org/spatialLIBD/reference/layer_stat_cor_plot.html

library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)
library(edgeR)

# mat_formula = ~ 0 + bayesSpace_harmony_9 + region + age + sex 
# foo <- function(sce, var_oi, covars, block_var = NULL) { #must specify in documentation that the second element of the formula is the cluster label, and the first element is zero
#   
#   terms <- attributes(terms(fo))$term.labels
#   pd <- colData(sce)[ , c(terms[!grepl(":", terms)], block_var) ]
#   
#   sce_pseudo <- aggregateAcrossCells(sce, pd)
#   
#   ## then model
#   
#   ## then re-arrange results
# }

var_oi = "bayesSpace_harmony_9"
covars = c("region","age","sex")

#cluster <-  attributes(terms(mat_formula))$term.labels[1]
mat_formula <- eval(str2expression(paste("~","0","+",var_oi,"+",paste(covars, collapse=" + "))))


## Pseudo-bulk for our current BayesSpace cluster results
spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    BayesSpace = colData(spe)[[var_oi]],
    sample_id = spe$sample_id
  )
)


###############################

#find a good expression cutoff using edgeR::filterByExpr https://rdrr.io/bioc/edgeR/man/filterByExpr.html
rowData(spe_pseudo)$filter_expr <- filterByExpr(spe_pseudo)
summary(rowData(spe_pseudo)$filter_expr)
# Mode   FALSE    TRUE 
# logical   21059    7857 

spe_pseudo <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr),]


#spe_pseudo <- logNormCounts(spe_pseudo,size.factors = NULL) #change this to use calcNormFactors from edgeR

#spe_pseudo <- logNormCounts(spe_pseudo,size.factors = NULL) #change this to use calcNormFactors from edgeR
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
## We don't need this 'x' object anymore
rm(x)


##### get mean expression  ####
mat_filter <- assays(spe_pseudo)$logcounts #make matrix of filtered just the log normalized counts


#####################
## Build a group model


#tell user in documentation to make sure their columns are converted to either factors or numerics as appropriate
#convert variables to factors 
# colData(spe_pseudo)[[cluster]] <- as.factor(colData(spe_pseudo)[[cluster]])
# colData(spe_pseudo)$region <- as.factor(colData(spe_pseudo)$region)
# colData(spe_pseudo)$age <- as.numeric(colData(spe_pseudo)$age)
# colData(spe_pseudo)$sex <- as.factor(colData(spe_pseudo)$sex)
# colData(spe_pseudo)$diagnosis <- as.factor(colData(spe_pseudo)$diagnosis)
# colData(spe_pseudo)$subject <- as.factor(colData(spe_pseudo)$subject)

### access different elements of formula and check to see if they're in colData(spe_pseudo)
terms <- attributes(terms(mat_formula))$term.labels
terms <- terms[!grepl(":", terms)]
for(i in seq_along(terms)){
  if(!terms[i] %in% colnames(colData(spe_pseudo))){
    stop("Error: formula term ",terms[i], " is not contained in colData()")
  }
}

#create matrix where the rownames are the sample:clusters and the columns are the other variables (spatial.cluster + region + age + sex)

mod<- model.matrix(mat_formula,
                   data = colData(spe_pseudo)) #binarizes factors 



## get duplicate correlation #http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/dupcor.html
corfit <- duplicateCorrelation(mat_filter, mod,
                               block = spe_pseudo$sample_id)

## Next for each layer test that layer vs the rest
cluster_idx <- splitit(colData(spe_pseudo)[,var_oi]) 

eb0_list_cluster <- lapply(cluster_idx, function(x) {
  res <- rep(0, ncol(spe_pseudo))
  res[x] <- 1
  res_formula <- paste("~","res","+",paste(covars, collapse=" + "))
  m <- with(colData(spe_pseudo),
            model.matrix(eval(str2expression(res_formula)))) 
  
  #josh suggested use top table as a wrapper because it makes the output of eBayes nicer

    eBayes(
    lmFit(
      mat_filter,
      design = m,
      block = spe_pseudo$sample_id,
      correlation = corfit$consensus.correlation
    )
  )
})


##########
## Extract the p-values
pvals0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cluster) = rownames(mat_filter)

t0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cluster) = rownames(mat_filter)
fdrs0_contrasts_cluster = apply(pvals0_contrasts_cluster, 2, p.adjust, 'fdr')

data.frame(
  'FDRsig' = colSums(fdrs0_contrasts_cluster < 0.05 &
                       t0_contrasts_cluster > 0),
  'Pval10-6sig' = colSums(pvals0_contrasts_cluster < 1e-6 &
                            t0_contrasts_cluster > 0),
  'Pval10-8sig' = colSums(pvals0_contrasts_cluster < 1e-8 &
                            t0_contrasts_cluster > 0)
)

# FDRsig Pval10.6sig Pval10.8sig
# 1  12608         178          71
# 2    498         171          68
# 3     52           6           2
# 4  18978       13168        6562
# 5  16379         182          67
# 6      5           0           0
# 7     61           9           3


f_merge <- function(p, fdr, t) { 
  colnames(p) <- paste0('p_value_', colnames(p))
  colnames(fdr) <- paste0('fdr_', colnames(fdr))
  colnames(t) <- paste0('t_stat_', colnames(t))
  res <- as.data.frame(cbind(t, p, fdr))
  res$ensembl <- rownames(res)
  ## Check it's all in order
  res <-merge(x=res,y=rowData(spe_pseudo)[, c("gene_id", "gene_name")],by.x = "ensembl",
              by.y = "gene_id")
  colnames(res)[11] <- "gene"
  rownames(res) <- res$ensembl
  return(res)
}

results_specificity <-
  f_merge(p = pvals0_contrasts_cluster, fdr = fdrs0_contrasts_cluster, t = t0_contrasts_cluster)
head(results_specificity)

#object to return
as.data.frame(results_specificity@listData)

#stop here

modeling_results = fetch_data(type = "modeling_results")

cor <- layer_stat_cor(
  test.df,
  modeling_results,
  model_type = names(modeling_results)[2],
  reverse = FALSE,
  top_n = NULL
)


### heatmap ### here can also use layer_stat_cor_plot() from spatialLIBD
theSeq = seq(-1.0, 1.0, by = 0.01)
my.col <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

dd = dist(1-cor_t_layer)
hc = hclust(dd)
cor_t_layer_toPlot = cor_t_layer[hc$order, c(1, 7:2)]
colnames(cor_t_layer_toPlot) = gsub("ayer", "", colnames(cor_t_layer_toPlot))
rownames(cor_t_layer_toPlot)[rownames(cor_t_layer_toPlot) == "Oligodendrocytes"] = "OLIGO" # does thismatter? 


#implement layer_stat_cor_plot
pdf(file = here::here("plots","07_spatial_registration",paste0("dlpfc_pseudobulked_bayesSpace_vs_mannual_annotations_k",k,".pdf")), width = 8)
print(
  levelplot(
    cor_t_layer_toPlot,
    aspect = "fill",
    at = theSeq,
    col.regions = my.col,
    ylab = "",
    xlab = "",
    scales = list(x = list(rot = 90, cex = 1.5), y = list(cex = 1.5))
  )
)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()