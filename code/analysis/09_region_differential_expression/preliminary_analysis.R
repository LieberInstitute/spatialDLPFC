library(SpatialExperiment)
library(spatialLIBD)
library(here)
library(edgeR)
library(scuttle)
#library(scRNAseq)
library(scater)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load spe object and clusters
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final.Rdata"),verbose = TRUE)
spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"), 
  prefix = ""
)

## Pseudo-bulk for our current BayesSpace cluster results
spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    BayesSpace = colData(spe)[[paste0("bayesSpace_harmony_",k)]],
    sample_id = spe$sample_id
  )
)
spe_pseudo$BayesSpace <- factor(spe_pseudo$BayesSpace)
#save here for differential expression analysis
save(spe_pseudo, file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("sce_pseudobulk_bayesSpace_k",k,".Rdata")))

#find a good expression cutoff using edgeR::filterByExpr https://rdrr.io/bioc/edgeR/man/filterByExpr.html
rowData(spe_pseudo)$low_expr <- filterByExpr(spe_pseudo) #add group by region keep <- filterByExpr(y, group=current$tomato)
summary(rowData(spe_pseudo)$low_expr)
# Mode   FALSE    TRUE 
# logical   21059    7857 

spe_pseudo <- spe_pseudo[which(!rowData(spe_pseudo)$low_expr),]
dim(spe_pseudo)
#[1] 21059    90

#log normalize the counts
# spe_pseudo <- logNormCounts(spe_pseudo,size.factors = NULL)
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
## We don't need this 'x' object anymore
rm(x)

dim(spe_pseudo)
#[1] 28916    90

#update this to indicate this version of the object is normalized and filtered
save(spe_pseudo, file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("sce_pseudobulk_bayesSpace_normalized_filtered_k",k,".Rdata")))

#run PCA
pca <- prcomp(t(assays(spe_pseudo)$logcounts))
jaffelab::getPcaVars(pca)[seq_len(50)]
#  [1] 14.600  4.520  3.150  1.390  1.220  1.080  0.942  0.929  0.889  0.873
# [11]  0.858  0.855  0.845  0.831  0.818  0.810  0.800  0.796  0.787  0.785
# [21]  0.774  0.768  0.765  0.754  0.748  0.742  0.737  0.725  0.720  0.716
# [31]  0.706  0.704  0.695  0.691  0.682  0.671  0.670  0.666  0.657  0.647
# [41]  0.641  0.628  0.625  0.615  0.606  0.605  0.595  0.584  0.578  0.575
pca_pseudo<- pca$x[, seq_len(50)]



reducedDims(spe_pseudo) <- list(PCA=pca_pseudo) 

# code adapted from: http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html#2_Diagnostic_plots_for_quality_control
###plot PCA###
pdf(file = here::here("plots","09_region_differential_expression",paste0("sce_pseudobulk_pca_k",k,".pdf")), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "subject", ncomponents = 12, point_size = 1) 
plotPCA(spe_pseudo, colour_by = "region", ncomponents = 12, point_size = 1) 
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 12, point_size = 1) 
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 12, point_size = 1) 
plotPCA(spe_pseudo, colour_by = "BayesSpace", ncomponents = 12, point_size = 1) 
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 12, point_size = 1)
dev.off()

####plot explanatory variables ####

#uses linear regression model
vars <- getVarianceExplained(spe_pseudo, 
                             variables=c("subject", "region", "sex", "age", "BayesSpace","sample_id")) 
head(vars)

# subject    region       sex         age
# ENSG00000243485  6.149209 0.8049553 0.9891498 0.684494271
# ENSG00000238009 22.059419 0.3825081 6.0221979 0.758922808
# ENSG00000239945 10.112360 2.2471910 0.7490637 0.006556772
# ENSG00000241860 15.910843 1.7697885 2.7770465 3.317220956
# ENSG00000229905 10.113068 0.5911937 1.0110245 0.266643458
# ENSG00000237491 17.499031 9.9500577 0.2762791 0.570564318

pdf(file = here::here("plots","09_region_differential_expression",paste0("plot_explanatory_vars_k",k,"_testLeo.pdf")))
plotExplanatoryVariables(vars)
dev.off()
