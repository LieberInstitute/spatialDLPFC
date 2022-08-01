library('zellkonverter')
library(SingleCellExperiment)
library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)
library(edgeR)

# define function to compute enrichment statistics 
computeEnrichment <- function(spe, var_oi, covars) {
  mat_formula <- eval(str2expression(paste("~", "0", "+", var_oi, "+", paste(covars, collapse = " + "))))
  
  ## Pseudo-bulk for our current BayesSpace cluster results
  message("Make psuedobulk object")
  spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
      BayesSpace = colData(spe)[[var_oi]],
      sample_id = spe$sample_id
    )
  )
  
  message("Filter lowly expressed genes")
  rowData(spe_pseudo)$filter_expr <- filterByExpr(spe_pseudo)
  summary(rowData(spe_pseudo)$filter_expr)
  # Mode   FALSE    TRUE
  # logical   31057    5544
  
  spe_pseudo <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr), ]
  
  message("Normalize expression")
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
  mat_filter <- assays(spe_pseudo)$logcounts # make matrix of filtered just the log normalized counts
  
  #####################
  ## Build a group model
  
  ### access different elements of formula and check to see if they're in colData(spe_pseudo)
  terms <- attributes(terms(mat_formula))$term.labels
  terms <- terms[!grepl(":", terms)]
  for (i in seq_along(terms)) {
    if (!terms[i] %in% colnames(colData(spe_pseudo))) {
      stop("Error: formula term ", terms[i], " is not contained in colData()")
    }
  }
  
  # create matrix where the rownames are the sample:clusters and the columns are the other variables (spatial.cluster + region + age + sex)
  message("Create model matrix")
  mod <- model.matrix(mat_formula,
                      data = colData(spe_pseudo)
  ) # binarizes factors
  
  
  ## get duplicate correlation #http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/dupcor.html
  message("Run dupllicateCorrelation()")
  corfit <- duplicateCorrelation(mat_filter, mod,
                                 block = spe_pseudo$sample_id
  )
  
  ## Next for each layer test that layer vs the rest
  cluster_idx <- splitit(colData(spe_pseudo)[, var_oi])
  
  message("Run enrichment statistics")
  eb0_list_cluster <- lapply(cluster_idx, function(x) {
    res <- rep(0, ncol(spe_pseudo))
    res[x] <- 1
    res_formula <- paste("~", "res", "+", paste(covars, collapse = " + "))
    m <- with(
      colData(spe_pseudo),
      model.matrix(eval(str2expression(res_formula)))
    )
    
    eBayes(
      lmFit(
        mat_filter,
        design = m,
        block = spe_pseudo$sample_id,
        correlation = corfit$consensus.correlation
      )
    )
  })
  
  
  message("extract and reformat enrichment results")
  ##########
  ## Extract the p-values
  pvals0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
    x$p.value[, 2, drop = FALSE]
  })
  rownames(pvals0_contrasts_cluster) <- rownames(mat_filter)
  
  t0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
    x$t[, 2, drop = FALSE]
  })
  rownames(t0_contrasts_cluster) <- rownames(mat_filter)
  
  fdrs0_contrasts_cluster <- apply(pvals0_contrasts_cluster, 2, p.adjust, "fdr")
  rownames(fdrs0_contrasts_cluster) <- rownames(mat_filter)
  
  data.frame(
    "FDRsig" = colSums(fdrs0_contrasts_cluster < 0.05 &
                         t0_contrasts_cluster > 0),
    "Pval10-6sig" = colSums(pvals0_contrasts_cluster < 1e-6 &
                              t0_contrasts_cluster > 0),
    "Pval10-8sig" = colSums(pvals0_contrasts_cluster < 1e-8 &
                              t0_contrasts_cluster > 0)
  )
  
  f_merge <- function(p, fdr, t) {
    colnames(p) <- paste0("p_value_", colnames(p))
    colnames(fdr) <- paste0("fdr_", colnames(fdr))
    colnames(t) <- paste0("t_stat_", colnames(t))
    res <- as.data.frame(cbind(t, p, fdr))
    res$ensembl <- rownames(res)
    ## Check it's all in order
    res <- merge(
      x = res, y = rowData(spe_pseudo)[, c("gene_id", "gene_name")], by.x = "ensembl",
      by.y = "gene_id"
    )
    colnames(res)[11] <- "gene"
    rownames(res) <- res$ensembl
    return(res)
  }
  
  results_specificity <-
    f_merge(p = pvals0_contrasts_cluster, fdr = fdrs0_contrasts_cluster, t = t0_contrasts_cluster)
  head(results_specificity)
  
  # object to return
  results_specificity <- as.data.frame(results_specificity@listData)
  rownames(results_specificity) <- results_specificity$ensembl
  message("return enrichment results")
  return(results_specificity)
}

#load annotated snRNAseq data
sce = readH5AD('/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version2/CMC/CMC-snRNAseq_annotated.h5ad')

# identify annotation/cluster labels
var_oi <- "anno"
covars <- c("demux_type")

rowData(sce)$gene_name <- rownames(rowData(sce)) #save gene name as column of rowData
rownames(sce) <- rowData(sce)$featureid # have to make row names of object the ensembl id instead of gene names
colData(sce)$anno <- as.factor(colData(sce)$anno)
colData(sce)$demux_type <- as.factor(colData(sce)$demux_type)
colnames(colData(sce))[9] <- "sample_id"
names(assays(sce)) <- "counts"

results_specificity <- computeEnrichment(sce, var_oi, covars)
save(results_specificity, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/12_spatial_registration_sc/results_specificity.RDS")

specificity_stats <- results_specificity[, grep("^t_stat", colnames(results_specificity))]
colnames(specificity_stats) <- gsub("^t_stat_", "", colnames(specificity_stats))

# vs manual annotations
modeling_results <- fetch_data(type = "modeling_results")
cor <- layer_stat_cor(
  specificity_stats,
  modeling_results,
  model_type = names(modeling_results)[2],
  reverse = FALSE,
  top_n = NULL
)


pdf("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/12_spatial_registration_sc/spatial_registration_plot_sc_v_manual.pdf")
layer_stat_cor_plot(cor, max = 0.7)
dev.off()


