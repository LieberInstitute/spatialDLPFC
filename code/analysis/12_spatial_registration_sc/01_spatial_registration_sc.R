library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)
library(edgeR)

computeEnrichment <- function(spe,var_oi,covars){
  mat_formula <- eval(str2expression(paste("~","0","+",var_oi,"+",paste(covars, collapse=" + "))))
  
  ## Pseudo-bulk for our current BayesSpace cluster results
  message("Make psuedobulk object")
  spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
      BayesSpace = colData(spe)[[var_oi]],
      sample_id = spe$sample_id
    )
  )

  # class: SingleCellExperiment 
  # dim: 36601 483 
  # metadata(1): Samples
  # assays(1): counts
  # rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
  # rowData names(7): source type ... gene_type binomial_deviance
  # colnames: NULL
  # colData names(33): sample_id Barcode ... sample_id ncells
  # reducedDimNames(4): GLMPCA_approx TSNE UMAP HARMONY
  # mainExpName: NULL
  # altExpNames(0):
 
  ###############################
  message("Filter lowly expressed genes")
  rowData(spe_pseudo)$filter_expr <- filterByExpr(spe_pseudo)
  summary(rowData(spe_pseudo)$filter_expr)
  # Mode   FALSE    TRUE 
  # logical   31057    5544 
  
  spe_pseudo <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr),]

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
  mat_filter <- assays(spe_pseudo)$logcounts #make matrix of filtered just the log normalized counts
  
  
  #####################
  ## Build a group model
  
  ### access different elements of formula and check to see if they're in colData(spe_pseudo)
  terms <- attributes(terms(mat_formula))$term.labels
  terms <- terms[!grepl(":", terms)]
  for(i in seq_along(terms)){
    if(!terms[i] %in% colnames(colData(spe_pseudo))){
      stop("Error: formula term ",terms[i], " is not contained in colData()")
    }
  }
  
  #create matrix where the rownames are the sample:clusters and the columns are the other variables (spatial.cluster + region + age + sex)
  message("Create model matrix")
  mod<- model.matrix(mat_formula,
                     data = colData(spe_pseudo)) #binarizes factors
  
  
  
  ## get duplicate correlation #http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/dupcor.html
  corfit <- duplicateCorrelation(mat_filter, mod,
                                 block = spe_pseudo$sample_id)
  
  ## Next for each layer test that layer vs the rest
  cluster_idx <- splitit(colData(spe_pseudo)[,var_oi])
  
  message("run enrichment statistics")
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
  
  
  message("extract and reformat enrichment results")
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
  rownames(fdrs0_contrasts_cluster) = rownames(mat_filter)
  
  data.frame(
    'FDRsig' = colSums(fdrs0_contrasts_cluster < 0.05 &
                         t0_contrasts_cluster > 0),
    'Pval10-6sig' = colSums(pvals0_contrasts_cluster < 1e-6 &
                              t0_contrasts_cluster > 0),
    'Pval10-8sig' = colSums(pvals0_contrasts_cluster < 1e-8 &
                              t0_contrasts_cluster > 0)
  )
  
  # FDRsig Pval10.6sig Pval10.8sig
  # Astro            992         187         107
  # Endo.Mural_01    448          71          36
  # Endo.Mural_02    909         106          41
  # Excit_01         382          19           6
  # Excit_02         205           4           1
  # Excit_03         133          10           3
  # Excit_04         131           9           1
  # Excit_05          20           0           0
  # Excit_06         293          32           9
  # Excit_07         190          16           7
  # Excit_08         253          22           2
  # Excit_09        1359         238          63
  # Excit_10        1354         392         224
  # Excit_11           0           0           0
  # Excit_12           0           0           0
  # Excit_13         234          13           3
  # Excit_14          17           0           0
  # Excit_15          14           0           0
  # Inhib_01          17           2           0
  # Inhib_02         287          37          14
  # Inhib_03          48           3           0
  # Inhib_04         381          73          35
  # Inhib_05         254          26           6
  # Inhib_06         183          32           6
  # Micro           1341         245         147
  # Oligo_01         439          11           8
  # Oligo_02        1291         296         178
  # Oligo_03         731         128          66
  # OPC              849         134          70
  
  
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
  results_specificity <-as.data.frame(results_specificity@listData)
  rownames(results_specificity) <- results_specificity$ensembl
  message("return enrichment results")
  return(results_specificity)
}

#load sc data
load(file = "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata")

var_oi = "cellType_hc"
covars = c("region","age","sex")

rownames(sce)<-rowData(sce)$gene_id #have to make row names of object the ensembl id instead of gene names
colData(sce)$region <- as.factor(colData(sce)$region)
colData(sce)$age <- as.numeric(colData(sce)$age)
colData(sce)$sex <- as.factor(colData(sce)$sex)
colnames(colData(sce))[1] <- "sample_id"

results_specificity <- computeEnrichment(sce,var_oi,covars)

specificity_stats <- results_specificity[, grep("^t_stat", colnames(results_specificity))]

modeling_results = fetch_data(type = "modeling_results")
cor <- layer_stat_cor(
  specificity_stats,
  modeling_results,
  model_type = names(modeling_results)[2],
  reverse = FALSE,
  top_n = NULL
)

pdf("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/12_spatial_registration_sc")
layer_stat_cor_plot(cor, max = 1)
dev.off()
