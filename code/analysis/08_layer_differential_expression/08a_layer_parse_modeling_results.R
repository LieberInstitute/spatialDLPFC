# library(sgejobs)
# sgejobs::job_loop(
#    loops = list(spetype = c(
#        "wholegenome", "targeted"
#     )),
#     name = "04_parse_model_results",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "15G")
# To execute the script builder, use: sh 04_parse_model_results.sh

# # Required libraries
# library("getopt")
# 
# ## Specify parameters
# spec <- matrix(c(
#   "spetype", "s", 2, "character", "SPE spetype: wholegenome or targeted",
#   "help", "h", 0, "logical", "Display help"
# ), byrow = TRUE, ncol = 5)
# opt <- getopt(spec = spec)
# 
# ## if help was asked for print a friendly message
# ## and exit with a non-zero error code
# if (!is.null(opt$help)) {
#   cat(getopt(spec, usage = TRUE))
#   q(status = 1)
# }
# 
# ## For testing
# if (FALSE) {
#   opt <- list(spetype = "wholegenome")
# }


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
load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("cluster_modeling_results_k",k,".Rdata")),verbose = TRUE)
#load(file = "processed-data/rdata/spe/08_layer_differential_expression/cluster_modeling_results_k9.Rdata", verbose = TRUE)

## load spe data
load(file = here::here("processed-data","rdata","spe","pseudo_bulked_spe",paste0("sce_pseudobulk_bayesSpace_normalized_filtered_k",k,".Rdata")))


## We don't want to model the pathology groups as integers / numeric
## so let's double check this
#stopifnot(is.factor(spe_pseudo$path_groups) || is.character(spe_pseudo$path_groups))

## Add APOe genotype info
#spe_pseudo$APOe <- c("Br3854" = "E3/E4", "Br3873" = "E3/E3", "Br3880" = "E3/E3", "Br3874" = "E2/E3")[spe_pseudo$subject]

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
#    Ab+             both          next_Ab+
# Mode :logical   Mode :logical   Mode :logical
# FALSE:27853     FALSE:27852     FALSE:27853
#                 TRUE :1
# next_both        next_pT+          none
# Mode :logical   Mode :logical   Mode :logical
# FALSE:27852     FALSE:27850     FALSE:27836
# TRUE :1         TRUE :3         TRUE :17
#    pT+
# Mode :logical
# FALSE:27764
# TRUE :89

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
#   t_stat_Ab+ t_stat_both t_stat_next_Ab+ t_stat_next_both t_stat_next_pT+ t_stat_none t_stat_pT+ p_value_Ab+ p_value_both p_value_next_Ab+ p_value_next_both p_value_next_pT+ p_value_none p_value_pT+   fdr_Ab+ fdr_both fdr_next_Ab+ fdr_next_both fdr_next_pT+  fdr_none   fdr_pT+         ensembl        gene
# 1 -0.4244225  -0.3072583       2.5699207       -0.4244000      -0.4244020 -0.42447838 -0.4244722   0.6732106    0.7600156       0.01343171         0.6732269       0.67322546    0.6731701  0.67317467 0.9707058        1    0.7741261     0.9950277    1.0000000 0.8585820 1.0000000 ENSG00000243485 MIR1302-2HG
# 2 -1.1039826   1.0202366       0.3772734       -0.6221515       0.6280650 -0.74566366  0.5061814   0.2752682    0.3128833       0.70768300         0.5368743       0.53302712    0.4596151  0.61511557 0.9707058        1    0.9710289     0.9950277    1.0000000 0.7952753 1.0000000 ENSG00000238009  AL627309.1
# 3 -0.3764208  -0.3554921      -0.3764600       -0.3762571      -0.3762774 -0.37680277  2.3400513   0.7083124    0.7238250       0.70828341         0.7084333       0.70841829    0.7080304  0.02361651 0.9707058        1    0.9710289     0.9950277    1.0000000 0.8585820 0.3793487 ENSG00000239945  AL627309.3
# 4 -0.1869257   0.3737276      -1.2173788       -0.6407296       1.2046181  0.07293914  0.4157680   0.8525302    0.7103018       0.22958845         0.5248365       0.23443258    0.9421667  0.67948807 1.0000000        1    0.9710289     0.9950277    1.0000000 1.0000000 1.0000000 ENSG00000241860  AL627309.5
# 5 -0.5067356   0.7722787      -0.4854839        0.7400436      -0.9161310  1.30626123 -0.8452318   0.6147296    0.4438482       0.62960767         0.4629854       0.36431449    0.1978740  0.40230009 0.9707058        1    0.9710289     0.9950277    1.0000000 0.7128183 1.0000000 ENSG00000237491   LINC01409
# 6 -0.4097255  -0.3869820      -0.4097330       -0.4096966       2.5711612 -0.40979604 -0.4097907   0.6838845    0.7005310       0.68387905         0.6839056       0.01338978    0.6838331  0.68383709 0.9707058        1    0.9710289     0.9950277    0.6833824 0.8585820 1.0000000 ENSG00000177757      FAM87B

## From
## https://github.com/LieberInstitute/HumanPilot/blob/879f11c7c57efd72334332b40feb3ad623e067c8/Analysis/Layer_Guesses/misc_numbers.R#L259-L264
pvals_contrasts <- eb_contrasts$p.value
fdrs_contrasts <- apply(pvals_contrasts, 2, p.adjust, "fdr")
dim(pvals_contrasts)
# [1] 27853    21
summary(fdrs_contrasts < 0.05)
# Abpos-both      Abpos-next_Abpos Abpos-next_both Abpos-next_pTpos
#  Mode :logical   Mode :logical    Mode :logical   Mode :logical
#  FALSE:27853     FALSE:27853      FALSE:27853     FALSE:27853
#
#  Abpos-none      Abpos-pTpos     both-next_Abpos both-next_both
#  Mode :logical   Mode :logical   Mode :logical   Mode :logical
#  FALSE:27852     FALSE:27850     FALSE:27853     FALSE:27853
#  TRUE :1         TRUE :3
#  both-next_pTpos both-none       both-pTpos      next_Abpos-next_both
#  Mode :logical   Mode :logical   Mode :logical   Mode :logical
#  FALSE:27853     FALSE:27850     FALSE:27852     FALSE:27853
#                  TRUE :3         TRUE :1
#  next_Abpos-next_pTpos next_Abpos-none next_Abpos-pTpos
#  Mode :logical         Mode :logical   Mode :logical
#  FALSE:27853           FALSE:27852     FALSE:27853
#                        TRUE :1
#  next_both-next_pTpos next_both-none  next_both-pTpos next_pTpos-none
#  Mode :logical        Mode :logical   Mode :logical   Mode :logical
#  FALSE:27853          FALSE:27852     FALSE:27850     FALSE:27852
#                       TRUE :1         TRUE :3         TRUE :1
#  next_pTpos-pTpos none-pTpos
#  Mode :logical    Mode :logical
#  FALSE:27850      FALSE:27850
#  TRUE :3          TRUE :3

results_pairwise <-
  f_merge(p = pvals_contrasts, fdr = fdrs_contrasts, t = eb_contrasts$t)
colnames(results_pairwise)
sort(colSums(fdrs_contrasts < 0.05))
#       Abpos-both      Abpos-next_Abpos       Abpos-next_both
#                0                     0                     0
# Abpos-next_pTpos       both-next_Abpos        both-next_both
#                0                     0                     0
#  both-next_pTpos  next_Abpos-next_both next_Abpos-next_pTpos
#                0                     0                     0
# next_Abpos-pTpos  next_both-next_pTpos            Abpos-none
#                0                     0                     1
#       both-pTpos       next_Abpos-none        next_both-none
#                1                     1                     1
#  next_pTpos-none           Abpos-pTpos             both-none
#                1                     3                     3
#  next_both-pTpos      next_pTpos-pTpos            none-pTpos
#                3                     3                     3

## From
## https://github.com/LieberInstitute/HumanPilot/blob/879f11c7c57efd72334332b40feb3ad623e067c8/Analysis/Layer_Guesses/misc_numbers.R#L205-L215
f_sig <- function(type, cut = 0.05) {
  cbind(
    "n" = addmargins(table(f_stats[[type]] < cut)),
    "ratio" = addmargins(table(f_stats[[type]] < cut)) / nrow(f_stats)
  )
}
f_sig("noWM_fdr")
#           n        ratio
# FALSE 27851 9.999282e-01
# TRUE      2 7.180555e-05
# Sum   27853 1.000000e+00

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
#   f_stat_noWM p_value_noWM  fdr_noWM noWM_AveExpr         ensembl        gene
# 1   0.9810021    0.4502729 0.9763246    0.6515295 ENSG00000243485 MIR1302-2HG
# 2   0.5485499    0.7681921 0.9921032    1.1652509 ENSG00000238009  AL627309.1
# 3   0.8062886    0.5708651 0.9763246    0.6407770 ENSG00000239945  AL627309.3
# 4   0.4809293    0.8186675 1.0000000    1.9572611 ENSG00000241860  AL627309.5
# 5   0.6429322    0.6952897 0.9801360    2.9180726 ENSG00000237491   LINC01409
# 6   0.9722365    0.4559056 0.9763246    0.6494224 ENSG00000177757      FAM87B

modeling_results <- list(
  "anova" = results_anova,
  "enrichment" = results_specificity,
  "pairwise" = results_pairwise
)

save(modeling_results,file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("parsed_modeling_results_k",k,".Rdata")))

which(modeling_results$anova$fdr_noWM < 0.05)
summary(modeling_results$anova$fdr_noWM)

length(which(modeling_results$enrichment$fdr_1 < 0.05))
length(which(modeling_results$enrichment$fdr_2 < 0.05))
length(which(modeling_results$enrichment$fdr_3 < 0.05))
length(which(modeling_results$enrichment$fdr_4 < 0.05))
length(which(modeling_results$enrichment$fdr_5 < 0.05))
length(which(modeling_results$enrichment$fdr_6 < 0.05))
length(which(modeling_results$enrichment$fdr_7 < 0.05))
length(which(modeling_results$enrichment$fdr_8 < 0.05))
length(which(modeling_results$enrichment$fdr_9 < 0.05))

cluster <- c(1:9)
genes <- c(9255,4051,1293,4739,2615,4422,3250,3332,417)
df <- data.frame(cluster, genes)
pdf(file = here::here("plots","08_layer_differential_expression","plot_enrichment_DEGs.pdf"))
plot(df$genes~df$cluster)
dev.off()


library(vioplot)
pdf(file = here::here("plots","08_layer_differential_expression","boxplot_num_enrichment_DEGs.pdf"))
x1 <- modeling_results$enrichment$fdr_1
x2 <- modeling_results$enrichment$fdr_2
x3 <- modeling_results$enrichment$fdr_3
x4 <- modeling_results$enrichment$fdr_4
x5 <- modeling_results$enrichment$fdr_5
x6 <- modeling_results$enrichment$fdr_6
x7 <- modeling_results$enrichment$fdr_7
x8 <- modeling_results$enrichment$fdr_8
x9 <- modeling_results$enrichment$fdr_9

vioplot(x1,x2,x3,x4,x5,x6,x7,x8,x9)
dev.off()


which(modeling_results$enrichment$p_value_1 < 0.05)
summary(modeling_results$pa$fdr_1)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()