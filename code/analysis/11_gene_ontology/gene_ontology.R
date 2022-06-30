################################################
# spatial_DG_lifespan project
# GO enrichment & plotting from modeling results
# Anthony Ramnauth, June 08 2022
################################################

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(SingleCellExperiment)
  library(dplyr)
  library(sessioninfo)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ggplot2)
})

modeling_results <- readRDS(file = here::here("processed-data","pseudobulk_spe","modeling_results.rds"))

enriched_model <- modeling_results$enrichment

## Create named vector of gene lists with ENTREZID (works best) for GO (ranking by t-statistic)

list_1 <- enriched_model %>%
  dplyr::select(gene, ensembl, fdr_1, t_stat_1)

list_2 <- enriched_model %>%
  dplyr::select(gene, ensembl, fdr_2, t_stat_2)

list_3 <- enriched_model %>%
  dplyr::select(gene, ensembl, fdr_3, t_stat_3)

list_4 <- enriched_model %>%
  dplyr::select(gene, ensembl, fdr_4, t_stat_4)

list_5 <- enriched_model %>%
  dplyr::select(gene, ensembl, fdr_5, t_stat_5)

list_6 <- enriched_model %>%
  dplyr::select(gene, ensembl, fdr_6, t_stat_6)

list_7 <- enriched_model %>%
  dplyr::select(gene, ensembl, fdr_7, t_stat_7)

list_8 <- enriched_model %>%
  dplyr::select(gene, ensembl, fdr_8, t_stat_8)

top_list_1 <- list_1 %>%
  filter(fdr_1 < 0.05) %>%
  dplyr::arrange(desc(t_stat_1)) %>%
  dplyr::select(gene, t_stat_1)

## feature 1: numeric vector
clust_1 = top_list_1[,2]
## feature 2: named vector
names(clust_1) = as.character(top_list_1[,1])
# final vector
clust_1 <- names(clust_1)
clust_1 = bitr(clust_1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(clust_1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   23.08% of input gene IDs are fail to map...

top_list_2 <- list_2 %>%
  filter(fdr_2 < 0.05) %>%
  dplyr::arrange(desc(t_stat_2)) %>%
  dplyr::select(gene, t_stat_2)

## feature 1: numeric vector
clust_2 = top_list_2[,2]
## feature 2: named vector
names(clust_2) = as.character(top_list_2[,1])
# final vector
clust_2 <- names(clust_2)
clust_2 = bitr(clust_2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(clust_2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.46% of input gene IDs are fail to map...

top_list_3 <- list_3 %>%
  filter(fdr_3 < 0.05) %>%
  dplyr::arrange(desc(t_stat_3)) %>%
  dplyr::select(gene, t_stat_3)

## feature 1: numeric vector
GCL = top_list_3[,2]
## feature 2: named vector
names(GCL) = as.character(top_list_3[,1])
# final vector
GCL <- names(GCL)
GCL = bitr(GCL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(GCL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.86% of input gene IDs are fail to map...

top_list_4 <- list_4 %>%
  filter(fdr_4 < 0.05) %>%
  dplyr::arrange(desc(t_stat_4)) %>%
  dplyr::select(gene, t_stat_4)

## feature 1: numeric vector
SGZ = top_list_4[,2]
## feature 2: named vector
names(SGZ) = as.character(top_list_4[,1])
# final vector
SGZ <- names(SGZ)
SGZ = bitr(SGZ, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

top_list_5 <- list_5 %>%
  filter(fdr_5 < 0.05) %>%
  dplyr::arrange(desc(t_stat_5)) %>%
  dplyr::select(gene, t_stat_5)

## feature 1: numeric vector
CA4 = top_list_5[,2]
## feature 2: named vector
names(CA4) = as.character(top_list_5[,1])
# final vector
CA4 <- names(CA4)
CA4 = bitr(CA4, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(CA4, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   1.79% of input gene IDs are fail to map...

top_list_6 <- list_6 %>%
  filter(fdr_6 < 0.05) %>%
  dplyr::arrange(desc(t_stat_6)) %>%
  dplyr::select(gene, t_stat_6)

## feature 1: numeric vector
CA3 = top_list_6[,2]
## feature 2: named vector
names(CA3) = as.character(top_list_6[,1])
# final vector
CA3 <- names(CA3)
CA3 = bitr(CA3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(CA3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.84% of input gene IDs are fail to map...

list_7 <- list_7 %>%
  arrange(desc(t_stat_7))
# BayesSpace cluster 7 (Molec. Layer) has only 1 gene that is FDR<0.05, so do top 20 instead
top_list_7 <- head(list_7, 20) %>%
  dplyr::select(gene, t_stat_7)

## feature 1: numeric vector
ML = top_list_7[,2]
## feature 2: named vector
names(ML) = as.character(top_list_7[,1])
# final vector
ML <- names(ML)
ML = bitr(ML, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(ML, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   5% of input gene IDs are fail to map...

top_list_8 <- list_8 %>%
  filter(fdr_8 < 0.05) %>%
  dplyr::arrange(desc(t_stat_8)) %>%
  dplyr::select(gene, t_stat_8)

## feature 1: numeric vector
clust_8 = top_list_8[,2]
## feature 2: named vector
names(clust_8) = as.character(top_list_8[,1])
# final vector
clust_8 <- names(clust_8)
clust_8 = bitr(clust_8, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(clust_8, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.48% of input gene IDs are fail to map...

# Make a list for comparing cluster ontologies
clust_compare <- list(clust_1$ENTREZID, clust_2$ENTREZID, GCL$ENTREZID, SGZ$ENTREZID,
                      CA4$ENTREZID, CA3$ENTREZID, ML$ENTREZID, clust_8$ENTREZID)

names(clust_compare) <- c("clust_1", "clust_2", "GCL", "SGZ", "CA4", "CA3", "ML", "clust_8")

## GO Gene Over Representation Analysis for enrichment in sub-ontologies

clust_1_CC <- enrichGO(gene = clust_1$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

clust_1_MF <- enrichGO(gene = clust_1$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

clust_1_BP <- enrichGO(gene = clust_1$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

clust_2_CC <- enrichGO(gene = clust_2$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

clust_2_MF <- enrichGO(gene = clust_2$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

clust_2_BP <- enrichGO(gene = clust_2$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

GCL_CC <- enrichGO(gene = GCL$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

GCL_MF <- enrichGO(gene = GCL$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

GCL_BP <- enrichGO(gene = GCL$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

SGZ_CC <- enrichGO(gene = SGZ$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

SGZ_MF <- enrichGO(gene = SGZ$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

SGZ_BP <- enrichGO(gene = SGZ$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

CA4_CC <- enrichGO(gene = CA4$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

CA4_MF <- enrichGO(gene = CA4$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

CA4_BP <- enrichGO(gene = CA4$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

CA3_CC <- enrichGO(gene = CA3$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

CA3_MF <- enrichGO(gene = CA3$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

CA3_BP <- enrichGO(gene = CA3$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

ML_CC <- enrichGO(gene = ML$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "fdr",
                  pvalueCutoff = 0.05,
                  qvalueCutoff  = 0.05,
                  readable = TRUE,
)

ML_MF <- enrichGO(gene = ML$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "fdr",
                  pvalueCutoff = 0.05,
                  qvalueCutoff  = 0.05,
                  readable = TRUE,
)

ML_BP <- enrichGO(gene = ML$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "fdr",
                  pvalueCutoff = 0.05,
                  qvalueCutoff  = 0.05,
                  readable = TRUE,
)

clust_8_CC <- enrichGO(gene = clust_8$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

clust_8_MF <- enrichGO(gene = clust_8$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

clust_8_BP <- enrichGO(gene = clust_8$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

# Create Comparative Gene Over Representation Analysis for enrichment in sub-ontologies

comp_CC <- compareCluster(clust_compare,
                          OrgDb = org.Hs.eg.db,
                          fun = enrichGO,
                          ont = "CC",
                          pAdjustMethod = "fdr",
                          pvalueCutoff=0.05,
                          qvalueCutoff  = 0.05,
                          readable = TRUE,
)

comp_MF <- compareCluster(clust_compare,
                          OrgDb = org.Hs.eg.db,
                          fun = enrichGO,
                          ont = "MF",
                          pAdjustMethod = "fdr",
                          pvalueCutoff=0.05,
                          qvalueCutoff  = 0.05,
                          readable = TRUE,
)

comp_BP <- compareCluster(clust_compare,
                          OrgDb = org.Hs.eg.db,
                          fun = enrichGO,
                          ont = "BP",
                          pAdjustMethod = "fdr",
                          pvalueCutoff=0.05,
                          qvalueCutoff  = 0.05,
                          readable = TRUE,
)

save(clust_1_CC, clust_1_MF, clust_1_BP, clust_2_CC, clust_2_MF, clust_2_BP,
     GCL_CC, GCL_MF, GCL_BP, SGZ_CC, SGZ_MF, SGZ_BP, CA4_CC, CA4_MF, CA4_BP,
     CA3_CC, CA3_MF, CA3_BP, ML_CC, ML_MF, ML_BP, clust_8_CC, clust_8_MF, clust_8_BP,
     comp_CC, comp_MF, comp_BP,
     file = here::here("processed-data","pseudobulk_spe", "gene_ontologies", "All_samples_enrichedGO.Rdata"))

## Barplots of subontologies

options(ggrepel.max.overlaps = Inf)

# For Cluster 1
pdf(file = here::here("plots","pseudobulked","clust_1_GO.pdf"), width = 16, height = 10)

barplot(clust_1_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_1 Cellular Compartment")
barplot(clust_1_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_1 Molecular Function")
barplot(clust_1_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_1 Biological Process")

clust_1_CCx <- setReadable(clust_1_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_1_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_1 Cellular Compartment Network")
clust_1_MFx <- setReadable(clust_1_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_1_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_1 Molecular Function Network")
clust_1_BPx <- setReadable(clust_1_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_1_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_1 Biological Process Network")

clust_1_CCp <- pairwise_termsim(clust_1_CC)
emapplot(clust_1_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_1 Cellular Compartment Modules")
clust_1_MFp <- pairwise_termsim(clust_1_MF)
emapplot(clust_1_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_1 Molecular Function Modules")
clust_1_BPp <- pairwise_termsim(clust_1_BP)
emapplot(clust_1_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_1 Biological Process Modules")

dev.off()

# For Cluster 2
pdf(file = here::here("plots","pseudobulked","clust_2_GO.pdf"), width = 16, height = 10)

barplot(clust_2_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_2 Cellular Compartment")
barplot(clust_2_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_2 Molecular Function")
barplot(clust_2_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_2 Biological Process")

clust_2_CCx <- setReadable(clust_2_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_2_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_2 Cellular Compartment Network")
clust_2_MFx <- setReadable(clust_2_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_2_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_2 Molecular Function Network")
clust_2_BPx <- setReadable(clust_2_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_2_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_2 Biological Process Network")

clust_2_CCp <- pairwise_termsim(clust_2_CC)
emapplot(clust_2_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_2 Cellular Compartment Modules")
clust_2_MFp <- pairwise_termsim(clust_2_MF)
emapplot(clust_2_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_2 Molecular Function Modules")
clust_2_BPp <- pairwise_termsim(clust_2_BP)
emapplot(clust_2_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_2 Biological Process Modules")

dev.off()

# For GCL

pdf(file = here::here("plots","pseudobulked","GCL_GO.pdf"), width = 16, height = 10)

barplot(GCL_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("GCL Cellular Compartment")
barplot(GCL_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("GCL Molecular Function")
barplot(GCL_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("GCL Biological Process")

GCL_CCx <- setReadable(GCL_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GCL_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("GCL Cellular Compartment Network")
GCL_MFx <- setReadable(GCL_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GCL_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("GCL Molecular Function Network")
GCL_BPx <- setReadable(GCL_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(GCL_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("GCL Biological Process Network")

GCL_CCp <- pairwise_termsim(GCL_CC)
emapplot(GCL_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("GCL Cellular Compartment Modules")
GCL_MFp <- pairwise_termsim(GCL_MF)
emapplot(GCL_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("GCL Molecular Function Modules")
GCL_BPp <- pairwise_termsim(GCL_BP)
emapplot(GCL_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("GCL Biological Process Modules")

dev.off()

# For SGZ

pdf(file = here::here("plots","pseudobulked","SGZ_GO.pdf"), width = 16, height = 10)

barplot(SGZ_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("SGZ Cellular Compartment")
barplot(SGZ_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("SGZ Molecular Function")
barplot(SGZ_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("SGZ Biological Process")

SGZ_CCx <- setReadable(SGZ_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(SGZ_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("SGZ Cellular Compartment Network")
SGZ_MFx <- setReadable(SGZ_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(SGZ_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("SGZ Molecular Function Network")
SGZ_BPx <- setReadable(SGZ_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(SGZ_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("SGZ Biological Process Network")

SGZ_CCp <- pairwise_termsim(SGZ_CC)
emapplot(SGZ_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("SGZ Cellular Compartment Modules")
SGZ_MFp <- pairwise_termsim(SGZ_MF)
emapplot(SGZ_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("SGZ Molecular Function Modules")
SGZ_BPp <- pairwise_termsim(SGZ_BP)
emapplot(SGZ_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("SGZ Biological Process Modules")

dev.off()

# For CA4

pdf(file = here::here("plots","pseudobulked","CA4_GO.pdf"), width = 16, height = 10)

barplot(CA4_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("CA4 Cellular Compartment")
barplot(CA4_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("CA4 Molecular Function")
barplot(CA4_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("CA4 Biological Process")

CA4_CCx <- setReadable(CA4_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(CA4_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("CA4 Cellular Compartment Network")
CA4_MFx <- setReadable(CA4_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(CA4_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("CA4 Molecular Function Network")
CA4_BPx <- setReadable(CA4_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(CA4_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("CA4 Biological Process Network")

CA4_CCp <- pairwise_termsim(CA4_CC)
emapplot(CA4_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("CA4 Cellular Compartment Modules")
CA4_MFp <- pairwise_termsim(CA4_MF)
emapplot(CA4_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("CA4 Molecular Function Modules")
CA4_BPp <- pairwise_termsim(CA4_BP)
emapplot(CA4_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("CA4 Biological Process Modules")

dev.off()

# For CA3

pdf(file = here::here("plots","pseudobulked","CA3_GO.pdf"), width = 16, height = 10)

barplot(CA3_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("CA3 Cellular Compartment")
barplot(CA3_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("CA3 Molecular Function")
barplot(CA3_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("CA3 Biological Process")

CA3_CCx <- setReadable(CA3_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(CA3_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("CA3 Cellular Compartment Network")
CA3_MFx <- setReadable(CA3_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(CA3_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("CA3 Molecular Function Network")
CA3_BPx <- setReadable(CA3_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(CA3_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("CA3 Biological Process Network")

CA3_CCp <- pairwise_termsim(CA3_CC)
emapplot(CA3_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("CA3 Cellular Compartment Modules")
CA3_MFp <- pairwise_termsim(CA3_MF)
emapplot(CA3_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("CA3 Molecular Function Modules")
CA3_BPp <- pairwise_termsim(CA3_BP)
emapplot(CA3_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("CA3 Biological Process Modules")

dev.off()

# For Molecular Layer

pdf(file = here::here("plots","pseudobulked","ML_GO.pdf"), width = 16, height = 10)

barplot(ML_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("ML Cellular Compartment")
barplot(ML_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("ML Molecular Function")
barplot(ML_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("ML Biological Process")

ML_CCx <- setReadable(ML_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(ML_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("ML Cellular Compartment Network")
ML_MFx <- setReadable(ML_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(ML_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("ML Molecular Function Network")
ML_BPx <- setReadable(ML_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(ML_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("ML Biological Process Network")

ML_CCp <- pairwise_termsim(ML_CC)
emapplot(ML_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("ML Cellular Compartment Modules")
ML_MFp <- pairwise_termsim(ML_MF)
emapplot(ML_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("ML Molecular Function Modules")
ML_BPp <- pairwise_termsim(ML_BP)
emapplot(ML_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("ML Biological Process Modules")

dev.off()

# For Cluster 8

pdf(file = here::here("plots","pseudobulked","clust_8_GO.pdf"), width = 16, height = 10)

barplot(clust_8_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_8 Cellular Compartment")
barplot(clust_8_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_8 Molecular Function")
barplot(clust_8_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_8 Biological Process")

clust_8_CCx <- setReadable(clust_8_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_8_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_8 Cellular Compartment Network")
clust_8_MFx <- setReadable(clust_8_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_8_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_8 Molecular Function Network")
clust_8_BPx <- setReadable(clust_8_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_8_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_8 Biological Process Network")

clust_8_CCp <- pairwise_termsim(clust_8_CC)
emapplot(clust_8_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_8 Cellular Compartment Modules")
clust_8_MFp <- pairwise_termsim(clust_8_MF)
emapplot(clust_8_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_8 Molecular Function Modules")
clust_8_BPp <- pairwise_termsim(clust_8_BP)
emapplot(clust_8_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_8 Biological Process Modules")

dev.off()

# For comparing clusters

pdf(file = here::here("plots","pseudobulked","all_comparative_GO.pdf"), width = 16, height = 14)

dotplot(comp_CC, showCategory=8, label_format = 50) +
  ggtitle("Top 8 Comparative GO Cellular Compartment")
dotplot(comp_MF, showCategory=8, label_format = 50) +
  ggtitle("Top 8 Comparative GO Molecular Function")
dotplot(comp_BP, showCategory=8, label_format = 50) +
  ggtitle("Top 8 Comparative GO Biological Process")

comp_CC <- pairwise_termsim(comp_CC)
emapplot(comp_CC, showCategory = 8, color = "p.adjust",
         pie = "count", cex_category = 3, label_format = 20) +
  ggtitle("Top 8 Comparative GO Cellular Compartment Modules")
comp_MF <- pairwise_termsim(comp_MF)
emapplot(comp_MF, showCategory = 8, color = "p.adjust",
         pie = "count", cex_category = 3, label_format = 20) +
  ggtitle("Top 8 Comparative GO Molecular Function Modules")
comp_BP <- pairwise_termsim(comp_BP)
emapplot(comp_BP, showCategory = 8, color = "p.adjust",
         pie = "count", cex_category = 3, label_format = 20) +
  ggtitle("Top 8 Comparative GO Biological Process Modules")

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()