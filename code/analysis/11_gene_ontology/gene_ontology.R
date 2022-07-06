################################################
# Gene Ontoloyg from https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/code/Pseudobulk/all_gene_ontology_pseudobulk.R
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

load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression","parsed_modeling_results_k9.Rdata"))

#load mito genes to drop
load(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/08_layer_differential_expression/mito_genes.rda")

enriched_model <- modeling_results$enrichment
rownames(enriched_model)<-enriched_model$ensembl
enriched_model <- enriched_model[!(rownames(enriched_model)%in%drop_mt),]

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

list_9 <- enriched_model %>%
  dplyr::select(gene, ensembl, fdr_9, t_stat_9)

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
clust_3 = top_list_3[,2]
## feature 2: named vector
names(clust_3) = as.character(top_list_3[,1])
# final vector
clust_3 <- names(clust_3)
clust_3 = bitr(clust_3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(clust_3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.86% of input gene IDs are fail to map...

top_list_4 <- list_4 %>%
  filter(fdr_4 < 0.05) %>%
  dplyr::arrange(desc(t_stat_4)) %>%
  dplyr::select(gene, t_stat_4)

## feature 1: numeric vector
clust_4 = top_list_4[,2]
## feature 2: named vector
names(clust_4) = as.character(top_list_4[,1])
# final vector
clust_4 <- names(clust_4)
clust_4 = bitr(clust_4, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

top_list_5 <- list_5 %>%
  filter(fdr_5 < 0.05) %>%
  dplyr::arrange(desc(t_stat_5)) %>%
  dplyr::select(gene, t_stat_5)

## feature 1: numeric vector
clust_5 = top_list_5[,2]
## feature 2: named vector
names(clust_5) = as.character(top_list_5[,1])
# final vector
clust_5 <- names(clust_5)
clust_5 = bitr(clust_5, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(clust_5, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   1.79% of input gene IDs are fail to map...

top_list_6 <- list_6 %>%
  filter(fdr_6 < 0.05) %>%
  dplyr::arrange(desc(t_stat_6)) %>%
  dplyr::select(gene, t_stat_6)

## feature 1: numeric vector
clust_6 = top_list_6[,2]
## feature 2: named vector
names(clust_6) = as.character(top_list_6[,1])
# final vector
clust_6 <- names(clust_6)
clust_6 = bitr(clust_6, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(clust_6, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   2.84% of input gene IDs are fail to map...

list_7 <- list_7 %>%
  arrange(desc(t_stat_7))
# BayesSpace cluster 7 (Molec. Layer) has only 1 gene that is FDR<0.05, so do top 20 instead
top_list_7 <- head(list_7, 20) %>%
  dplyr::select(gene, t_stat_7)

## feature 1: numeric vector
clust_7 = top_list_7[,2]
## feature 2: named vector
names(clust_7) = as.character(top_list_7[,1])
# final vector
clust_7 <- names(clust_7)
clust_7 = bitr(clust_7, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(clust_7, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
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

top_list_9 <- list_9 %>%
  filter(fdr_9 < 0.05) %>%
  dplyr::arrange(desc(t_stat_9)) %>%
  dplyr::select(gene, t_stat_9)

## feature 1: numeric vector
clust_9 = top_list_9[,2]
## feature 2: named vector
names(clust_9) = as.character(top_list_9[,1])
# final vector
clust_9 <- names(clust_9)
clust_9 = bitr(clust_9, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")



# Make a list for comparing cluster ontologies
clust_compare <- list(clust_1$ENTREZID, clust_2$ENTREZID, clust_3$ENTREZID, clust_4$ENTREZID,
                      clust_5$ENTREZID, clust_6$ENTREZID, clust_7$ENTREZID, clust_8$ENTREZID, clust_9$ENTREZID)

names(clust_compare) <- c("clust_1", "clust_2", "clust_3", "clust_4", "clust_5", "clust_6", "clust_7", "clust_8","clust_9")

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

clust_3_CC <- enrichGO(gene = clust_3$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_3_MF <- enrichGO(gene = clust_3$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_3_BP <- enrichGO(gene = clust_3$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_4_CC <- enrichGO(gene = clust_4$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_4_MF <- enrichGO(gene = clust_4$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_4_BP <- enrichGO(gene = clust_4$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_5_CC <- enrichGO(gene = clust_5$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_5_MF <- enrichGO(gene = clust_5$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_5_BP <- enrichGO(gene = clust_5$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_6_CC <- enrichGO(gene = clust_6$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_6_MF <- enrichGO(gene = clust_6$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_6_BP <- enrichGO(gene = clust_6$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   readable = TRUE,
)

clust_7_CC <- enrichGO(gene = clust_7$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "fdr",
                  pvalueCutoff = 0.05,
                  qvalueCutoff  = 0.05,
                  readable = TRUE,
)

clust_7_MF <- enrichGO(gene = clust_7$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "fdr",
                  pvalueCutoff = 0.05,
                  qvalueCutoff  = 0.05,
                  readable = TRUE,
)

clust_7_BP <- enrichGO(gene = clust_7$ENTREZID,
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

clust_9_CC <- enrichGO(gene = clust_9$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

clust_9_MF <- enrichGO(gene = clust_9$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE,
)

clust_9_BP <- enrichGO(gene = clust_9$ENTREZID,
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
     clust_3_CC, clust_3_MF, clust_3_BP, clust_4_CC, clust_4_MF, clust_4_BP, clust_5_CC, clust_5_MF, clust_5_BP,
     clust_6_CC, clust_6_MF, clust_6_BP, clust_7_CC, clust_7_MF, clust_7_BP, clust_8_CC, clust_8_MF, clust_8_BP,
     comp_CC, comp_MF, comp_BP,
     file = here::here("processed-data","rdata","spe","11_gene_ontology", "All_samples_enrichedGO.Rdata"))

## Barplots of subontologies

options(ggrepel.max.overlaps = Inf)

# For Cluster 1
pdf(file = here::here("plots","11_gene_ontology","clust_1_GO.pdf"), width = 16, height = 10)

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
pdf(file = here::here("plots","11_gene_ontology","clust_2_GO.pdf"), width = 16, height = 10)

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

# For clust_3

pdf(file = here::here("plots","11_gene_ontology","clust_3_GO.pdf"), width = 16, height = 10)

barplot(clust_3_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_3 Cellular Compartment")
barplot(clust_3_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_3 Molecular Function")
barplot(clust_3_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_3 Biological Process")

clust_3_CCx <- setReadable(clust_3_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_3_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_3 Cellular Compartment Network")
clust_3_MFx <- setReadable(clust_3_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_3_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_3 Molecular Function Network")
clust_3_BPx <- setReadable(clust_3_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_3_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_3 Biological Process Network")

clust_3_CCp <- pairwise_termsim(clust_3_CC)
emapplot(clust_3_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_3 Cellular Compartment Modules")
clust_3_MFp <- pairwise_termsim(clust_3_MF)
emapplot(clust_3_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_3 Molecular Function Modules")
clust_3_BPp <- pairwise_termsim(clust_3_BP)
emapplot(clust_3_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_3 Biological Process Modules")

dev.off()

# For clust_4

pdf(file = here::here("plots","11_gene_ontology","clust_4_GO.pdf"), width = 16, height = 10)

barplot(clust_4_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_4 Cellular Compartment")
barplot(clust_4_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_4 Molecular Function")
barplot(clust_4_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_4 Biological Process")

clust_4_CCx <- setReadable(clust_4_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_4_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_4 Cellular Compartment Network")
clust_4_MFx <- setReadable(clust_4_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_4_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_4 Molecular Function Network")
clust_4_BPx <- setReadable(clust_4_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_4_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_4 Biological Process Network")

clust_4_CCp <- pairwise_termsim(clust_4_CC)
emapplot(clust_4_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_4 Cellular Compartment Modules")
clust_4_MFp <- pairwise_termsim(clust_4_MF)
emapplot(clust_4_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_4 Molecular Function Modules")
clust_4_BPp <- pairwise_termsim(clust_4_BP)
emapplot(clust_4_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_4 Biological Process Modules")

dev.off()

# For clust_5

pdf(file = here::here("plots","11_gene_ontology","clust_5_GO.pdf"), width = 16, height = 10)

barplot(clust_5_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_5 Cellular Compartment")
barplot(clust_5_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_5 Molecular Function")
barplot(clust_5_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_5 Biological Process")

clust_5_CCx <- setReadable(clust_5_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_5_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_5 Cellular Compartment Network")
clust_5_MFx <- setReadable(clust_5_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_5_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_5 Molecular Function Network")
clust_5_BPx <- setReadable(clust_5_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_5_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_5 Biological Process Network")

clust_5_CCp <- pairwise_termsim(clust_5_CC)
emapplot(clust_5_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_5 Cellular Compartment Modules")
clust_5_MFp <- pairwise_termsim(clust_5_MF)
emapplot(clust_5_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_5 Molecular Function Modules")
clust_5_BPp <- pairwise_termsim(clust_5_BP)
emapplot(clust_5_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_5 Biological Process Modules")

dev.off()

# For clust_6

pdf(file = here::here("plots","11_gene_ontology","clust_6_GO.pdf"), width = 16, height = 10)

barplot(clust_6_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_6 Cellular Compartment")
barplot(clust_6_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_6 Molecular Function")
barplot(clust_6_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_6 Biological Process")

clust_6_CCx <- setReadable(clust_6_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_6_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_6 Cellular Compartment Network")
clust_6_MFx <- setReadable(clust_6_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_6_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_6 Molecular Function Network")
clust_6_BPx <- setReadable(clust_6_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_6_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_6 Biological Process Network")

clust_6_CCp <- pairwise_termsim(clust_6_CC)
emapplot(clust_6_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_6 Cellular Compartment Modules")
clust_6_MFp <- pairwise_termsim(clust_6_MF)
emapplot(clust_6_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_6 Molecular Function Modules")
clust_6_BPp <- pairwise_termsim(clust_6_BP)
emapplot(clust_6_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_6 Biological Process Modules")

dev.off()

# For Molecular Layer

pdf(file = here::here("plots","11_gene_ontology","clust_7_GO.pdf"), width = 16, height = 10)

barplot(clust_7_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_7 Cellular Compartment")
barplot(clust_7_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_7 Molecular Function")
barplot(clust_7_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("clust_7 Biological Process")

clust_7_CCx <- setReadable(clust_7_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_7_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_7 Cellular Compartment Network")
clust_7_MFx <- setReadable(clust_7_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_7_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_7 Molecular Function Network")
clust_7_BPx <- setReadable(clust_7_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_7_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("clust_7 Biological Process Network")

clust_7_CCp <- pairwise_termsim(clust_7_CC)
emapplot(clust_7_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_7 Cellular Compartment Modules")
clust_7_MFp <- pairwise_termsim(clust_7_MF)
emapplot(clust_7_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_7 Molecular Function Modules")
clust_7_BPp <- pairwise_termsim(clust_7_BP)
emapplot(clust_7_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("clust_7 Biological Process Modules")

dev.off()

# For Cluster 8

pdf(file = here::here("plots","11_gene_ontology","clust_8_GO.pdf"), width = 16, height = 10)

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

# For Cluster 9

pdf(file = here::here("plots","11_gene_ontology","clust_9_GO.pdf"), width = 16, height = 10)

barplot(clust_9_CC, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_9 Cellular Compartment")
barplot(clust_9_MF, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_9 Molecular Function")
barplot(clust_9_BP, showCategory=20, x= "GeneRatio") +
  ggtitle("cluster_9 Biological Process")

clust_9_CCx <- setReadable(clust_9_CC, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_9_CCx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_9 Cellular Compartment Network")
clust_9_MFx <- setReadable(clust_9_MF, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_9_MFx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_9 Molecular Function Network")
clust_9_BPx <- setReadable(clust_9_BP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(clust_9_BPx, showCategory = 20, colorEdge = TRUE, cex_category = 0.2,
         cex_gene = 0.1, cex_label_category = 0.4, cex_label_gene = 0.2,
         shadow_text = "category", layout = "kk") +
  ggtitle("cluster_9 Biological Process Network")

clust_9_CCp <- pairwise_termsim(clust_9_CC)
emapplot(clust_9_CCp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_9 Cellular Compartment Modules")
clust_9_MFp <- pairwise_termsim(clust_9_MF)
emapplot(clust_9_MFp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_9 Molecular Function Modules")
clust_9_BPp <- pairwise_termsim(clust_9_BP)
emapplot(clust_9_BPp, showCategory = 20, color = "p.adjust",
         cex_label_category = 0.5, layout = "kk") +
  ggtitle("cluster_9 Biological Process Modules")

dev.off()

# For comparing clusters

pdf(file = here::here("plots","11_gene_ontology","all_comparative_GO.pdf"), width = 18, height = 14)

dotplot(comp_CC, showCategory=9, label_format = 50) +
  ggtitle("Top 9 Comparative GO Cellular Compartment")
dotplot(comp_MF, showCategory=9, label_format = 50) +
  ggtitle("Top 9 Comparative GO Molecular Function")
dotplot(comp_BP, showCategory=9, label_format = 50) +
  ggtitle("Top 9 Comparative GO Biological Process")

comp_CC <- pairwise_termsim(comp_CC)
emapplot(comp_CC, showCategory = 9, color = "p.adjust",
         pie = "count", cex_category = 3, label_format = 20) +
  ggtitle("Top 9 Comparative GO Cellular Compartment Modules")
comp_MF <- pairwise_termsim(comp_MF)
emapplot(comp_MF, showCategory = 9, color = "p.adjust",
         pie = "count", cex_category = 3, label_format = 20) +
  ggtitle("Top 9 Comparative GO Molecular Function Modules")
comp_BP <- pairwise_termsim(comp_BP)
emapplot(comp_BP, showCategory = 9, color = "p.adjust",
         pie = "count", cex_category = 3, label_format = 20) +
  ggtitle("Top 9 Comparative GO Biological Process Modules")

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()