###
library('readxl')
library('limma')
library('sessioninfo')
library('parallel')
library('jaffelab')
library('janitor')
library('lattice')
library('org.Hs.eg.db')
library('GenomicFeatures')
library('scran')
library('here')
library('RColorBrewer')
library('ggplot2')
library('fields')

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## load sce pseudobulked object
readRDS(file = here::here("processed-data","rdata","spe","pseudo_bulked_spe","spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k9.RDS"))


###################
## load modeling outputs
#load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("parsed_modeling_results_k",k,".Rdata")))
load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("cluster_modeling_results_k",k,".Rdata")))

## Extract the p-values
logFC0_contrasts <- sapply(eb0_list, function(x) {
  x$coef[, 2, drop = FALSE]
})
rownames(logFC0_contrasts) = rownames(eb_contrasts)
## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
  x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) = rownames(eb_contrasts)
fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
  x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) = rownames(eb_contrasts)

# ## Expand https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/layer_specificity.R#L1445-L1457
do.call(rbind, lapply(seq_len(ncol(fdrs0_contrasts)), function(i) {
  data.frame(
    Layer = colnames(fdrs0_contrasts)[i],
    FDR5_anyT = sum(fdrs0_contrasts[, i] < 0.05),
    FDR5_positiveT = sum(t0_contrasts[, i] > 0 & fdrs0_contrasts[, i] < 0.05),
    FDR10_positiveT = sum(t0_contrasts[, i] > 0 & fdrs0_contrasts[, i] < 0.1)
  )
}))

# Layer FDR5_anyT FDR5_positiveT FDR10_positiveT
# 1     1     10906            612             658
# 2     2      2770           1010            1203
# 3     3      1485           1069            1551
# 4     4      1481           1228            1921
# 5     5      1265            956            1632
# 6     6      5192           1638            1858
# 7     7       694            642             966
# 8     8      2326           1856            2675
# 9     9       700            609             777


## Total genes: 22331


#########################
## load in gene sets ####
#########################

##################################
## Satterstrom et al, Cell 2020 ##    
##################################
asd_exome = read_excel("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/gene_sets/1-s2.0-S0092867419313984-mmc2.xlsx", 
                       sheet = 2)
asd_exome = as.data.frame(asd_exome)

## get ensembl IDs
asd_exome_geneList = apply(asd_exome[,
                                     c("ASC33_2014",
                                       "SSC27_2014",
                                       "ASC65_2015",
                                       "ASC102_2018",
                                       "ASD53",
                                       "DDID49")], 2,
                           function(x)
                             asd_exome$ensembl_gene_id[x == 1])
names(asd_exome_geneList) = gsub("_", ".", names(asd_exome_geneList))
names(asd_exome_geneList) = paste0("Gene_Satterstrom_",
                                   names(asd_exome_geneList))

###############
### SFARI #####
###############

asd_sfari = read.csv("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/gene_sets/SFARI-Gene_genes_01-03-2020release_02-04-2020export.csv",
                     as.is = TRUE)
asd_sfari_geneList = list(
  Gene_SFARI_all = asd_sfari$ensembl.id,
  Gene_SFARI_high = asd_sfari$ensembl.id[asd_sfari$gene.score < 3],
  Gene_SFARI_syndromic = asd_sfari$ensembl.id[asd_sfari$syndromic == 1]
)

# #################
# ## harmonizome ##
# #################

# harmonizome = read.delim(
# "gene_sets/Harmonizome_CTD Gene-Disease Associations Dataset.txt",
# as.is = TRUE,
# skip = 1
# )
# ## add ensembl
# ens = select(org.Hs.eg.db,
# columns = c("ENSEMBL", "ENTREZID"),
# keys = as.character(unique(harmonizome$GeneID)))
# harmonizome$ensemblID = ens$ENSEMBL[match(harmonizome$GeneID, ens$ENTREZID)]

# ## split by dx
# harmonizome_geneList = split(harmonizome$ensemblID, harmonizome$Disease)

# ## filter by set size
# harmonizome_geneList = harmonizome_geneList[lengths(harmonizome_geneList) >= 100]
# names(harmonizome_geneList) = gsub(" ", ".", names(harmonizome_geneList))
# names(harmonizome_geneList) = paste0("Harmonizome_",
# names(harmonizome_geneList))

####################
### birnbaum sets ##
####################

birnbaum = read_excel("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/gene_sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx",
                      sheet = 1)
ens2 = select(org.Hs.eg.db,
              columns = c("ENSEMBL", "ENTREZID"),
              keys = as.character(unique(birnbaum$`EntrezGene ID`)))
birnbaum$ensemblID = ens2$ENSEMBL[match(birnbaum$`EntrezGene ID`, ens2$ENTREZID)] #need to add ensembl ID to birnbaum data

birnbaum_geneList = split(birnbaum$ensemblID, birnbaum$`Gene Set`)
names(birnbaum_geneList) = gsub(" ", ".", names(birnbaum_geneList))
names(birnbaum_geneList) = gsub("-", ".", names(birnbaum_geneList))
names(birnbaum_geneList) = paste0("Gene_Birnbaum_",
                                  names(birnbaum_geneList))
birnbaum_geneList = birnbaum_geneList[rev(seq(along=birnbaum_geneList))]

######################
## psychENCODE DEGs ##
######################

psychENCODE = as.data.frame(read_excel("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/gene_sets/aat8127_Table_S1.xlsx", sheet = "DGE"))

pe_geneList = with(
  psychENCODE,
  list(
    DE_PE_ASD.Up = ensembl_gene_id[ASD.t.value > 0 & ASD.fdr < 0.05],
    DE_PE_ASD.Down = ensembl_gene_id[ASD.t.value < 0 & ASD.fdr < 0.05],
    DE_PE_BD.Up = ensembl_gene_id[BD.t.value > 0 & BD.fdr < 0.05],
    DE_PE_BD.Down = ensembl_gene_id[BD.t.value < 0 & BD.fdr < 0.05],
    DE_PE_SCZ.Up = ensembl_gene_id[SCZ.t.value > 0 & SCZ.fdr < 0.05],
    DE_PE_SCZ.Down = ensembl_gene_id[SCZ.t.value < 0 & SCZ.fdr < 0.05]
  )
)

#################
## brainseq  ####  
#################

## DLPFC RiboZero
load(
  "/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda"
)

#Error in readChar(con, 5L, useBytes = TRUE) : cannot open the connection
bs2_geneList = with(outGene,
                    list(DE_BS2_SCZ.Up = ensemblID[logFC > 0 & adj.P.Val < 0.05],
                         DE_BS2_SCZ.Down = ensemblID[logFC < 0 & adj.P.Val < 0.05]))


##############################
### Sestan DS Neuron 2017? ###

ds = read_excel("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/gene_sets/1-s2.0-S0896627316000891-mmc4.xlsx",skip=2)
ds = clean_names(ds)
ds = as.data.frame(ds)
ens3 = select(org.Hs.eg.db,
              columns = c("ENSEMBL", "ENTREZID"),
              keys = as.character(unique(ds$geneid)))
ds$ensemblID = ens3$ENSEMBL[match(ds$geneid, ens3$ENTREZID)]
ds$fold_difference_log2 = as.numeric(ds$fold_difference_log2)
ds$p_value = readr::parse_number(ds$p_value)
ds$qval = readr::parse_number(ds$qval)

ds_geneList = list(DE_DS_DS.Up = ds$ensemblID[ds$fold_difference_log2 > 0 & ds$qval < 0.05],
                   DE_DS_DS.Down = ds$ensemblID[ds$fold_difference_log2 < 0 & ds$qval < 0.05])

#############################
## various TWAS sets ########
#############################

## brainseq 2
load("/dcl01/ajaffe/data/lab/dg_hippo_paper/rdas/tt_objects_gene.Rdata")
tt_dlpfc=  as.data.frame(tt[tt$region == "DLPFC",])
tt_dlpfc$ensemblID = ss(tt_dlpfc$geneid, "\\.")

## PE 
twas_sczd = as.data.frame(read_excel("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/gene_sets/aat8127_Table_S4.xlsx", sheet = "SCZ.TWAS"))
twas_sczd$TWAS.FDR = p.adjust(twas_sczd$TWAS.P, "fdr")
twas_asd = as.data.frame(read_excel("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/gene_sets/aat8127_Table_S4.xlsx", sheet = "ASD.TWAS"))
twas_asd$TWAS.FDR = p.adjust(twas_asd$TWAS.P, "fdr")
twas_bpdscz = as.data.frame(read_excel("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/gene_sets/aat8127_Table_S4.xlsx", sheet = "BD.SCZ"))
twas_bpdscz$TWAS.FDR = p.adjust(twas_bpdscz$TWAS.P, "fdr")

twas_geneList = list(TWAS_BS2_SCZ.Up = tt_dlpfc$ensemblID[tt_dlpfc$TWAS.Z > 0 & tt_dlpfc$TWAS.FDR < 0.05],
                     TWAS_BS2_SCZ.Down = tt_dlpfc$ensemblID[tt_dlpfc$TWAS.Z < 0 & tt_dlpfc$TWAS.FDR < 0.05],
                     TWAS_PE_SCZ.Up = twas_sczd$GeneID[twas_sczd$TWAS.Z > 0 & twas_sczd$TWAS.FDR < 0.05],
                     TWAS_PE_SCZ.Down = twas_sczd$GeneID[twas_sczd$TWAS.Z < 0 & twas_sczd$TWAS.FDR < 0.05],
                     TWAS_PE_ASD.Up = twas_asd$GeneID[twas_asd$TWAS.Z > 0 & twas_asd$TWAS.FDR < 0.05],
                     TWAS_PE_ASD.Down = twas_asd$GeneID[twas_asd$TWAS.Z < 0 & twas_asd$TWAS.FDR < 0.05],
                     TWAS_PE_SCZBD.Up = twas_bpdscz$ID[twas_bpdscz$TWAS.Z > 0 & twas_bpdscz$TWAS.FDR < 0.05],
                     TWAS_PE_SCZBD.Down = twas_bpdscz$ID[twas_bpdscz$TWAS.Z < 0 & twas_bpdscz$TWAS.FDR < 0.05])


##############
#### Nagy sn_rna_seq in dlpfc in MDD
####https://www.nature.com/articles/s41593-020-0621-y 
####sup table 6 is marker genes
####sup table 32 is DEGs
#### file = /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/10_clinical_gene_set_enrichment/41593_2020_621_MOESM3_ESM.xlsx
#############
mdd = as.data.frame(read_excel("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/10_clinical_gene_set_enrichment/41593_2020_621_MOESM3_ESM.xlsx", sheet = "Supplementary Table 32", skip = 2, ))

# need to somehow get gene_id or ensemblID
ens4 = select(org.Hs.eg.db,
              columns = c("ENSEMBL", "ENTREZID","SYMBOL"),
              keytypes = "SYMBOL",
              keys = as.character(unique(mdd$Gene)))

mdd_geneList = with(
  mdd,
  list(
    # DE_PE_ASD.Up = ensembl_gene_id[ASD.t.value > 0 & ASD.fdr < 0.05],
    # DE_PE_ASD.Down = ensembl_gene_id[ASD.t.value < 0 & ASD.fdr < 0.05],
    # DE_PE_BD.Up = ensembl_gene_id[BD.t.value > 0 & BD.fdr < 0.05],
    # DE_PE_BD.Down = ensembl_gene_id[BD.t.value < 0 & BD.fdr < 0.05],
    # DE_PE_SCZ.Up = ensembl_gene_id[SCZ.t.value > 0 & SCZ.fdr < 0.05],
    # DE_PE_SCZ.Down = ensembl_gene_id[SCZ.t.value < 0 & SCZ.fdr < 0.05]
  )
)

##############
### snRNAseq ASD
### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7678724/#SD5
### data S4 is DEGS
### data S3 is marker genes
### file = /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/10_clinical_gene_set_enrichment/NIHMS1053005-supplement-Data_S4.xls
#############

asd_rnaseq = as.data.frame(read_excel("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/10_clinical_gene_set_enrichment/NIHMS1053005-supplement-Data_S4.xls", sheet = "ASD_DEGs" ))
asd_rnaseq = clean_names(asd_rnaseq)
asdRNA_geneList = list(DE_ASD_RNA.Up = asd_rnaseq$gene_id[asd_rnaseq$fold_change > 0 & asd_rnaseq$q_value < 0.05],
                   DE_ASD_RNA.Down = asd_rnaseq$gene_id[asd_rnaseq$fold_chang < 0 & asd_rnaseq$q_value < 0.05])


############
###https://www.medrxiv.org/content/10.1101/2020.11.06.20225342v1.full
### snRNAseq SCZ
### supplementary table 6 is DEGS, can't figure out how to download
##########

##########
###snRNAseq and spatial SCZ
###https://www.biorxiv.org/content/10.1101/2020.11.17.386458v2
### supplementary table 2 is marker genes
### supplementary table 4 is DEGs, can't figure out how to download
###########



###############
### combine ###
###############

## gene list ##
geneList = c(
  birnbaum_geneList,
  asd_sfari_geneList,
  asd_exome_geneList,
  pe_geneList,
  bs2_geneList,
  ds_geneList,
  twas_geneList,
  asdRNA_geneList
)

## filter for those present in spatial data
geneList_present = lapply(geneList, function(x) {
  x = x[!is.na(x)]
  x[x %in% rownames(t0_contrasts)]
})

## do enrichment
enrich_stat_list = eb0_list
for (i in seq(along = eb0_list)) {
  layer = t0_contrasts[, i] > 0 & fdrs0_contrasts[, i] < 0.1
  tabList = mclapply(geneList_present, function(g) {
    tt = table(Set = factor(names(layer) %in% g, c(FALSE, TRUE)),
               Layer = factor(layer, c(FALSE, TRUE)))
  }, mc.cores = 8)
  enrichList = lapply(tabList,fisher.test)
  
  o = data.frame(
    OR = sapply(enrichList, "[[", "estimate"),
    Pval = sapply(enrichList, "[[", "p.value"),
    NumSig = sapply(tabList, function(x) x[2,2])
  )
  rownames(o) = gsub(".odds ratio", "", rownames(o))
  enrich_stat_list[[i]] = o
}
enrichTab = do.call("cbind", enrich_stat_list)

#  name
enrichTab$Type = ss(rownames(enrichTab), "_", 1)
enrichTab$Type[enrichTab$Group == "Birnbaum"] = "Birnbaum"
enrichTab$Type[enrichTab$Type == "Gene"] = "ASD"
enrichTab$Group = ss(rownames(enrichTab), "_", 2)
enrichTab$Set = ss(rownames(enrichTab), "_", 3)
enrichTab$ID = rownames(enrichTab)
enrichTab$SetSize = sapply(geneList_present, length)

### save a copy as a supp table
#enrichTabOut = enrichTab[,c(25, 22:24,26, 1:21)]
enrichTabOut = enrichTab
write.csv(enrichTabOut, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/10_clinical_gene_set_enrichment/SupplementaryTableXX_clinical_enrichment.csv",row.names=FALSE)

## look at enrichment
pMat = enrichTab[, grep("Pval", colnames(enrichTab))]
orMat = enrichTab[, grep("OR", colnames(enrichTab))]
colnames(pMat) = ss(colnames(pMat), "\\.")
colnames(orMat) = ss(colnames(orMat), "\\.")
pMat < 0.05 / nrow(pMat)
pMat < 0.001
round(-log10(pMat),1)

# 1     2     3     4     5     6     7     8
# Gene_Birnbaum_SCZ.SNV           FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_SCZ.PGC.GWAS      FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_SCZ.Meta.analysis FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_SCZ.CNV           FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_Neurodegenerative FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_NDD               FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_ID                FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_BPAD.GWAS         FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_ASD.DATABASE      FALSE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE
# Gene_Birnbaum_ASD.CNV           FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_SFARI_all                   TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE
# Gene_SFARI_high                 FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
# Gene_SFARI_syndromic            FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
# DE_PE_ASD.Up                     TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE
# DE_PE_ASD.Down                   TRUE  TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE
# DE_PE_BD.Up                     FALSE FALSE  TRUE FALSE  TRUE  TRUE FALSE FALSE
# DE_PE_BD.Down                    TRUE FALSE  TRUE FALSE  TRUE  TRUE FALSE FALSE
# DE_PE_SCZ.Up                     TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE
# DE_PE_SCZ.Down                  FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE
# DE_DS_DS.Up                     FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE
# DE_DS_DS.Down                   FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE
# TWAS_BS2_SCZ.Up                 FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_BS2_SCZ.Down               FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE
# TWAS_PE_SCZ.Up                  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_SCZ.Down                FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_ASD.Up                  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_ASD.Down                FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_SCZBD.Up                FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_SCZBD.Down              FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# DE_ASD_RNA.Up                   FALSE FALSE  TRUE  TRUE  TRUE FALSE  TRUE FALSE
# DE_ASD_RNA.Down                 FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE
# 9
# Gene_Birnbaum_SCZ.SNV           FALSE
# Gene_Birnbaum_SCZ.PGC.GWAS      FALSE
# Gene_Birnbaum_SCZ.Meta.analysis FALSE
# Gene_Birnbaum_SCZ.CNV           FALSE
# Gene_Birnbaum_Neurodegenerative FALSE
# Gene_Birnbaum_NDD               FALSE
# Gene_Birnbaum_ID                FALSE
# Gene_Birnbaum_BPAD.GWAS         FALSE
# Gene_Birnbaum_ASD.DATABASE      FALSE
# Gene_Birnbaum_ASD.CNV           FALSE
# Gene_SFARI_all                  FALSE
# Gene_SFARI_high                 FALSE
# Gene_SFARI_syndromic            FALSE
# DE_PE_ASD.Up                    FALSE
# DE_PE_ASD.Down                   TRUE
# DE_PE_BD.Up                      TRUE
# DE_PE_BD.Down                    TRUE
# DE_PE_SCZ.Up                     TRUE
# DE_PE_SCZ.Down                   TRUE
# DE_DS_DS.Up                     FALSE
# DE_DS_DS.Down                   FALSE
# TWAS_BS2_SCZ.Up                 FALSE
# TWAS_BS2_SCZ.Down               FALSE
# TWAS_PE_SCZ.Up                  FALSE
# TWAS_PE_SCZ.Down                FALSE
# TWAS_PE_ASD.Up                  FALSE
# TWAS_PE_ASD.Down                FALSE
# TWAS_PE_SCZBD.Up                FALSE
# TWAS_PE_SCZBD.Down              FALSE
# DE_ASD_RNA.Up                    TRUE
# DE_ASD_RNA.Down                 FALSE
# > pMat < 0.001
# 1     2     3     4     5     6     7     8
# Gene_Birnbaum_SCZ.SNV           FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_SCZ.PGC.GWAS      FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_SCZ.Meta.analysis FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_SCZ.CNV           FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_Neurodegenerative FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_NDD               FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_ID                FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_BPAD.GWAS         FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_Birnbaum_ASD.DATABASE      FALSE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE
# Gene_Birnbaum_ASD.CNV           FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Gene_SFARI_all                   TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE
# Gene_SFARI_high                 FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
# Gene_SFARI_syndromic            FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
# DE_PE_ASD.Up                     TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE
# DE_PE_ASD.Down                   TRUE  TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE
# DE_PE_BD.Up                     FALSE FALSE  TRUE FALSE  TRUE  TRUE FALSE FALSE
# DE_PE_BD.Down                    TRUE FALSE  TRUE FALSE  TRUE  TRUE FALSE FALSE
# DE_PE_SCZ.Up                     TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE
# DE_PE_SCZ.Down                  FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE
# DE_DS_DS.Up                     FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE
# DE_DS_DS.Down                   FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_BS2_SCZ.Up                 FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_BS2_SCZ.Down               FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_SCZ.Up                  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_SCZ.Down                FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_ASD.Up                  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_ASD.Down                FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_SCZBD.Up                FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# TWAS_PE_SCZBD.Down              FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# DE_ASD_RNA.Up                   FALSE FALSE  TRUE  TRUE  TRUE FALSE  TRUE FALSE
# DE_ASD_RNA.Down                 FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE
# 9
# Gene_Birnbaum_SCZ.SNV           FALSE
# Gene_Birnbaum_SCZ.PGC.GWAS      FALSE
# Gene_Birnbaum_SCZ.Meta.analysis FALSE
# Gene_Birnbaum_SCZ.CNV           FALSE
# Gene_Birnbaum_Neurodegenerative FALSE
# Gene_Birnbaum_NDD               FALSE
# Gene_Birnbaum_ID                FALSE
# Gene_Birnbaum_BPAD.GWAS         FALSE
# Gene_Birnbaum_ASD.DATABASE      FALSE
# Gene_Birnbaum_ASD.CNV           FALSE
# Gene_SFARI_all                  FALSE
# Gene_SFARI_high                 FALSE
# Gene_SFARI_syndromic            FALSE
# DE_PE_ASD.Up                    FALSE
# DE_PE_ASD.Down                  FALSE
# DE_PE_BD.Up                      TRUE
# DE_PE_BD.Down                    TRUE
# DE_PE_SCZ.Up                     TRUE
# DE_PE_SCZ.Down                   TRUE
# DE_DS_DS.Up                     FALSE
# DE_DS_DS.Down                   FALSE
# TWAS_BS2_SCZ.Up                 FALSE
# TWAS_BS2_SCZ.Down               FALSE
# TWAS_PE_SCZ.Up                  FALSE
# TWAS_PE_SCZ.Down                FALSE
# TWAS_PE_ASD.Up                  FALSE
# TWAS_PE_ASD.Down                FALSE
# TWAS_PE_SCZBD.Up                FALSE
# TWAS_PE_SCZBD.Down              FALSE
# DE_ASD_RNA.Up                    TRUE
# DE_ASD_RNA.Down                 FALSE
# > round(-log10(pMat),1)
# 1    2    3    4    5    6    7    8    9
# Gene_Birnbaum_SCZ.SNV            0.1  0.3  0.1  0.3  0.8  0.9  0.6  0.3  0.4
# Gene_Birnbaum_SCZ.PGC.GWAS       0.0  0.0  0.2  2.1  0.1  0.6  0.5  0.2  0.3
# Gene_Birnbaum_SCZ.Meta.analysis  0.2  0.6  2.3  0.1  2.7  0.2  2.1  0.5  0.2
# Gene_Birnbaum_SCZ.CNV            0.4  1.5  0.1  0.8  0.4  0.5  0.0  0.1  1.0
# Gene_Birnbaum_Neurodegenerative  0.0  0.0  0.2  1.2  0.1  0.2  1.3  1.0  0.6
# Gene_Birnbaum_NDD                1.2  0.5  0.8  0.5  0.4  0.4  0.3  0.8  0.0
# Gene_Birnbaum_ID                 0.1  0.0  1.6  1.5  0.4  0.6  1.2  0.2  0.2
# Gene_Birnbaum_BPAD.GWAS          0.0  0.1  0.8  0.1  0.9  0.9  0.5  0.1  1.7
# Gene_Birnbaum_ASD.DATABASE       0.3  0.0  6.7  1.6  4.9  0.2  4.8  0.3  0.6
# Gene_Birnbaum_ASD.CNV            1.0  1.4  1.0  1.0  0.9  0.8  0.7  0.3  0.9
# Gene_SFARI_all                   3.7  0.5 11.3  1.0  5.8  0.1  8.9  0.9  1.4
# Gene_SFARI_high                  1.4  1.0  4.8  0.3  1.5  2.6  4.5  1.5  0.6
# Gene_SFARI_syndromic             2.3  0.0  1.5  0.6  4.0  0.6  1.5  0.4  0.0
# DE_PE_ASD.Up                    20.9 51.9  3.7  8.8 11.1  0.5  5.6 22.0  0.3
# DE_PE_ASD.Down                   3.5 18.6  2.1 53.7 37.0  0.3 25.6 67.4  2.8
# DE_PE_BD.Up                      0.8  0.6 13.0  0.1  5.8  6.1  0.1  1.7  4.7
# DE_PE_BD.Down                    9.6  0.4  3.7  0.0  4.1  4.5  0.1  2.1  5.1
# DE_PE_SCZ.Up                    16.3 45.5  8.1  8.5  1.0 15.3  6.4  4.7  8.6
# DE_PE_SCZ.Down                   0.0 12.8  2.6  1.9  1.7 30.1  9.0  0.1 37.6
# DE_DS_DS.Up                      2.7  9.8  1.3  1.1  3.3  2.0  0.1  2.1  0.5
# DE_DS_DS.Down                    2.5  3.9  0.0  2.1  0.8  0.6  2.8  1.0  0.2
# TWAS_BS2_SCZ.Up                  0.5  0.0  0.2  0.3  0.0  0.5  0.3  0.1  0.0
# TWAS_BS2_SCZ.Down                2.3  2.7  0.1  0.7  2.4  1.4  1.0  3.0  0.5
# TWAS_PE_SCZ.Up                   0.4  0.3  0.2  0.8  0.3  1.5  0.0  0.2  0.6
# TWAS_PE_SCZ.Down                 2.1  2.1  0.8  0.6  1.3  1.9  1.2  0.8  1.3
# TWAS_PE_ASD.Up                   0.0  0.0  0.9  0.2  0.3  0.0  0.5  0.0  0.0
# TWAS_PE_ASD.Down                 0.0  0.0  0.3  1.3  0.2  0.0  0.4  0.9  0.5
# TWAS_PE_SCZBD.Up                 0.1  0.1  0.3  0.8  1.1  0.3  1.2  1.2  0.8
# TWAS_PE_SCZBD.Down               0.6  0.1  0.2  0.0  0.6  0.6  0.3  0.2  0.3
# DE_ASD_RNA.Up                    0.1  1.6 12.0  3.5  7.7  0.9 17.8  0.9  3.0
# DE_ASD_RNA.Down                  0.6  0.9  1.5  3.4  2.9  1.4  2.5  1.8  1.4

# #######################
# # FDR < 0.05 version ##
# #######################

# # do enrichment
# enrich_stat_list_05 = eb0_list
# for (i in seq(along = eb0_list)) {
# layer = t0_contrasts[, i] > 0 & fdrs0_contrasts[, i] < 0.05
# tabList = mclapply(geneList_present, function(g) {
# tt = table(Set = factor(names(layer) %in% g, c(FALSE, TRUE)),
# Layer = factor(layer, c(FALSE, TRUE)))
# }, mc.cores = 8)
# enrichList = lapply(tabList,fisher.test)

# o = data.frame(
# OR = sapply(enrichList, "[[", "estimate"),
# Pval = sapply(enrichList, "[[", "p.value"),
# NumSig = sapply(tabList, function(x) x[2,2])
# )
# rownames(o) = gsub(".odds ratio", "", rownames(o))
# enrich_stat_list_05[[i]] = o
# }
# enrichTab_05 = do.call("cbind", enrich_stat_list_05)

# name
# enrichTab_05$Type = ss(rownames(enrichTab_05), "_", 1)
# enrichTab_05$Type[enrichTab_05$Group == "Birnbaum"] = "Birnbaum"
# enrichTab_05$Type[enrichTab_05$Type == "Gene"] = "ASD"
# enrichTab_05$Group = ss(rownames(enrichTab_05), "_", 2)
# enrichTab_05$Set = ss(rownames(enrichTab_05), "_", 3)
# enrichTab_05$ID = rownames(enrichTab_05)
# enrichTab_05$SetSize = sapply(geneList_present, length)


######################
## pull out results ##
######################

## summary stats from genes
enrichTab["Gene_SFARI_all",]
enrichTab["Gene_Satterstrom_ASC102.2018",]
enrichTab["Gene_Satterstrom_ASD53",]
enrichTab["Gene_Satterstrom_DDID49",]

## Satterstrom deep dive
# not applicable right now
sat_102_l2= which(t0_contrasts[,"Layer2"] > 0 & fdrs0_contrasts[,"Layer2"] < 0.1 & 
                    rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASC102.2018)
rowData(sce_layer)$gene_name[sat_102_l2]
sat_102_l5= which(t0_contrasts[,"Layer5"] > 0 & fdrs0_contrasts[,"Layer5"] < 0.1 & 
                    rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASC102.2018)
rowData(sce_layer)$gene_name[sat_102_l5]

sat_49_l2= which(t0_contrasts[,"Layer2"] > 0 & fdrs0_contrasts[,"Layer2"] < 0.1 & 
                   rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_DDID49)
cat(rowData(sce_layer)$gene_name[sat_49_l2], sep=", ")

sat_53_l5= which(t0_contrasts[,"Layer5"] > 0 & fdrs0_contrasts[,"Layer5"] < 0.1 & 
                   rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASD53)
cat(rowData(sce_layer)$gene_name[sat_53_l5], sep=", ")

## case control - asd
enrichTab["DE_PE_ASD.Up",]
enrichTab["DE_PE_ASD.Down",]

## case control - sczd
enrichTab[c("DE_PE_SCZ.Up","DE_BS2_SCZ.Up"),]
enrichTab[c("DE_PE_SCZ.Down","DE_BS2_SCZ.Down"),]


################
## make plots ##
################

## make long
enrichLong = reshape2::melt(enrichTab[,c(seq(1,25,by=3),28:32)]) #Using Type, Group, Set, ID as id variables
colnames(enrichLong)[5:6] = c("Layer", "OR")
enrichLong_P = reshape2::melt(enrichTab[,c(seq(2,26,by=3),28:32)])
identical(enrichLong$ID, enrichLong_P$ID)
enrichLong$P = enrichLong_P$value
enrichLong$Layer = ss(as.character(enrichLong$Layer), "\\.")
enrichLong$ID = factor(enrichLong$ID, levels=rev(rownames(enrichTab)))
enrichLong$Set = factor(enrichLong$Set, levels=unique(rev(enrichTab$Set)))
enrichLong$FDR = p.adjust(enrichLong$P, "fdr")

## what p-value controls FDR?
enrichLongSort = enrichLong[order(enrichLong$P),]
max(enrichLongSort$P[enrichLongSort$FDR < 0.05] )
# [1] 0.01193861

## overall ##
enrichLong$P_thresh = enrichLong$P
enrichLong$P_thresh[enrichLong$P_thresh < 2.2e-16] = 2.2e-16

### ASD focus  
enrichLong_ASD = enrichLong[enrichLong$ID %in% 
                              c("Gene_SFARI_all", "Gene_Satterstrom_ASC102.2018",
                                "Gene_Satterstrom_ASD53", "Gene_Satterstrom_DDID49",
                                "DE_PE_ASD.Down", "DE_PE_ASD.Up",
                                "TWAS_PE_ASD.Up", "TWAS_PE_ASD.Down",
                                "DE_ASD_RNA.Up","DE_ASD_RNA.Down"),] #why doesn't this list include the BirnBaum asd sets?
enrichLong_ASD$ID2 =  as.character(droplevels(enrichLong_ASD$Set))
enrichLong_ASD$ID2[enrichLong_ASD$ID2 == "all"] = "SFARI"
enrichLong_ASD$ID2[enrichLong_ASD$ID2 == "ASC102.2018"] = "ASC102"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "DE_PE_ASD.Up"] = "DE.Up"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "DE_PE_ASD.Down"] = "DE.Down"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "TWAS_PE_ASD.Up"] = "TWAS.Up"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "TWAS_PE_ASD.Down"] = "TWAS.Down"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "DE_ASD_RNA.Up"] = "Velmeshev.Up"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "DE_ASD_RNA.Down"] = "Velmeshev.Down"
enrichLong_ASD$ID2 = factor(enrichLong_ASD$ID2, unique(enrichLong_ASD$ID2))

#need to get these two lines working
enrichLong_ASD$LayerFac = factor(as.character(enrichLong_ASD$Layer), 
                                 c(paste0("cluster", 1:9)))
enrichLong_ASD = enrichLong_ASD[order(enrichLong_ASD$ID2, enrichLong_ASD$LayerFac),]

### custom heatmap
midpoint = function(x) x[-length(x)] + diff(x)/2

customLayerEnrichment = function(enrichTab , groups, xlabs, 
                                 Pthresh = 12, ORcut = 3, enrichOnly = FALSE,
                                 layerHeights = c(0,40,55,75,85,110,120,135,145,155),
                                 mypal = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(50)), ...) {
  
  wide_p = -log10( enrichTab[groups,grep("Pval", colnames(enrichTab))])
  wide_p[wide_p > Pthresh] = Pthresh
  wide_p = t(round(wide_p[,
                          c("1.Pval", "2.Pval", "3.Pval", "4.Pval", "5.Pval","6.Pval", "7.Pval","8.Pval","9.Pval")],2))
  
  wide_or = enrichTab[groups,grep("OR", colnames(enrichTab))]
  wide_or= round(t(wide_or[,
                           c("1.OR", "2.OR", "3.OR", "4.OR", "5.OR", "6.OR", "7.OR","8.OR","9.OR")]),1)
  if(enrichOnly) wide_p[wide_or < 1] = 0
  wide_or[wide_p < ORcut] = ""
  
  image.plot(x = seq(0,ncol(wide_p),by=1), y = layerHeights, z = as.matrix(t(wide_p)),
             col = mypal,xaxt="n", yaxt="n",xlab = "", ylab="", ...)
  axis(2, c(1:9), at = midpoint(layerHeights),las=1)
  axis(1, rep("", ncol(wide_p)), at = seq(0.5,ncol(wide_p)-0.5))
  text(x = seq(0.5,ncol(wide_p)-0.5), y=-1*max(nchar(xlabs))/2, xlabs,
       xpd=TRUE, srt=45,cex=2,adj= 1)
  abline(h=layerHeights,v=0:ncol(wide_p))
  text(x = rep(seq(0.5,ncol(wide_p)-0.5),each = nrow(wide_p)), 
       y = rep(midpoint(layerHeights), ncol(wide_p)),
       as.character(wide_or),cex=1.5,font=2)
}

pdf("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/10_clinical_gene_set_enrichment/asd_geneSet_heatmap.pdf",w=6)
par(mar=c(8,4.5,2.5,1), cex.axis=2,cex.lab=2)
groups = unique(as.character(enrichLong_ASD$ID))[1:10]
xlabs  = as.character(enrichLong_ASD$ID2[match(groups, enrichLong_ASD$ID)])
customLayerEnrichment(enrichTab, groups,xlabs, enrichOnly=TRUE)
#abline(v=4,lwd=3)
text(x = 4, y = 165, c("ASD"), xpd=TRUE,cex=2.5,font=2)

dev.off()


pdf("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/10_clinical_gene_set_enrichment/sczd_geneSet_heatmap.pdf",w=8)
par(mar=c(8,4.5,2.5,1), cex.axis=2,cex.lab=2)
groups =c("DE_PE_SCZ.Up", "DE_PE_SCZ.Down",
          "DE_BS2_SCZ.Up", "DE_BS2_SCZ.Down",
          "TWAS_BS2_SCZ.Up", "TWAS_BS2_SCZ.Down", "TWAS_PE_SCZ.Up",
          "TWAS_PE_SCZ.Down")
# groups =c("DE_PE_SCZ.Up", "DE_PE_SCZ.Down", 
#           "TWAS_BS2_SCZ.Up", "TWAS_BS2_SCZ.Down", "TWAS_PE_SCZ.Up",
#           "TWAS_PE_SCZ.Down")
xlabs = ss(gsub("_SCZ", "", groups), "_", 2)
customLayerEnrichment(enrichTab, groups,xlabs, enrichOnly=TRUE)
abline(v=4,lwd=3)
text(x = c(2,6), y = 162, c("SCZD-DE", "SCZD-TWAS"), xpd=TRUE,cex=2.5,font=2)
dev.off()


pdf("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/10_clinical_gene_set_enrichment/suppXX_birnbaum_geneSet_heatmap.pdf",w=8)
par(mar=c(12,5.5,2.5,1), cex.axis=2,cex.lab=2)
groups =grep(enrichTab$ID, pattern = "Birnbaum", value=TRUE)
xlabs = ss(groups, "_", 3)
customLayerEnrichment(enrichTab, groups,xlabs, enrichOnly=TRUE,
                      breaks = seq(0,12,len = 51))
customLayerEnrichment(enrichTab, groups,xlabs, enrichOnly=TRUE)
dev.off()


#############################
### GSEA #####################
#############################

## do enrichment
gst_tab = apply(t0_contrasts, 2, function(tt) {
  sapply(geneList_present, function(g) {
    geneSetTest(index = rownames(t0_contrasts) %in% g,
                statistics = tt, alternative = "up")
  })
})
round(-log10(gst_tab),1)

## check densities 
mypar(ncol(t0_contrasts),1)
g = rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASC102.2018 
g_asd = rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASD53
g_dd = rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_DDID49
for(i in 1:ncol(t0_contrasts)) {
  plot(density(t0_contrasts[!g,i]), lwd=2,col="black",xlab="",
       main=colnames(t0_contrasts)[i],xlim = c(-10,15))
  lines(density(t0_contrasts[g,i]), lwd=2,col="red")
  lines(density(t0_contrasts[g_asd,i]), lwd=2,col="red",lty=2)
  lines(density(t0_contrasts[g_dd,i]), lwd=2,col="red",lty=3)
  abline(v=0,lty=2)
}

pdf("pdf/ASD_genes_layer_density.pdf",h=4,useDingbats=FALSE)
par(mar=c(5,6,1,1),cex.axis= 1.4,cex.lab=1.8)
for(i in 1:ncol(t0_contrasts)) {
  layer = t0_contrasts[, i] > 0 & fdrs0_contrasts[, i] < 0.1
  plot(density(t0_contrasts[!g,i]), lwd=3,col="black",
       xlab=paste0(colnames(t0_contrasts)[i], ": Specificity T-stats"),
       sub = "", main="",xlim = c(-8,8))
  lines(density(t0_contrasts[g,i]), lwd=3,col="red")
  lines(density(t0_contrasts[g_asd,i]), lwd=3,col="red",lty=2)
  lines(density(t0_contrasts[g_dd,i]), lwd=3,col="red",lty=3)
  abline(v=0,lty=2)
  abline(v=	min(t0_contrasts[layer,i]))
  
  ll = ifelse(i == 1, "topright", "topleft")
  legend(ll, c("BG", "102 All", "53 ASD", "49 DDID"), bty="n",
         col = c("black","red","red","red"),	lty = c(1,1,2,3),cex=1.5,lwd=4)
}
dev.off()

diag(cor(t(-log10(gst_tab)),t(-log10(pMat))))