## Adapted from
## https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/check_clinical_gene_sets.R
library('readxl')
library("org.Hs.eg.db")
library('janitor')
library("readr")
library('jaffelab')
library('here')
library('sessioninfo')


## output directory
dir_rdata <- here::here(
    "processed-data",
    "rdata",
    "spe",
    "10_Clinical_Gene_Set_Enrichment"
)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

#########################
## load in gene sets ####
#########################

dir_pilot <- "/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/"

##################################
## Satterstrom et al, Cell 2020 ##
##################################
asd_exome = read_excel(file.path(dir_pilot, "/gene_sets/1-s2.0-S0092867419313984-mmc2.xlsx"),
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

asd_sfari = read.csv(file.path(dir_pilot, "gene_sets/SFARI-Gene_genes_01-03-2020release_02-04-2020export.csv"),
    as.is = TRUE)
asd_sfari_geneList = list(
    Gene_SFARI_all = asd_sfari$ensembl.id,
    Gene_SFARI_high = asd_sfari$ensembl.id[asd_sfari$gene.score < 3],
    Gene_SFARI_syndromic = asd_sfari$ensembl.id[asd_sfari$syndromic == 1]
)

####################
### birnbaum sets ##
####################

birnbaum = read_excel(file.path(dir_pilot, "gene_sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx"),
    sheet = 1)
ens2 = select(org.Hs.eg.db,
    columns = c("ENSEMBL", "ENTREZID"),
    keys = as.character(unique(birnbaum$`EntrezGene ID`)))
birnbaum$ensemblID = ens2$ENSEMBL[match(birnbaum$`EntrezGene ID`, ens2$ENTREZID)]

birnbaum_geneList = split(birnbaum$ensemblID, birnbaum$`Gene Set`)
names(birnbaum_geneList) = gsub(" ", ".", names(birnbaum_geneList))
names(birnbaum_geneList) = gsub("-", ".", names(birnbaum_geneList))
names(birnbaum_geneList) = paste0("Gene_Birnbaum_",
    names(birnbaum_geneList))
birnbaum_geneList = birnbaum_geneList[rev(seq(along=birnbaum_geneList))]

######################
## psychENCODE DEGs ##
######################

psychENCODE = as.data.frame(read_excel(file.path(dir_pilot, "gene_sets/aat8127_Table_S1.xlsx"), sheet = "DGE"))

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

bs2_geneList = with(outGene,
    list(DE_BS2_SCZ.Up = ensemblID[logFC > 0 & adj.P.Val < 0.05],
        DE_BS2_SCZ.Down = ensemblID[logFC < 0 & adj.P.Val < 0.05]))


##############################
### Sestan DS Neuron 2017? ###

ds = read_excel(file.path(dir_pilot, "gene_sets/1-s2.0-S0896627316000891-mmc4.xlsx"),skip=2)
ds = janitor::clean_names(ds)
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
twas_sczd = as.data.frame(read_excel(file.path(dir_pilot, "gene_sets/aat8127_Table_S4.xlsx"), sheet = "SCZ.TWAS"))
twas_sczd$TWAS.FDR = p.adjust(twas_sczd$TWAS.P, "fdr")
twas_asd = as.data.frame(read_excel(file.path(dir_pilot, "gene_sets/aat8127_Table_S4.xlsx"), sheet = "ASD.TWAS"))
twas_asd$TWAS.FDR = p.adjust(twas_asd$TWAS.P, "fdr")
twas_bpdscz = as.data.frame(read_excel(file.path(dir_pilot, "gene_sets/aat8127_Table_S4.xlsx"), sheet = "BD.SCZ"))
twas_bpdscz$TWAS.FDR = p.adjust(twas_bpdscz$TWAS.P, "fdr")

twas_geneList = list(TWAS_BS2_SCZ.Up = tt_dlpfc$ensemblID[tt_dlpfc$TWAS.Z > 0 & tt_dlpfc$TWAS.FDR < 0.05],
			TWAS_BS2_SCZ.Down = tt_dlpfc$ensemblID[tt_dlpfc$TWAS.Z < 0 & tt_dlpfc$TWAS.FDR < 0.05],
			TWAS_PE_SCZ.Up = twas_sczd$GeneID[twas_sczd$TWAS.Z > 0 & twas_sczd$TWAS.FDR < 0.05],
			TWAS_PE_SCZ.Down = twas_sczd$GeneID[twas_sczd$TWAS.Z < 0 & twas_sczd$TWAS.FDR < 0.05],
			TWAS_PE_ASD.Up = twas_asd$GeneID[twas_asd$TWAS.Z > 0 & twas_asd$TWAS.FDR < 0.05],
			TWAS_PE_ASD.Down = twas_asd$GeneID[twas_asd$TWAS.Z < 0 & twas_asd$TWAS.FDR < 0.05],
			TWAS_PE_SCZBD.Up = twas_bpdscz$ID[twas_bpdscz$TWAS.Z > 0 & twas_bpdscz$TWAS.FDR < 0.05],
			TWAS_PE_SCZBD.Down = twas_bpdscz$ID[twas_bpdscz$TWAS.Z < 0 & twas_bpdscz$TWAS.FDR < 0.05])

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
	twas_geneList
)

save(geneList, file = here(dir_rdata, "gene_sets_HumanPilot.Rdata"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
