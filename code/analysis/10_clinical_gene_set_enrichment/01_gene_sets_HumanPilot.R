## Adapted from
## https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/check_clinical_gene_sets.R
library("readxl")
library("org.Hs.eg.db")
library("janitor")
library("readr")
library("jaffelab")
library("here")
library("sessioninfo")


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

dir_pilot <-
    "/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/"

##################################
## Satterstrom et al, Cell 2020 ##
##################################
asd_exome <- read_excel(file.path(dir_pilot, "/gene_sets/1-s2.0-S0092867419313984-mmc2.xlsx"),
    sheet = 2
)
asd_exome <- as.data.frame(asd_exome)

## get ensembl IDs
asd_exome_geneList <- apply(
    asd_exome[
        ,
        c(
            "ASC33_2014",
            "SSC27_2014",
            "ASC65_2015",
            "ASC102_2018",
            "ASD53",
            "DDID49"
        )
    ], 2,
    function(x) {
        asd_exome$ensembl_gene_id[x == 1]
    }
)
names(asd_exome_geneList) <- gsub("_", ".", names(asd_exome_geneList))
names(asd_exome_geneList) <- paste0(
    "Gene_Satterstrom_",
    names(asd_exome_geneList)
)

###############
### SFARI #####
###############

asd_sfari <- read.csv(
    file.path(
        dir_pilot,
        "gene_sets/SFARI-Gene_genes_01-03-2020release_02-04-2020export.csv"
    ),
    as.is = TRUE
)
asd_sfari_geneList <- list(
    Gene_SFARI_all = asd_sfari$ensembl.id,
    Gene_SFARI_high = asd_sfari$ensembl.id[asd_sfari$gene.score < 3],
    Gene_SFARI_syndromic = asd_sfari$ensembl.id[asd_sfari$syndromic == 1]
)

####################
### birnbaum sets ##
####################

birnbaum <- read_excel(
    file.path(
        dir_pilot,
        "gene_sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx"
    ),
    sheet = 1
)
ens2 <- select(org.Hs.eg.db,
    columns = c("ENSEMBL", "ENTREZID"),
    keys = as.character(unique(birnbaum$`EntrezGene ID`))
)
birnbaum$ensemblID <- ens2$ENSEMBL[match(birnbaum$`EntrezGene ID`, ens2$ENTREZID)]

birnbaum_geneList <- split(birnbaum$ensemblID, birnbaum$`Gene Set`)
names(birnbaum_geneList) <- gsub(" ", ".", names(birnbaum_geneList))
names(birnbaum_geneList) <- gsub("-", ".", names(birnbaum_geneList))
names(birnbaum_geneList) <- paste0(
    "Gene_Birnbaum_",
    names(birnbaum_geneList)
)
birnbaum_geneList <- birnbaum_geneList[rev(seq(along = birnbaum_geneList))]

######################
## psychENCODE DEGs ##
######################

psychENCODE <- as.data.frame(read_excel(
    file.path(dir_pilot, "gene_sets/aat8127_Table_S1.xlsx"),
    sheet = "DGE"
))

pe_geneList <- with(
    psychENCODE,
    list(
        DE_PE_ASD.Up = ensembl_gene_id[ASD.t.value > 0 & ASD.fdr < 0.05],
        DE_PE_ASD.Down = ensembl_gene_id[ASD.t.value < 0 &
            ASD.fdr < 0.05],
        DE_PE_BD.Up = ensembl_gene_id[BD.t.value > 0 &
            BD.fdr < 0.05],
        DE_PE_BD.Down = ensembl_gene_id[BD.t.value < 0 &
            BD.fdr < 0.05],
        DE_PE_SCZ.Up = ensembl_gene_id[SCZ.t.value > 0 & SCZ.fdr < 0.05],
        DE_PE_SCZ.Down = ensembl_gene_id[SCZ.t.value < 0 &
            SCZ.fdr < 0.05]
    )
)

#################
## brainseq  ####
#################

## DLPFC RiboZero
load(
    "/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda"
)

bs2_geneList <- with(
    outGene,
    list(
        DE_BS2_SCZ.Up = ensemblID[logFC > 0 & adj.P.Val < 0.05],
        DE_BS2_SCZ.Down = ensemblID[logFC < 0 & adj.P.Val < 0.05]
    )
)


##############################
### Sestan DS Neuron 2017? ###

ds <- read_excel(file.path(dir_pilot, "gene_sets/1-s2.0-S0896627316000891-mmc4.xlsx"),
    skip = 2
)
ds <- janitor::clean_names(ds)
ds <- as.data.frame(ds)
ens3 <- select(org.Hs.eg.db,
    columns = c("ENSEMBL", "ENTREZID"),
    keys = as.character(unique(ds$geneid))
)
ds$ensemblID <- ens3$ENSEMBL[match(ds$geneid, ens3$ENTREZID)]
ds$fold_difference_log2 <- as.numeric(ds$fold_difference_log2)
ds$p_value <- readr::parse_number(ds$p_value)
ds$qval <- readr::parse_number(ds$qval)

ds_geneList <- list(
    DE_DS_DS.Up = ds$ensemblID[ds$fold_difference_log2 > 0 &
        ds$qval < 0.05],
    DE_DS_DS.Down = ds$ensemblID[ds$fold_difference_log2 < 0 &
        ds$qval < 0.05]
)

#############################
## various TWAS sets ########
#############################

## brainseq 2
load("/dcl01/ajaffe/data/lab/dg_hippo_paper/rdas/tt_objects_gene.Rdata")
tt_dlpfc <- as.data.frame(tt[tt$region == "DLPFC", ])
tt_dlpfc$ensemblID <- ss(tt_dlpfc$geneid, "\\.")

## PE
twas_sczd <- as.data.frame(read_excel(
    file.path(dir_pilot, "gene_sets/aat8127_Table_S4.xlsx"),
    sheet = "SCZ.TWAS"
))
twas_sczd$TWAS.FDR <- p.adjust(twas_sczd$TWAS.P, "fdr")
twas_asd <- as.data.frame(read_excel(
    file.path(dir_pilot, "gene_sets/aat8127_Table_S4.xlsx"),
    sheet = "ASD.TWAS"
))
twas_asd$TWAS.FDR <- p.adjust(twas_asd$TWAS.P, "fdr")
twas_bpdscz <- as.data.frame(read_excel(
    file.path(dir_pilot, "gene_sets/aat8127_Table_S4.xlsx"),
    sheet = "BD.SCZ"
))
twas_bpdscz$TWAS.FDR <- p.adjust(twas_bpdscz$TWAS.P, "fdr")

twas_geneList <- list(
    TWAS_BS2_SCZ.Up = tt_dlpfc$ensemblID[tt_dlpfc$TWAS.Z > 0 &
        tt_dlpfc$TWAS.FDR < 0.05],
    TWAS_BS2_SCZ.Down = tt_dlpfc$ensemblID[tt_dlpfc$TWAS.Z < 0 &
        tt_dlpfc$TWAS.FDR < 0.05],
    TWAS_PE_SCZ.Up = twas_sczd$GeneID[twas_sczd$TWAS.Z > 0 &
        twas_sczd$TWAS.FDR < 0.05],
    TWAS_PE_SCZ.Down = twas_sczd$GeneID[twas_sczd$TWAS.Z < 0 &
        twas_sczd$TWAS.FDR < 0.05],
    TWAS_PE_ASD.Up = twas_asd$GeneID[twas_asd$TWAS.Z > 0 &
        twas_asd$TWAS.FDR < 0.05],
    TWAS_PE_ASD.Down = twas_asd$GeneID[twas_asd$TWAS.Z < 0 &
        twas_asd$TWAS.FDR < 0.05],
    TWAS_PE_SCZBD.Up = twas_bpdscz$ID[twas_bpdscz$TWAS.Z > 0 &
        twas_bpdscz$TWAS.FDR < 0.05],
    TWAS_PE_SCZBD.Down = twas_bpdscz$ID[twas_bpdscz$TWAS.Z < 0 &
        twas_bpdscz$TWAS.FDR < 0.05]
)

###############
### combine ###
###############

## gene list ##
geneList <- c(
    birnbaum_geneList,
    asd_sfari_geneList,
    asd_exome_geneList,
    pe_geneList,
    bs2_geneList,
    ds_geneList,
    twas_geneList
)


## Save for later
save(geneList, file = here(dir_rdata, "gene_sets_HumanPilot.Rdata"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.2 Patched (2022-11-23 r83388)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2022-11-29
#  pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date (UTC) lib source
#  AnnotationDbi        * 1.60.0   2022-11-01 [2] Bioconductor
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.2.1)
#  Biobase              * 2.58.0   2022-11-01 [1] Bioconductor
#  BiocGenerics         * 0.44.0   2022-11-01 [2] Bioconductor
#  Biostrings             2.66.0   2022-11-01 [2] Bioconductor
#  bit                    4.0.5    2022-11-15 [2] CRAN (R 4.2.2)
#  bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.2.1)
#  bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.2.1)
#  blob                   1.2.3    2022-04-10 [2] CRAN (R 4.2.1)
#  cachem                 1.0.6    2021-08-19 [2] CRAN (R 4.2.1)
#  cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.2.1)
#  cli                    3.4.1    2022-09-23 [2] CRAN (R 4.2.1)
#  colorout               1.2-2    2022-11-02 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.2.1)
#  crayon                 1.5.2    2022-09-29 [2] CRAN (R 4.2.1)
#  DBI                    1.1.3    2022-06-18 [2] CRAN (R 4.2.1)
#  DelayedArray           0.24.0   2022-11-01 [2] Bioconductor
#  digest                 0.6.30   2022-10-18 [2] CRAN (R 4.2.1)
#  dplyr                  1.0.10   2022-09-01 [2] CRAN (R 4.2.1)
#  ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.2.1)
#  fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.2.1)
#  fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.2.1)
#  fs                     1.5.2    2021-12-08 [2] CRAN (R 4.2.1)
#  gargle                 1.2.1    2022-09-08 [2] CRAN (R 4.2.1)
#  generics               0.1.3    2022-07-05 [2] CRAN (R 4.2.1)
#  GenomeInfoDb           1.34.3   2022-11-10 [1] Bioconductor
#  GenomeInfoDbData       1.2.9    2022-09-29 [2] Bioconductor
#  GenomicRanges          1.50.1   2022-11-06 [2] Bioconductor
#  ggplot2                3.4.0    2022-11-04 [2] CRAN (R 4.2.2)
#  glue                   1.6.2    2022-02-24 [2] CRAN (R 4.2.1)
#  googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.2.1)
#  gtable                 0.3.1    2022-09-01 [2] CRAN (R 4.2.1)
#  here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.2.1)
#  hms                    1.1.2    2022-08-19 [2] CRAN (R 4.2.1)
#  htmltools              0.5.3    2022-07-18 [2] CRAN (R 4.2.1)
#  htmlwidgets            1.5.4    2021-09-08 [2] CRAN (R 4.2.1)
#  httpuv                 1.6.6    2022-09-08 [2] CRAN (R 4.2.1)
#  httr                   1.4.4    2022-08-17 [2] CRAN (R 4.2.1)
#  IRanges              * 2.32.0   2022-11-01 [2] Bioconductor
#  jaffelab             * 0.99.32  2022-11-02 [1] Github (LieberInstitute/jaffelab@7b7afe3)
#  janitor              * 2.1.0    2021-01-05 [1] CRAN (R 4.2.2)
#  jsonlite               1.8.3    2022-10-21 [2] CRAN (R 4.2.2)
#  KEGGREST               1.38.0   2022-11-01 [2] Bioconductor
#  later                  1.3.0    2021-08-18 [2] CRAN (R 4.2.1)
#  lattice                0.20-45  2021-09-22 [3] CRAN (R 4.2.2)
#  lifecycle              1.0.3    2022-10-07 [2] CRAN (R 4.2.1)
#  limma                  3.54.0   2022-11-01 [1] Bioconductor
#  lubridate              1.9.0    2022-11-06 [2] CRAN (R 4.2.2)
#  magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.2.1)
#  MASS                   7.3-58.1 2022-08-03 [3] CRAN (R 4.2.2)
#  Matrix                 1.5-3    2022-11-11 [2] CRAN (R 4.2.2)
#  MatrixGenerics         1.10.0   2022-11-01 [1] Bioconductor
#  matrixStats            0.63.0   2022-11-18 [2] CRAN (R 4.2.2)
#  memoise                2.0.1    2021-11-26 [2] CRAN (R 4.2.1)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.2.1)
#  nlme                   3.1-160  2022-10-10 [2] CRAN (R 4.2.1)
#  org.Hs.eg.db         * 3.16.0   2022-09-28 [2] Bioconductor
#  pillar                 1.8.1    2022-08-19 [2] CRAN (R 4.2.1)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.2.1)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 4.2.1)
#  promises               1.2.0.1  2021-02-11 [2] CRAN (R 4.2.1)
#  purrr                  0.3.5    2022-10-06 [2] CRAN (R 4.2.1)
#  R6                     2.5.1    2021-08-19 [2] CRAN (R 4.2.1)
#  rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.2.2)
#  RColorBrewer           1.1-3    2022-04-03 [2] CRAN (R 4.2.1)
#  Rcpp                   1.0.9    2022-07-08 [2] CRAN (R 4.2.1)
#  RCurl                  1.98-1.9 2022-10-03 [2] CRAN (R 4.2.1)
#  readr                * 2.1.3    2022-10-01 [2] CRAN (R 4.2.1)
#  readxl               * 1.4.1    2022-08-17 [2] CRAN (R 4.2.1)
#  rlang                  1.0.6    2022-09-24 [2] CRAN (R 4.2.1)
#  rmote                  0.3.4    2022-11-02 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.3    2022-04-02 [2] CRAN (R 4.2.1)
#  RSQLite                2.2.19   2022-11-24 [2] CRAN (R 4.2.2)
#  S4Vectors            * 0.36.0   2022-11-01 [1] Bioconductor
#  scales                 1.2.1    2022-08-20 [2] CRAN (R 4.2.1)
#  segmented              1.6-1    2022-11-08 [1] CRAN (R 4.2.2)
#  servr                  0.25     2022-11-04 [1] CRAN (R 4.2.2)
#  sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.2.1)
#  snakecase              0.11.0   2019-05-25 [1] CRAN (R 4.2.2)
#  stringi                1.7.8    2022-07-11 [2] CRAN (R 4.2.1)
#  stringr                1.4.1    2022-08-20 [2] CRAN (R 4.2.1)
#  SummarizedExperiment   1.28.0   2022-11-01 [2] Bioconductor
#  tibble                 3.1.8    2022-07-22 [2] CRAN (R 4.2.1)
#  tidyselect             1.2.0    2022-10-10 [2] CRAN (R 4.2.1)
#  timechange             0.1.1    2022-11-04 [2] CRAN (R 4.2.2)
#  tzdb                   0.3.0    2022-03-28 [2] CRAN (R 4.2.1)
#  utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.2.1)
#  vctrs                  0.5.1    2022-11-16 [2] CRAN (R 4.2.2)
#  xfun                   0.35     2022-11-16 [2] CRAN (R 4.2.2)
#  XVector                0.38.0   2022-11-01 [2] Bioconductor
#  zlibbioc               1.44.0   2022-11-01 [1] Bioconductor
#
#  [1] /users/lcollado/R/4.2.x
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/R/4.2.x/lib64/R/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────
