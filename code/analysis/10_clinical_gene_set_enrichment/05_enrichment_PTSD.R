library("jaffelab")
library("here")
library("purrr")
library("dplyr")
library("spatialLIBD")
library("sessioninfo")

## Input dir
dir_input <- here::here(
    "processed-data",
    "rdata",
    "spe",
    "07_layer_differential_expression"
)

## Output directories
dir_rdata <- here::here(
    "processed-data",
    "rdata",
    "spe",
    "10_clinical_gene_set_enrichment",
    "PTSD"
)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully

dir_plots <- here::here(
    "plots",
    "10_clinical_gene_set_enrichment",
    "05_enrichment_PTSD"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## Load the gene sets from Singel Nucleus datasets
de_gene <- mget(load(here(dir_rdata, "genes.rda"), verbose = TRUE))
de_protein <- mget(load(here(dir_rdata, "proteins.rda"), verbose = TRUE))

names(de_gene)
names(de_protein)

split_region <- function(list) {
    n <- stringr::str_split(names(list), "_")
    n <- map_chr(n, ~ tail(.x, 1))
    list2 <- map(splitit(n), ~ list[.x])
    return(list2)
}

de_gene2 <- split_region(flatten(de_gene))
names(de_gene2) <- paste0("gene_", names(de_gene2))

de_protein2 <- split_region(flatten(de_protein))
names(de_protein2) <- paste0("protein_", names(de_protein2))

ptsd_genes <- c(de_gene2, de_protein2)
names(ptsd_gene_list)

## convert to list and filter
min_gene <- 10

ptsd_gene_list <- map_depth(ptsd_genes, 2, ~ .$gene)
ptsd_gene_list <- map(ptsd_gene_list, ~ .x[map_int(.x, length) > 20])
ptsd_gene_list <- map_depth(ptsd_gene_list, 2, ~ ss(.x, "\\."))
map(ptsd_gene_list, ~ map_int(.x, length))

# $gene_DLPFC
# MDD_DLPFC      PTSD_MDD_DLPFC      MDD_Male_DLPFC PTSD_MDD_Male_DLPFC
# 315                 414                  34                  22
#
# $gene_mPFC
# MDD_mPFC                 PTSD_mPFC             PTSD_MDD_mPFC      MDD_ChildTrauma_mPFC
# 1264                      1498                      2322                      1070
# PTSD_ChildTrauma_mPFC PTSD_MDD_ChildTrauma_mPFC          PTSD_Female_mPFC      PTSD_MDD_Female_mPFC
# 408                       726                       211                       105
# MDD_Male_mPFC        PTSD_MDD_Male_mPFC          MDD_suicide_mPFC
# 410                       226                       681
#
# $protein_mPFC
# MDD_mPFC                 PTSD_mPFC             PTSD_MDD_mPFC      MDD_ChildTrauma_mPFC
# 56                        46                       100                        78
# PTSD_MDD_ChildTrauma_mPFC           MDD_Female_mPFC          PTSD_Female_mPFC             MDD_Male_mPFC
# 46                        36                        39                        23

## Specify what k's we want to look at
k_list <- c(2, 7, 9, 16, 28)
names(k_list) <- paste0("k", sprintf("%02d", k_list))

## Load the modeling results from the BayesSpace models
bayesSpace_registration_fn <-
    map(k_list, ~ here(
        dir_input,
        paste0(
            "modeling_results_BayesSpace_k",
            sprintf("%02d", .x),
            ".Rdata"
        )
    ))
bayesSpace_registration <-
    lapply(bayesSpace_registration_fn, function(x) {
        get(load(x))
    })

## Read in the spatial registration labels
bayes_anno <-
    read.csv(
        file = here(
            "processed-data",
            "rdata",
            "spe",
            "08_spatial_registration",
            "bayesSpace_layer_annotations.csv"
        )
    ) |>
    select(layer_combo,
        test = cluster,
        Annotation = bayesSpace
    )

## Takes 2-3 min to run
enriched <-
    map(ptsd_gene_list, function(gl) {
        map(
            bayesSpace_registration,
            ~ gene_set_enrichment(
                gene_list = gl,
                modeling_results = .x,
                model_type = "enrichment"
            ) |>
                left_join(bayes_anno, by = "test") |>
                mutate(
                    test = factor(layer_combo, levels = bayes_anno$layer_combo[bayes_anno$layer_combo %in% layer_combo])
                ) |>
                select(-c(layer_combo, Annotation, fdr_cut, model_type))
        )
    })

enriched$gene_DLPFC$k09

## Save for later
save(enriched, file = here(dir_rdata, "enriched_PTSD.Rdata"))

#### Prep for plotting ####
## source updated enrichment plot
source(here("code", "analysis", "10_clinical_gene_set_enrichment", "gene_set_enrichment_plot_complex.R"))

## how many genes for each domain or DE set,
gene_enrichment_count <- map(bayesSpace_registration, function(r) {
    en_count <- get_gene_enrichment_count(r)
    ## reorder to match layer_combo
    rownames(en_count) <- bayes_anno$layer_combo[match(rownames(en_count), bayes_anno$test)]
    layer_order <- bayes_anno$layer_combo[bayes_anno$layer_combo %in% rownames(en_count)]
    return(en_count[layer_order, , drop = FALSE])
})

gene_list_count <- map(ptsd_gene_list, get_gene_list_count)


#### Plot Enrichments ####
pdf(here(dir_plots, "Enrich_PTSD_DLPFC_genes_k09.pdf"), height = 8, width = 10)
gene_set_enrichment_plot_complex(enriched$gene_DLPFC$k09,
    gene_count_col = gene_list_count[["gene_DLPFC"]],
    gene_count_row = gene_enrichment_count[["k09"]],
    anno_title_col = "n DE Genes",
    anno_title_row = "n Domain\nGenes"
)
dev.off()

walk2(enriched, names(enriched), function(enriched, ds_name) {
    pdf(here(dir_plots, paste0("Enrich_PTSD_", ds_name, ".pdf")), height = 8, width = 9)
    map2(enriched, names(enriched), function(x, k) {
        message(ds_name, " - ", k)

        print(gene_set_enrichment_plot_complex(x,
            gene_count_col = gene_list_count[[ds_name]],
            gene_count_row = gene_enrichment_count[[k]],
            anno_title_col = "n DE Genes",
            anno_title_row = "n Domain\nGenes"
        ))
    })
    dev.off()
})

## Select gene set plot @ k09
enriched_select <- rbind(
    enriched$gene_DLPFC$k09 |> filter(ID %in% c("MDD_DLPFC", "PTSD_MDD_DLPFC")),
    enriched$gene_mPFC$k09 |> filter(ID %in% c("MDD_mPFC", "PTSD_MDD_mPFC"))
)

gene_list_count_select <- rbind(
    gene_list_count$gene_DLPFC[c("MDD_DLPFC", "PTSD_MDD_DLPFC"), , drop = FALSE],
    gene_list_count$gene_mPFC[c("MDD_mPFC", "PTSD_MDD_mPFC"), , drop = FALSE]
)

pdf(here(dir_plots, "Enrich_PTSD_select_k09.pdf"), height = 8, width = 5)
gene_set_enrichment_plot_complex(enriched_select,
    gene_count_col = gene_list_count_select,
    gene_count_row = gene_enrichment_count[["k09"]],
    anno_title_col = "n DE Genes",
    anno_title_row = "n Domain\nGenes"
)
dev.off()

## select gene sets for PTSD team
# Sp09 resolution, RNA and Protein for mPFC for the MDD, PTSD, PTSD_MDD

enriched_select_ptsd <- list(
    gene = enriched$gene_mPFC$k09 |> filter(ID %in% c("MDD_mPFC", "PTSD_mPFC", "PTSD_MDD_mPFC")),
    protein = enriched$protein_mPFC$k09 |> filter(ID %in% c("MDD_mPFC", "PTSD_mPFC", "PTSD_MDD_mPFC"))
)

gene_list_count_select_ptsd <- list(
    gene = gene_list_count$gene_mPFC[c("MDD_mPFC", "PTSD_mPFC", "PTSD_MDD_mPFC"), , drop = FALSE],
    protein = gene_list_count$protein_mPFC[c("MDD_mPFC", "PTSD_mPFC", "PTSD_MDD_mPFC"), , drop = FALSE]
)


pdf(here(dir_plots, "Enrich_PTSD_select_mPFC_k09.pdf"), height = 8, width = 5)
map2(enriched_select_ptsd, gene_list_count_select_ptsd
~ gene_set_enrichment_plot_complex(.x,
        gene_count_col = .y,
        gene_count_row = gene_enrichment_count[["k09"]],
        anno_title_col = "n DE Genes",
        anno_title_row = "n Domain\nGenes"
    ))
dev.off()


# independent color scale
# pal <- c(
#   "white",
#          grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(5)
# )
#
# lgd2 = Legend(col_fun = circlize::colorRamp2(seq(0, 15, by = 1.6),
#                                              c("white", RColorBrewer::brewer.pal(9, "YlOrRd"))),
#               title = "-log10(p-val)",
#               direction = "horizontal")
#
# pdf("enrich_legend.pdf", height = 1, width = 2)
# draw(lgd2)
# dev.off()

#### Interesting gene sets ####

colnames(bayesSpace_registration$k09$enrichment)

gene_check_Sp09D01 <- bayesSpace_registration$k09$enrichment |>
    select(ends_with("Sp09D01"), ensembl, gene) |>
    # filter(fdr_Sp09D01 < 0.05) |>
    mutate(MDD_DLPFC = ensembl %in% ptsd_gene_list$gene_DLPFC$MDD_DLPFC) |>
    dplyr::count(fdr_Sp09D01 < 0.05, MDD_DLPFC)


enrichment_check <- function(domain = "Sp09D01",
    gene_list = ptsd_gene_list$gene_DLPFC$MDD_DLPFC,
    enrichment = bayesSpace_registration$k09$enrichment) {
    enrichment2 <- enrichment |>
        select(ends_with(domain), ensembl, gene) |>
        rename_all(~ stringr::str_replace(., paste0("_", domain), "")) |>
        mutate(
            sig = fdr < 0.1 & t_stat > 0,
            deg = ensembl %in% gene_list
        )

    gene_list <- enrichment2 |>
        filter(deg, sig) |>
        pull(gene)

    gene_tab <- addmargins(table(enrichment2$deg, enrichment2$sig))


    return(list(gl = gene_list, tab = gene_tab))
}

enriched$gene_DLPFC$k09
#           OR         Pval         test NumSig SetSize                  ID
# 1  3.3464902 5.373078e-08 Sp09D01 ~ L1     32     206           MDD_DLPFC
# 5  2.3947302 4.512138e-06 Sp09D02 ~ L1     42     206           MDD_DLPFC

enriched$gene_DLPFC$k09 |> filter(test == "Sp09D01 ~ L1")
enriched$gene_DLPFC$k09 |> filter(ID == "MDD_DLPFC")

enrichment_check(
    domain = "Sp09D01",
    gene_list = ptsd_gene_list$gene_DLPFC$MDD_DLPFC,
    enrichment = bayesSpace_registration$k09$enrichment
)

enrichment_check(domain = "Sp09D01", gene_list = ptsd_gene_list$gene_DLPFC$MDD_DLPFC)
# $gl
# [1] "TIE1"     "GBP4"     "RHOC"     "TXNIP"    "S100A6"   "RPS27"    "FCGR2A"   "SEMA3G"   "ABCG2"    "CARMN"
# [11] "CDKN1A"   "RPS12"    "UTRN"     "GIMAP6"   "ADGRA2"   "RPS6"     "IFITM2"   "RPL27A"   "RERGL"    "PTPRB"
# [21] "LRP10"    "CRISPLD2" "GADD45B"  "ISYNA1"   "LSR"      "RPS19"    "HSPA12B"  "EDN3"     "JAM2"     "TIMP1"
# [31] "MSN"      "RPL10"
#
# $tab
#
# FALSE  TRUE   Sum
# FALSE 11393   626 12019
# TRUE    174    32   206
# Sum   11567   658 12225
enrichment_check(domain = "Sp09D02", gene_list = ptsd_gene_list$gene_DLPFC$MDD_DLPFC)
# $gl
# [1] "HMGN2"      "TGFBR3"     "RHOC"       "RPS27"      "GAS5"       "RGS8"       "DOCK10"     "PROS1"
# [9] "ALDH1L1"    "EPHB1"      "RGS12"      "TRPC3"      "PCDH18"     "HMGB2"      "RHOBTB3"    "ALDH7A1"
# [17] "RPS12"      "SLC29A4"    "DLX6-AS1"   "TMEM176A"   "CLU"        "CRH"        "RPS6"       "RPL7A"
# [25] "FZD8"       "BAG3"       "IFITM2"     "CD81"       "RPL27A"     "GIHCG"      "LRP10"      "RPL4"
# [33] "RLBP1"      "MT1X"       "KCNJ16"     "DLGAP1-AS1" "NCAN"       "RPS19"      "FTL"        "MSN"
# [41] "AFF2"       "RPL10"
#
# $tab
#
# FALSE  TRUE   Sum
# FALSE 10858  1161 12019
# TRUE    164    42   206
# Sum   11022  1203 12225
enrichment_check(domain = "Sp09D01", gene_list = ptsd_gene_list$gene_DLPFC$PTSD_MDD_DLPFC)
# $gl
# [1] "TNFRSF1B" "TIE1"     "GBP4"     "RHOC"     "TXNIP"    "S100A6"   "SEMA3G"   "TNFSF10"  "ABCG2"    "FKBP5"
# [11] "CDKN1A"   "PNRC1"    "UTRN"     "GIMAP7"   "GIMAP6"   "ADGRA2"   "RPS6"     "MARVELD1" "IFITM2"   "TNFRSF1A"
# [21] "A2M"      "RERGL"    "PTPRB"    "LRP10"    "MT2A"     "CRISPLD2" "RAB34"    "STAT3"    "TUBB6"    "LSR"
# [31] "EMP3"     "HSPA12B"  "EDN3"     "JAM2"     "TIMP1"    "MSN"
#
# $tab
#
# FALSE  TRUE   Sum
# FALSE 11333   622 11955
# TRUE    234    36   270
# Sum   11567   658 12225
enrichment_check(domain = "Sp09D02", gene_list = ptsd_gene_list$gene_DLPFC$PTSD_MDD_DLPFC)
# $gl
# [1] "CAMK2N1"    "HMGN2"      "TGFBR3"     "RHOC"       "GAS5"       "RGS8"       "DOCK10"     "SLC6A11"
# [9] "PROS1"      "ALDH1L1"    "EPHB1"      "RGS12"      "KIT"        "PPM1K"      "TRPC3"      "SMAD1"
# [17] "HMGB2"      "MYO10"      "RHOBTB3"    "ALDH7A1"    "DDR1"       "PNRC1"      "TNS3"       "DLX6-AS1"
# [25] "GAL3ST4"    "TMEM176B"   "TMEM176A"   "CLU"        "CRH"        "RPS6"       "RPL7A"      "FZD8"
# [33] "BAG3"       "IFITM2"     "CD81"       "GSTP1"      "TNFRSF1A"   "GIHCG"      "ACSS3"      "LRP10"
# [41] "ALDH6A1"    "RPL4"       "RLBP1"      "CHD2"       "IQCK"       "MT2A"       "MT1M"       "MT1X"
# [49] "TPPP3"      "RAB34"      "STAT3"      "KCNJ16"     "DLGAP1-AS1" "FTL"        "MSN"
#
# $tab
#
# FALSE  TRUE   Sum
# FALSE 10807  1148 11955
# TRUE    215    55   270
# Sum   11022  1203 12225

# $gene_DLPFC
# n
# MDD_DLPFC           315
# PTSD_MDD_DLPFC      414
# MDD_Male_DLPFC       34
# PTSD_MDD_Male_DLPFC  22


# sgejobs::job_single('05_enrichment_PTSD', create_shell = TRUE, memory = '10G', command = "Rscript 05_enrichment_PTSD.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
