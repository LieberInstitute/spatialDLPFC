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

split_region <- function(list){
  n <- stringr::str_split(names(list),"_")
  n <- map_chr(n, ~tail(.x, 1))
  list2 <- map(splitit(n),~list[.x])
  return(list2)
}

de_gene2 <- split_region(flatten(de_gene))
names(de_gene2) <- paste0("gene_",names(de_gene2))

de_protein2 <- split_region(flatten(de_protein))
names(de_protein2) <- paste0("protein_",names(de_protein2))

ptsd_genes <- c(de_gene2, de_protein2)
names(ptsd_gene_list)

## convert to list and filter
min_gene <- 10

ptsd_gene_list <- map_depth(ptsd_genes, 2, ~.$gene)
ptsd_gene_list <- map(ptsd_gene_list, ~.x[map_int(.x,length) >20])
ptsd_gene_list <- map_depth(ptsd_gene_list,2, ~ss(.x,"\\."))
map(ptsd_gene_list,~map_int(.x,length))

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


## Find the position for the text title on the plots
# y_text <-
#     map_int(enriched, ~ length(unique(.x$test))) * 15 + c(3, 12, 15, 25, 45)
#

# load(
#   here(
#     "processed-data",
#     "rdata",
#     "spe",
#     "12_spatial_registration_sn",
#     "05_velm_correlation_annotation",
#     "Velmeshev_k9_annotation_details.Rdata"
#   ),
#   verbose = TRUE
# )

pdf(here(dir_plots, "Enrich_PTSD_DLPFC_genes_k09.pdf"), height = 8, width = 10)
gene_set_enrichment_plot_complex(enriched$gene_DLPFC$k09,
                                 gene_count_col = gene_list_count[["gene_DLPFC"]],
                                 gene_count_row = gene_enrichment_count[["k09"]],
                                 anno_title_col = "n ASD DE Genes",
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

#### Interesting gene sets ####

colnames(bayesSpace_registration$k09$enrichment)

gene_check_Sp09D01 <- bayesSpace_registration$k09$enrichment |>
  select(ends_with("Sp09D01"), ensembl, gene) |>
  # filter(fdr_Sp09D01 < 0.05) |>
  mutate(MDD_DLPFC = ensembl %in% ptsd_gene_list$gene_DLPFC$MDD_DLPFC) |>
  dplyr::count(fdr_Sp09D01 < 0.05, MDD_DLPFC)


enrichment_check <- function(domain = "Sp09D01", 
                             gene_list = ptsd_gene_list$gene_DLPFC$MDD_DLPFC, 
                             enrichment = bayesSpace_registration$k09$enrichment){
  
  
  enrichment2 <- enrichment |>
    select(ends_with(domain), ensembl, gene)  |>
    rename_all(~stringr::str_replace(.,paste0("_", domain),"")) |>
    mutate(sig = fdr < 0.01 & t_stat > 0,
           deg = ensembl %in% gene_list)
  
  gene_list <- enrichment2 |>
    filter(deg, sig) |>
    pull(gene)
  
  gene_tab <- table(enrichment2$deg, enrichment2$sig)
  
  
  return(list(gl = gene_list, tab = gene_tab))
}

enriched$gene_DLPFC$k09
#           OR         Pval         test NumSig SetSize                  ID
# 1  3.3464902 5.373078e-08 Sp09D01 ~ L1     32     206           MDD_DLPFC
# 5  2.3947302 4.512138e-06 Sp09D02 ~ L1     42     206           MDD_DLPFC

enriched$gene_DLPFC$k09 |> filter(test == "Sp09D01 ~ L1")
enriched$gene_DLPFC$k09 |> filter(ID == "MDD_DLPFC")

enrichment_check(domain = "Sp09D01")
enrichment_check(domain = "Sp09D02")
enrichment_check(domain = "Sp09D01", gene_list = ptsd_gene_list$gene_DLPFC$PTSD_MDD_DLPFC)
enrichment_check(domain = "Sp09D02", gene_list = ptsd_gene_list$gene_DLPFC$PTSD_MDD_DLPFC)


# $gene_DLPFC
# n
# MDD_DLPFC           315
# PTSD_MDD_DLPFC      414
# MDD_Male_DLPFC       34
# PTSD_MDD_Male_DLPFC  22

gene_enrichment_count$k09
#                 n
# Sp09D01 ~ L1  658
# Sp09D02 ~ L1 1203
# Sp09D03 ~ L2 1551
# Sp09D05 ~ L3 1632
# Sp09D08 ~ L4 2675
# Sp09D04 ~ L5 1921
# Sp09D07 ~ L6  966
# Sp09D06 ~ WM 1858
# Sp09D09 ~ WM  777

gene_list_count

# sgejobs::job_single('05_enrichment_PTSD', create_shell = TRUE, memory = '10G', command = "Rscript 05_enrichment_PTSD.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
