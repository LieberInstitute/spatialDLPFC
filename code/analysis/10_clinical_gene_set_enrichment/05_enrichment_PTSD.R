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
map_depth(ptsd_gene_list,2, length)


map_int(degs_dlpfc_list, length)
# MDD_DLPFC      PTSD_MDD_DLPFC      MDD_Male_DLPFC PTSD_MDD_Male_DLPFC 
# 315                 414                  34                  22

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

# sgejobs::job_single('05_enrichment_PTSD', create_shell = TRUE, memory = '10G', command = "Rscript 05_enrichment_PTSD.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
