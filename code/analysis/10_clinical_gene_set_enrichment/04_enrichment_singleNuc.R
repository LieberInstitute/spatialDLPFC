# library(sgejobs)
# sgejobs::job_single(
#     name = "02_enrichment_HumanPilot_sets",
#     create_shell = TRUE,
#     queue = "bluejay",
#     memory = "5G"
# )
# To execute the script builder, use: qsub 02_enrichment_HumanPilot_sets.sh

library("here")
library("purrr")
library("dplyr")
## Due to this recent change
## https://github.com/LieberInstitute/spatialLIBD/commit/cefc7db61a16e80b16c14f0df30b40701cfc788c
stopifnot(packageVersion("spatialLIBD") >= "1.10.1")
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
    "10_clinical_gene_set_enrichment"
)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully

dir_plots <- here::here(
    "plots",
    "10_clinical_gene_set_enrichment",
    "04_enrichment_singleNuc"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## Load the gene sets from Singel Nucleus datasets
load(here(dir_rdata, "gene_sets_SingleNuc.Rdata"), verbose = TRUE)
names(geneList)


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
  map(geneList, function(gl){
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


## Save for later
save(enriched, file = here(dir_rdata, "enriched_singleNuc.Rdata"))

## Some sets are much larger than others
map_depth(enriched,2 ,~ filter(.x, test == unique(test)[1])[, c("SetSize", "ID")])


## source updated enrichment plot 
source(here("code", "analysis","10_clinical_gene_set_enrichment", "gene_set_enrichment_plot_complex.R"))

## how many genes for each domain or DE set, 
gene_enrichment_count <- map(bayesSpace_registration, function(r){
  en_count <- get_gene_enrichment_count(r)
  ##reorder to match layer_combo
  rownames(en_count) <- bayes_anno$layer_combo[match(rownames(en_count), bayes_anno$test)]
  layer_order <- bayes_anno$layer_combo[bayes_anno$layer_combo %in% rownames(en_count)]
  return(en_count[layer_order,,drop = FALSE])
  })

gene_list_count <- map(geneList, get_gene_list_count)


## Find the position for the text title on the plots
# y_text <-
#     map_int(enriched, ~ length(unique(.x$test))) * 15 + c(3, 12, 15, 25, 45)
# 

load(here("processed-data",
          "rdata",
          "spe",
          "12_spatial_registration_sn", 
          "05_velm_correlation_annotation", 
          "Velmeshev_k9_annotation_details.Rdata"),
     verbose = TRUE)

pdf(here(dir_plots, "Enrich_Velmeshev_k09.pdf"), height = 8, width = 10)
  gene_set_enrichment_plot_complex(enriched$Velmeshev$k09,
                           gene_count_col = gene_list_count[["Velmeshev"]],
                           gene_count_row = gene_enrichment_count[["k09"]],
                           anno_title_col = "n ASD DE Genes",
                           anno_title_row = "n Domain\nGenes",
                           column_order = velm_col_order_k9)
dev.off()

pdf(here(dir_plots, "Enrich_Velmeshev_k09anno.pdf"), height = 8, width = 10)
gene_set_enrichment_plot_complex(enriched$Velmeshev$k09,
                                 gene_count_col = gene_list_count[["Velmeshev"]],
                                 gene_count_row = gene_enrichment_count[["k09"]],
                                 anno_title_col = "n ASD DE Genes",
                                 anno_title_row = "n Domain\nGenes",
                                 column_order = velm_col_order_k9,
                                 anno_add = anno_matrix_k9)
dev.off()


walk2(enriched, names(enriched), function(enriched, ds_name){
  
  pdf(here(dir_plots, paste0("Enrich_", ds_name, "_FDR05.pdf")), height = 8, width = 9)
  map2(enriched, names(enriched), function(x, k) {
          message(ds_name, " - ", k)
    
          print(gene_set_enrichment_plot_complex(x,
               gene_count_col = gene_list_count[[ds_name]],
               gene_count_row = gene_enrichment_count[[k]],
               anno_title_col = "n DE Genes",
               anno_title_row = "n Domain\nGenes"))
    
  })
  dev.off()
  
})


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

