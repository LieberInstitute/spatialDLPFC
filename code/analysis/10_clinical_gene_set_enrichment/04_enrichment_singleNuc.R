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
          test = factor(layer_combo, levels = rev(bayes_anno$layer_combo[bayes_anno$layer_combo %in% layer_combo]))
        ) |>
        select(-c(layer_combo, Annotation, fdr_cut, model_type))
    )
  })


## Save for later
save(enriched, file = here(dir_rdata, "enriched_singleNuc.Rdata"))

## Some sets are much larger than others
map_depth(enriched,2 ,~ filter(.x, test == unique(test)[1])[, c("SetSize", "ID")])

## Find the position for the text title on the plots
y_text <-
    map_int(enriched, ~ length(unique(.x$test))) * 15 + c(3, 12, 15, 25, 45)


walk2(enriched, names(enriched), function(enriched, ds_name){
  
  pdf(here(dir_plots, paste0(ds_name, "_FDR05.pdf")), height = 8, width = 9)
  walk2(enriched, y_text, function(x, ypos) {

    m <- match(unique(x$ID), x$ID)
    gene_set_enrichment_plot(x, xlabs = x$ID[m])

    # text(
    #   x = 3,
    #   y = ypos,
    #   c(paste0(ds_name, " et al.")),
    #   xpd = TRUE,
    #   cex = 2.5,
    #   font = 2
    # )
  })
  dev.off()
  
})

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

