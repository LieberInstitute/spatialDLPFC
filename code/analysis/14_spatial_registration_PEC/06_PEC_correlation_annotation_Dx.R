library("spatialLIBD")
library("tidyverse")
library("xlsx")
library("jaffelab")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

## plot dir
plot_dir <- here("plots", "14_spatial_registration_PEC", "06_PEC_correlation_annotation_Dx")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <-
    here(
        "processed-data",
        "rdata",
        "spe",
        "14_spatial_registration_PEC"
    )

load(here(data_dir, "pec_cell_type_tb.Rdata"), verbose = TRUE)

#### Load Layer and k Data  ####
layer_modeling_results <- fetch_data(type = "modeling_results")

paths <-
    list(k09 = "modeling_results_BayesSpace_k09.Rdata", k16 = "modeling_results_BayesSpace_k16.Rdata")

modeling_results <- lapply(paths, function(x) {
    get(load(
        here(
            "processed-data",
            "rdata",
            "spe",
            "07_layer_differential_expression",
            x
        )
    ))
})

modeling_results <-
    c(list(layer = layer_modeling_results), modeling_results)
names(modeling_results)


# layer_anno[[1]] |> left_join(pec_cell_type_tb) |> filter(is.na(cell_type))

correlate_and_annotate <- function(dataset, make_cor_plot = FALSE) {
    #### Load Registration Stats ####
    message(Sys.time(), " - ", dataset)
    registration_stats <- readRDS(here(
        "processed-data",
        "rdata",
        "spe",
        "14_spatial_registration_PEC",
        paste0("registration_stats_Dx_", dataset, ".rds")
    ))

    cor_top100 <- map2(registration_stats, names(registration_stats), function(dx_stats, dx_name) {
        registration_t_stats <-
            dx_stats[, grep("^t_stat", colnames(dx_stats))]

        colnames(registration_t_stats) <-
            paste0(gsub("^t_stat_", "", colnames(registration_t_stats)), "-", dx_name)

        # Fix Cell Types

        #### Correlate with modeling results ####
        cor_top100 <-
            map(
                modeling_results,
                ~ layer_stat_cor(
                    registration_t_stats,
                    .x,
                    model_type = "enrichment",
                    reverse = FALSE,
                    top_n = 100
                )
            )

        return(cor_top100)
    })


    #### Annotate Layers ####
    layer_anno <- map2(cor_top100, names(cor_top100), function(dx_cor, dx) {
        layer_anno <- map2(dx_cor, names(dx_cor), function(cor, name) {
            anno <- annotate_registered_clusters(
                cor_stats_layer = cor,
                confidence_threshold = 0.25,
                cutoff_merge_ratio = 0.10
            )
            colnames(anno) <- gsub("layer", name, colnames(anno))
            return(anno)
        })

        layer_anno <- purrr::reduce(layer_anno, left_join, by = "cluster")

        return(layer_anno)
    })

    layer_anno2 <- do.call("bind_rows", layer_anno) |>
        mutate(Dataset = dataset, .before = 1)

    return(list(cor_top100 = cor_top100, layer_anno = layer_anno2))
}

datasets <-
    datasets <-
    c(
        "CMC",
        "DevBrain-snRNAseq",
        "IsoHuB",
        "LIBD",
        "MultiomeBrain-DLPFC",
        # "PTSDBrainomics",
        "SZBDMulti-Seq",
        "UCLA-ASD"
    )
names(datasets) <- datasets

## Calculate correlations and annotations for each dataset
pe_correlation_annotation <- map(datasets, correlate_and_annotate)

save(pe_correlation_annotation,
    file = here(data_dir, "PEC_correlation_annotation_Dx.Rdata")
)

# load(here(data_dir, "PEC_correlation_annotation_Dx.Rdata"), verbose = TRUE)

#### Save Output to XLSX sheet ####
# key <- data.frame(
#     data = c("annotation", paste0("cor_", names(modeling_results))),
#     description = c(
#         "Annotations of spatial registration",
#         "Correlation values with manual layer annotations",
#         "Correlation values with k09 domains",
#         "Correlation values with k16 domains"
#     )
# )
#
# ## Clear file and write key
# annotation_xlsx <- here(data_dir, "PEC_spatial_annotations.xlsx")
# write.xlsx(
#     key,
#     file = annotation_xlsx,
#     sheetName = "Key",
#     append = FALSE,
#     row.names = FALSE
# )
#
# ## write annotations
# walk2(
#     pe_correlation_annotation,
#     names(pe_correlation_annotation),
#     ~ write.xlsx(
#         .x$layer_anno,
#         file = annotation_xlsx,
#         sheetName = paste0("annotation_", .y),
#         append = TRUE,
#         row.names = FALSE
#     )
# )
#
# ## write correlations
# walk2(pe_correlation_annotation, names(pe_correlation_annotation), function(data, name) {
#     name <- paste0("cor_", name)
#     # message(name)
#     walk2(
#         data$cor_top100,
#         names(data$cor_top100),
#         ~ write.xlsx(
#             t(.x),
#             file = annotation_xlsx,
#             sheetName = paste0(name, "_", .y),
#             append = TRUE,
#             row.names = TRUE
#         )
#     )
# })


# sgejobs::job_single('06_PEC_correlation_annotation_Dx', create_shell = TRUE, memory = '10G', command = "Rscript 06_PEC_correlation_annotation_Dx.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
