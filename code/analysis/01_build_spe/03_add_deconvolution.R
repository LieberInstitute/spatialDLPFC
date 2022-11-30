library("spatialLIBD")
library("here")
library("sessioninfo")

Sys.time()
load(
    here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final.Rdata"
    ),
    verbose = TRUE
)
Sys.time()

spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "rdata", "spe", "clustering_results"),
    prefix = ""
)
save(
    spe,
    file = here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters.Rdata"
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
