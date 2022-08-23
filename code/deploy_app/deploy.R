library("rsconnect")
library("here")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
# load(here("code", "deploy_app", ".deploy_info.Rdata"), verbose = TRUE)
# rsconnect::setAccountInfo(
#     name = deploy_info$name,
#     token = deploy_info$token,
#     secret = deploy_info$secret
# )

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Deploy the app, that is, upload it to shinyapps.io
rsconnect::deployApp(
    appDir = here("code", "deploy_app"),
    appFiles = c(
        "app.R",
        "spe_subset.Rdata",
        "spe_pseudobulk_bayesSpace_normalized_filtered_cluster_k9.RDS",
        "parsed_modeling_results_k9.Rdata",
        "sig_genes_subset.Rdata"
    ),
    appName = "k9_spatialDLPFC_Spangler2022",
    account = "libd",
    server = "shinyapps.io"
)
