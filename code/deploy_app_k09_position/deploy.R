library("rsconnect")
library("here")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
# load(here("code", "deploy_app_k09", ".deploy_info.Rdata"), verbose = TRUE)
# rsconnect::setAccountInfo(
#     name = deploy_info$name,
#     token = deploy_info$token,
#     secret = deploy_info$secret
# )

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Deploy the app, that is, upload it to shinyapps.io
rsconnect::deployApp(
    appDir = here("code", "deploy_app_k09"),
    appFiles = c(
        "app.R",
        "spe_subset_for_spatialLIBD.rds",
        "sce_pseudo_BayesSpace_k09.rds",
        "modeling_results_BayesSpace_k09.Rdata",
        "sig_genes_subset_k09.Rdata",
        withr::with_dir(here("code", "deploy_app_k09"), dir("www", full.names = TRUE))
    ),
    appName = "spatialDLPFC_Visium_Sp09",
    account = "libd",
    server = "shinyapps.io"
)
