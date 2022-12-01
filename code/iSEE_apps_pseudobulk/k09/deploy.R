library("rsconnect")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
# load(here("code", "deploy_app_k09", ".deploy_info.Rdata"), verbose = TRUE)
# rsconnect::setAccountInfo(
#     name = deploy_info$name,
#     token = deploy_info$token,
#     secret = deploy_info$secret
# )

options(repos = BiocManager::repositories())

rsconnect::deployApp(
    appFiles = c("app.R", "sce_pseudo_BayesSpace_k09.rdsS", "initial.R"),
    appName = "spatialDLPFC_Visium_Sp09_pseudobulk",
    account = "libd",
    server = "shinyapps.io"
)
