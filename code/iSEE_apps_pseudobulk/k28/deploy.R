
library("rsconnect")

source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "spe_pseudobulk_bayesSpace_normalized_filtered_region_k28.RDS", "initial.R"),
    appName = "spatialDLPFC_Visium_Sp28_pseudobulk",
    account = "libd",
    server = "shinyapps.io"
)
