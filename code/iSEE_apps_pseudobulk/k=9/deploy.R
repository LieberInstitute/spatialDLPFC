
library("rsconnect")

source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "spe_pseudobulk_bayesSpace_normalized_filtered_region_k9.RDS", "initial.R"),
    appName = "Spangler2022_pseudobulk_cluster_k09",
    account = "libd",
    server = "shinyapps.io"
)
