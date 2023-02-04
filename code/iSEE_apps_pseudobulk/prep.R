## Create soft links to the sce_pseudo files
for (k in c(9, 16, 28)) {
    withr::with_dir(
        here("code", "iSEE_apps_pseudobulk", paste0("k", sprintf("%02d", k))),
        system(
            paste0(
                "ln -s ../../../processed-data/rdata/spe/07_layer_differential_expression/sce_pseudo_BayesSpace_k",
                sprintf("%02d", k),
                ".rds sce_pseudo_BayesSpace_k",
                sprintf("%02d", k),
                ".rds"
            )
        )
    )
}
