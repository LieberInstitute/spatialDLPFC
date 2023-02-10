library(purrr)


# Analysis Parameters ---------------------------------------------------
prmt <- expand.grid(
    deconv_method = c("tangram", "cell2location", "SPOTlight")
)


# Helper Function for Setting Up Job --------------------------------------
start.sim <- function(deconv_method) {
    # Compose job name
    job.name <- paste0("Cell_comp_", deconv_method)

    # NOTE:
    ## Job name has to be unique for each of your simulation settings
    ## DO NOT USE GENERIC JOB NAME FOR CONVENIENCE
    job.flag <- paste0("-N ", job.name)

    err.flag <- paste0("-e /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/analysis/15_cell_composition/logs/", job.name, ".txt")

    out.flag <- paste0("-o /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/analysis/15_cell_composition/logs/", job.name, ".txt")

    # Pass simulation parameters to jobs using export flag
    arg.flag <- paste0("-v deconv_method=", deconv_method)

    # Create Jobs
    system(
        paste(
            "qsub", job.flag, err.flag, out.flag, arg.flag,
            "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/analysis/15_cell_composition/01_1_cell_comp_batch_config.sh"
        )
    )
}



# Set up job for all simulation settings ----------------------------------

# for(i in 1:nrow(prmt))
#     do.call(start.sim, prmt[i,,drop =FALSE])

pwalk(
    .l = prmt,
    .f = function(deconv_method, ...) start.sim(deconv_method)
)
