library("basilisk")
library("zellkonverter")
library("spatialLIBD")
library("here")
library("sessioninfo")

spe_in <- here(
    "spagcn", "raw-data", "04-spatial_dlpfc", "spe_final.Rdata"
)

anndata_out <- here(
    "spagcn", "processed-data", "04-spatial_dlpfc", "spe_anndata.h5ad"
)

load(spe_in, verbose = TRUE)

#  Append 'spatialCoords' and 'spatialData' slots to 'colData', since in
#  conversion we're treating the spatialExperiment object as if it is a
#  singleCellExperiment, which doesn't have those additional slots.
colData(spe) <- cbind(colData(spe), spatialData(spe), spatialCoords(spe))

write_anndata <- function(sce, out_path) {
    invisible(
        basiliskRun(
            fun = function(sce, filename) {
                library("zellkonverter")
                library("reticulate")

                # Convert SCE to AnnData:
                adata <- SCE2AnnData(sce)

                #  Write AnnData object to disk
                adata$write(filename = filename)

                return()
            },
            env = zellkonverterAnnDataEnv,
            sce = sce,
            filename = out_path
        )
    )
}

write_anndata(spe, anndata_out)

session_info()
