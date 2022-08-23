library("basilisk")
library("zellkonverter")
library("spatialLIBD")
library("here")
library("sessioninfo")

anndata_out <- here(
    "spagcn", "processed-data", "02-our_data_tutorial", "spe_anndata.h5ad"
)

spe <- fetch_data("spe")

#  Append 'spatialCoords' slot to 'colData', since in
#  conversion we're treating the spatialExperiment object as if it is a
#  singleCellExperiment, which doesn't have that additional slot. Note that
#  in this case, the 'spatialData' slot only contains duplicate information,
#  which is why we don't append it
colData(spe) <- cbind(colData(spe), spatialCoords(spe))

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
