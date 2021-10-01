library('basilisk')
library('zellkonverter')
library('spatialLIBD')
library('here')
library('sessioninfo')

anndata_out = here(
    'spagcn', 'processed-data', '02-our_data', 'spe_anndata.h5ad'
)

spe = fetch_data("spe")

write_anndata = function(sce, out_path) {
    invisible(
        basiliskRun(
            fun = function(sce, filename) {
                library('zellkonverter')
                library('reticulate')
                
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
