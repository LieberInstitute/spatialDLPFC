
library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")

load("sce_for_iSEE_LS.rda", verbose = TRUE)

stopifnot(all(unique(sce.ls.small$cellType.final) %in% names(cell_cols.clean)))

## Don't run this on app.R since we don't want to run this every single time
# lobstr::obj_size(sce.ls.small)
# 876.33 MB

source("initial.R", print.eval = TRUE)

## From https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/shiny_apps/tran2021_AMY/app.R#L10-L14
## Related to https://github.com/iSEE/iSEE/issues/568
colData(sce.ls.small) <- cbind(
  colData(sce.ls.small)[, !colnames(colData(sce.ls.small)) %in% c("Sample", "cellType.final")],
  colData(sce.ls.small)[, c("cellType.final", "Sample")]
)

sce.ls.small$Sample <- as.factor(sce.ls.small$Sample)

sce.ls.small <- registerAppOptions(sce.ls.small, color.maxlevels = length(cell_cols.clean))
iSEE(
    sce.ls.small,
    appTitle = "mm_LS_2022",
    initial = initial,
    colormap = ExperimentColorMap(colData = list(
        Sample = function(n) {
            cols <- paletteer::paletteer_d(
                palette = "RColorBrewer::Dark2",
                n = length(unique(sce.ls.small$Sample))
            )
            cols <- as.vector(cols)
            names(cols) <- levels(sce.ls.small$Sample)
            return(cols)
        },
        cellType.final = function(n) {
            return(cell_cols.clean)
        }
    ))
)
