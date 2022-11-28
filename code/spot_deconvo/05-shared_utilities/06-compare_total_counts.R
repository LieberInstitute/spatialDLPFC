suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("cowplot"))

spe_IF_in <- here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)
spe_nonIF_in <- here(
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
)

plot_dir <- here("plots", "spot_deconvo", "05-shared_utilities")

#   "Ground-truth" cell counts from cellpose + trained classification tree
actual_paths <- here(
    "processed-data", "spot_deconvo", "02-cellpose", "{sample_id}", "clusters.csv"
)

################################################################################
#   Functions
################################################################################

plot_counts <- function(spe_IF, spe_nonIF, plot_name) {
    plot_list <- list()
    i <- 1

    for (sample_id in unique(spe_IF$sample_id)) {
        plot_list[[i]] <- vis_grid_gene(
            spe_nonIF[, spe_nonIF$sample_id == sample_id],
            geneid = "count",
            return_plots = TRUE
        )[[1]] +
            labs(title = paste0(sample_id, " (non-IF): VistoSeg counts"))

        plot_list[[i + 1]] <- vis_grid_gene(
            spe_IF[, spe_IF$sample_id == sample_id],
            geneid = "counts_cellpose",
            return_plots = TRUE
        )[[1]] +
            labs(title = paste0(sample_id, " (IF): cellpose counts"))

        i <- i + 2
    }

    pdf(
        file.path(plot_dir, plot_name),
        width = 14, height = 7 * length(unique(spe_IF$sample_id))
    )
    print(plot_grid(plotlist = plot_list, ncol = 2))
    dev.off()
}

################################################################################
#   Main
################################################################################

#   Load objects
spe_IF <- readRDS(spe_IF_in)
load(spe_nonIF_in, verbose = TRUE)
spe_nonIF <- spe
rm(spe)
gc()

#-------------------------------------------------------------------------------
#   Add cellpose counts to IF SPE
#-------------------------------------------------------------------------------

#   Read in cellpose counts
counts_list <- list()
for (sample_id in unique(spe_IF$sample_id)) {
    #   Read in the "ground-truth" counts for this sample
    actual_path <- sub("\\{sample_id\\}", sample_id, actual_paths)
    counts_list[[sample_id]] <- read.csv(actual_path)
}

counts <- do.call(rbind, counts_list)

stopifnot(all(counts$key %in% spe_IF$key))
stopifnot(all(spe_IF$key %in% counts$key))

#   Add to IF SPE
spe_IF$counts_cellpose <- counts$n_cells[match(spe_IF$key, counts$key)]

#-------------------------------------------------------------------------------
#   Plot cell counts for adjacent sections
#-------------------------------------------------------------------------------

#   Use the same naming convention for sample names
spe_nonIF$sample_id <- tolower(spe_nonIF$sample_id)
spe_IF$sample_id <- tolower(sub("_IF", "", spe_IF$sample_id))

#   Subset nonIF samples to those taken from the same sections as from the IF
#   samples
spe_nonIF <- spe_nonIF[, spe_nonIF$sample_id %in% unique(spe_IF$sample_id)]

plot_counts(spe_IF, spe_nonIF, "cellpose_v_vistoseg_counts_raw.pdf")

#   Cap VistoSeg at 15 cells, since there are regions throwing off the scale of
#   gigantic counts

a <- spe_nonIF
a$count[a$count > 15] <- 15

plot_counts(spe_IF, a, "cellpose_v_vistoseg_counts_trimmed.pdf")
