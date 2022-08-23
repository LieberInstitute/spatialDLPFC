library(spatialLIBD)
library(here)

load(file = here::here("processed-data", "rdata", "spe", "08_layer_differential_expression", paste0("parsed_modeling_results_k", k, ".Rdata")))

t0_contrasts_cell <- modeling_results$enrichment
rownames(t0_contrasts_cell) <- t0_contrasts_cell$ensembl
t0_contrasts_cell <- t0_contrasts_cell[, grep("t_stat", colnames(t0_contrasts_cell))]
colnames(t0_contrasts_cell) <- gsub("^t_stat_", "", colnames(t0_contrasts_cell))

ground_truth <- spatialLIBD::fetch_data("modeling_results")

cor_stats_layer <- layer_stat_cor(
    t0_contrasts_cell,
    modeling_results = ground_truth,
    model_type = "enrichment",
    top_n = 100
)

layer_matrix_plot_AS <-
    function(matrix_values,
    matrix_labels = NULL,
    xlabs = NULL,
    layerHeights = NULL,
    mypal = c(
        "white",
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(50)
    ),
    breaks = NULL,
    axis.args = NULL,
    srt = 0,
    mar = c(8, 4 + (max(nchar(rownames(matrix_values))) %/% 3) * 0.5, 4, 2) + 0.1,
    cex = 1.2) {
        ## Create some default values in case the user didn't specify them
        if (is.null(xlabs)) {
            if (is.null(colnames(matrix_values))) {
                xlabs <- paste0("V", seq_len(ncol(matrix_values)))
            } else {
                xlabs <- colnames(matrix_values)
            }
        }

        if (is.null(layerHeights)) {
            layerHeights <- c(0, seq_len(nrow(matrix_values))) * 15
        }

        if (is.null(matrix_labels)) {
            ## Make an empty matrix of labels if none were specified
            matrix_labels <-
                matrix(
                    "",
                    ncol = ncol(matrix_values),
                    nrow = nrow(matrix_values),
                    dimnames = dimnames(matrix_values)
                )
        }

        ## Check inputs
        stopifnot(length(layerHeights) == nrow(matrix_values) + 1)
        stopifnot(length(xlabs) == ncol(matrix_values))
        stopifnot(layerHeights[1] == 0)



        ## For the y-axis labels
        midpoint <- function(x) {
            x[-length(x)] + diff(x) / 2
        }

        ## Make the plot
        par(mar = mar)
        fields::image.plot(
            x = seq(0, ncol(matrix_values), by = 1),
            y = layerHeights,
            z = as.matrix(t(matrix_values)),
            col = mypal,
            xaxt = "n",
            yaxt = "n",
            xlab = "", # x label for plot
            ylab = "",
            breaks = breaks,
            nlevel = length(mypal),
            axis.args = axis.args
        )
        axis(2,
            rownames(matrix_labels),
            at = midpoint(layerHeights),
            las = 1,
            cex.axis = 2 # size of y axis test
        )
        axis(1, rep("", ncol(matrix_values)), at = seq(0.5, ncol(matrix_values) - 0.5))
        text(
            x = seq(0.5, ncol(matrix_values) - 0.5),
            y = -1 * max(nchar(xlabs)) / 2,
            xlabs,
            xpd = TRUE,
            srt = srt,
            cex = 2.5, # size of xaxis text. 1.6 for k = 16, 1.8 for k = 9, 1.0 for k = 28
            adj = c(0.5, 2) # moves x axis number labels (x direction, y direction)
        )
        abline(h = layerHeights, v = c(0, seq_len(ncol(matrix_values))))
        text(
            x = rep(seq(0.5, ncol(matrix_values) - 0.5), each = nrow(matrix_values)),
            y = rep(midpoint(layerHeights), ncol(matrix_values)),
            as.character(matrix_labels),
            cex = cex,
            # cex = cex * 3 / 4,
            font = 4
        )
    }

layer_stat_cor_plot_AS <-
    function(cor_stats_layer,
    max = 0.81,
    min = -max,
    layerHeights = NULL,
    cex = 1.2) {
        ## From https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/dlpfc_snRNAseq_annotation.R
        theSeq <- seq(min, max, by = 0.01)
        my.col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq))

        ## Subset values
        cor_stats_layer[cor_stats_layer <= min] <- min
        cor_stats_layer[cor_stats_layer >= max] <- max

        ## Re-shape the matrix
        mat_vals <- t(cor_stats_layer)

        ## Re-order and shorten names if they match our data
        if (all(rownames(mat_vals) %in% c("WM", paste0("Layer", seq_len(6))))) {
            rownames(mat_vals) <- gsub("ayer", "", rownames(mat_vals))
            mat_vals <- mat_vals[c("WM", paste0("L", rev(seq_len(6)))), , drop = FALSE]

            ## Use our default layer heights also
            if (is.null(layerHeights)) {
                layerHeights <- c(0, 40, 55, 75, 85, 110, 120, 135)
            }
        }

        ## From fields:::imagePlotInfo
        midpoints <- seq(min, max, length.out = length(my.col))
        delta <- (midpoints[2] - midpoints[1]) / 2
        breaks <- c(midpoints[1] - delta, midpoints + delta)

        legend_cuts <- seq(-1, 1, by = 1)
        legend_cuts <- legend_cuts[legend_cuts >= min & legend_cuts <= max]
        axis.args <- list(
            at = legend_cuts,
            labels = legend_cuts,
            cex.axis = 2
        )

        layer_matrix_plot_AS(
            matrix_values = mat_vals,
            matrix_labels = NULL,
            xlabs = NULL,
            layerHeights = layerHeights,
            mypal = my.col,
            breaks = breaks,
            axis.args = axis.args,
            srt = 0, # keep x axis labels normal
            cex = cex
        )
    }

## plot output directory
dir_plots <-
    here::here("plots", "07_spatial_registration")
dir.create(dir_plots, showWarnings = FALSE)

# http://research.libd.org/spatialLIBD/reference/layer_stat_cor_plot.html newer function for plotting


pdf(
    file = here::here(
        "plots",
        "07_spatial_registration",
        paste0(
            "dlpfc_pseudobulked_bayesSpace_vs_mannual_annotations_k",
            k,
            ".pdf"
        )
    ),
    width = 8
)
layer_stat_cor_plot_AS(
    cor_stats_layer,
    max = 1,
    cex = 2.5
)

dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
