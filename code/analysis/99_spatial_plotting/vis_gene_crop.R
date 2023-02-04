vis_gene_crop <-
    function(spe,
    sampleid,
    geneid = "SCGB2A2; ENSG00000110484",
    spatial = TRUE,
    assayname = "logcounts",
    minCount = 0,
    color_scale = "plasma",
    frame_lim_df,
    image_id = "lowres",
    alpha = NA,
    cont_colors = c("aquamarine4", "springgreen", "goldenrod", "red"),
    point_size = 2,
    legend_overlap = FALSE,
    ...) {
        if (color_scale == "viridis") {
            cont_colors <- viridisLite::viridis(21)
        } else if (color_scale == "plasma") {
            cont_colors <- viridisLite::plasma(21)
        }

        spe_sub <- spe[, spe$sample_id == sampleid]
        d <- as.data.frame(cbind(colData(spe_sub), SpatialExperiment::spatialCoords(spe_sub)), optional = TRUE)

        if (geneid %in% colnames(colData(spe_sub))) {
            d$COUNT <- colData(spe_sub)[[geneid]]
        } else if (geneid %in% rowData(spe_sub)$gene_search) {
            d$COUNT <-
                assays(spe_sub)[[assayname]][which(rowData(spe_sub)$gene_search == geneid), ]
        } else if (geneid %in% rownames(spe_sub)) {
            d$COUNT <- assays(spe_sub)[[assayname]][which(rownames(spe_sub) == geneid), ]
        } else {
            stop("Could not find the 'geneid' ", geneid, call. = FALSE)
        }
        d$COUNT[d$COUNT <= minCount] <- NA

        legend_title <- paste(
            if (!geneid %in% colnames(colData(spe_sub))) {
                assayname
            } else {
                NULL
            },
            "min >", minCount
        )

        p <- vis_gene_p_crop(
            spe = spe_sub,
            d = d,
            sampleid = sampleid,
            spatial = spatial,
            frame_lim_df = frame_lim_df,
            legend_title = legend_title,
            title = paste(
                sampleid,
                geneid,
                ...
            ),
            viridis = viridis,
            image_id = image_id,
            alpha = alpha,
            cont_colors = cont_colors,
            point_size = point_size,
            legend_overlap = legend_overlap
        )

        return(p)
    }

vis_gene_p_crop <-
    function(spe,
    d,
    sampleid,
    spatial,
    title,
    viridis = TRUE,
    image_id = "lowres",
    frame_lim_df,
    legend_title = "Test",
    alpha = NA,
    cont_colors = if (viridis) viridisLite::viridis(21) else c("aquamarine4", "springgreen", "goldenrod", "red"),
    point_size = 2,
    legend_overlap = FALSE) {
        ## Some variables
        pxl_row_in_fullres <- pxl_col_in_fullres <- key <- COUNT <- NULL
        # stopifnot(all(c("pxl_col_in_fullres", "pxl_row_in_fullres", "COUNT", "key") %in% colnames(d)))

        scale_factor <- SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id)
        frame_lims <- ceiling(frame_lim_df[match(sampleid, frame_lim_df$sample_id), -1] * scale_factor)

        ## call image
        img <- SpatialExperiment::imgRaster(spe, sample_id = sampleid, image_id = image_id)

        ## keep crop in frame
        frame_lims$y_min <- ifelse(frame_lims$y_min < 1, 1, frame_lims$y_min)
        frame_lims$y_max <- ifelse(frame_lims$y_max > nrow(img), nrow(img), frame_lims$y_max)
        frame_lims$x_min <- ifelse(frame_lims$x_min < 1, 1, frame_lims$x_min)
        frame_lims$x_max <- ifelse(frame_lims$x_max > ncol(img), ncol(img), frame_lims$x_max)
        ## crop image
        img <- img[frame_lims$y_min:frame_lims$y_max, frame_lims$x_min:frame_lims$x_max]

        p <-
            ggplot(
                d,
                aes(
                    x = (pxl_col_in_fullres * scale_factor) - frame_lims$x_min,
                    y = (pxl_row_in_fullres * scale_factor) - frame_lims$y_min,
                    fill = COUNT,
                    color = COUNT,
                    key = key
                )
            )

        if (spatial) {
            grob <- grid::rasterGrob(img, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
            p <-
                p + geom_spatial(
                    data = tibble::tibble(grob = list(grob)),
                    aes(grob = grob),
                    x = 0.5,
                    y = 0.5
                )
        }

        ## From https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/layer_marker_genes_plots.R
        # add.alpha('black', 0.175)
        # black
        # "#0000002D"

        p <- p +
            geom_point(
                shape = 21,
                size = point_size,
                stroke = 0,
                alpha = alpha
            ) +
            coord_cartesian(expand = FALSE)

        p <- p + scale_fill_gradientn(
            name = legend_title,
            colors = cont_colors,
            na.value = c("black" = "#0000002D")
        ) +
            scale_color_gradientn(
                name = legend_title,
                colors = cont_colors,
                na.value = c("black" = "#0000002D")
            )



        p <- p +
            coord_fixed() +
            xlim(0, ncol(img)) +
            ylim(nrow(img), 0) +
            xlab("") + ylab("") +
            labs(fill = NULL, color = NULL) +
            # ggtitle(title) +
            theme_set(theme_bw(base_size = 20)) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.position = "bottom",
                legend.box.spacing = unit(0, "pt")
            )

        if (legend_overlap) {
            p <- p + theme(
                legend.position = c(0.5, 0.08),
                legend.box = "horizontal",
                legend.direction = "horizontal"
            )
        }

        return(p)
    }
