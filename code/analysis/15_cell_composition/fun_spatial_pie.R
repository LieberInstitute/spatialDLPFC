# Helper function to generate vis_spatial_pie
vis_scatter_pie <- function(spe,
    plot_df,
    sampleid,
    image_id,
    cell_names,
    auto_crop = TRUE,
    spatial = FALSE) {
    pxl_row_in_fullres <- spatialCoords(spe)[, "pxl_row_in_fullres"]
    pxl_col_in_fullres <- spatialCoords(spe)[, "pxl_col_in_fullres"]
    img <- SpatialExperiment::imgRaster(spe,
        sample_id = sampleid,
        image_id = image_id
    )

    ## Crop the image if needed
    if (auto_crop) {
        frame_lims <- spatialLIBD::frame_limits(spe,
            sampleid = sampleid,
            image_id = image_id
        )
        img <- img[frame_lims$y_min:frame_lims$y_max, frame_lims$x_min:frame_lims$x_max]
        adjust <- list(x = frame_lims$x_min, y = frame_lims$y_min)
    } else {
        adjust <- list(x = 0, y = 0)
    }

    plot_df$scaled_x <- pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id) - adjust$x
    plot_df$scaled_y <- pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id) - adjust$y

    ret_plot <- ggplot2::ggplot(
        # data = plot_df,
        # ggplot2::aes(
        #     x = scaled_x,
        #     y = scaled_y
        # )
    )



    if (spatial) {
        grob <- grid::rasterGrob(img, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
        ret_plot <-
            ret_plot + spatialLIBD::geom_spatial(
                data = tibble::tibble(grob = list(grob)),
                aes(grob = grob),
                x = 0.5,
                y = 0.5
            )

        # ggsave("~/tmp.pdf",
        #        ret_plot,
        #        height = 7.5,
        #        width = 9)
    }

    ret_plot <- ret_plot +
        scatterpie::geom_scatterpie(
            data = plot_df,
            ggplot2::aes(
                x = scaled_x,
                y = -1 * scaled_y
            ),
            col = cell_names,
            color = NA,
            pie_scale = 0.35
        ) +
        # coord_fixed(expand = FALSE) +
        # coord_fixed(expand = FALSE, ratio = 1) +
        ggplot2::coord_fixed(ratio = 1) +
        ggplot2::scale_fill_manual(
            limits = fill_vals,
            labels = fill_labs,
            values = cell_type_colors
        ) +
        # xlim(0, ncol(img)) +
        # ylim(nrow(img), 0) +
        xlab("") + ylab("") +
        labs(fill = NULL) +
        # theme_bw()
        theme_set(theme_bw(base_size = 20)) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.box.spacing = unit(0, "pt")
        )

    return(ret_plot)
}


# Function to generate spatial scatter pie plot
vis_spatial_pie <- function(spe,
    res, # Cell Resolution
    deconv, # Deconvolution method
    image_id,
    sample_id,
    spatial = FALSE) {
    plot_df <- colData(spe) |>
        data.frame() |>
        dplyr::select(
            starts_with(paste0(res, "_", deconv)), # cell_names
            row, col # Spatial Coordinates
        ) |>
        rename_with(.fn = str_remove, pattern = paste0(res, "_", deconv, "_"))
    all_cell_names <- grep(
        pattern = paste0(res, "_", deconv),
        x = names(colData(spe)),
        value = TRUE
    ) |>
        str_remove(pattern = paste0(res, "_", deconv, "_"))

    ret_plot <- vis_scatter_pie(
        spe,
        plot_df = plot_df,
        sampleid = sample_id,
        image_id = image_id,
        cell_names = all_cell_names,
        spatial = FALSE
    )
}
