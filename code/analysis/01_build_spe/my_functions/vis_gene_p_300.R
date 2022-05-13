vis_gene_p_300 <-
  function(spe,
           d,
           sampleid,
           spatial,
           title,
           viridis = TRUE,
           image_id = "lowres",
           alpha = 1,
           cont_colors = if (viridis) viridisLite::viridis(21) else c("aquamarine4", "springgreen", "goldenrod", "red"),
           point_size = 1.25) {
    
    ## Some variables
    pxl_row_in_fullres <- pxl_col_in_fullres <- key <- COUNT <- NULL
    # stopifnot(all(c("pxl_col_in_fullres", "pxl_row_in_fullres", "COUNT", "key") %in% colnames(d)))
    img <- SpatialExperiment::imgRaster(spe, sample_id = sampleid, image_id = image_id)
    
    p <-
      ggplot(
        d,
        aes(
          x = pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id),
          y = pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id),
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
      colors = cont_colors,
      na.value = c("black" = "#0000002D"),
      limits=c(0,300)
    ) +
      scale_color_gradientn(
        colors = cont_colors,
        na.value = c("black" = "#0000002D"),
        limits=c(0,300)
      )

    p <- p +
      xlim(0, ncol(img)) +
      ylim(nrow(img), 0) +
      xlab("Testing") + ylab("") +
      labs(fill = NULL, color = NULL) +
      ggtitle(title) +
      theme_set(theme_bw(base_size = 20)) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.95, 0.10)
      )
    return(p)
  }

