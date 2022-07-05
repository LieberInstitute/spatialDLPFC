library("here")
library("sessioninfo")
library("spatialLIBD")
library("RColorBrewer")
library(ggplot2)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#load spe object and clusters
load(file = here::here("processed-data","rdata","spe","01_build_spe","spe_filtered_final_with_clusters.Rdata"),verbose = TRUE)

mycolors <- Polychrome::palette36.colors(k)
names(mycolors) <- sort(unique(colData(spe)[[paste0("bayesSpace_harmony_",k)]]))
bayesSpace_name <- paste0("bayesSpace_harmony_", k)
sample_ids <- unique(colData(spe)$sample_id)

vis_clus_p_AS <-
  function(spe,
           d,
           clustervar,
           sampleid,
           colors,
           spatial,
           title,
           image_id = "lowres",
           alpha = 1,
           point_size = 1.25) {
    
    ## Some variables
    pxl_row_in_fullres <- pxl_col_in_fullres <- key <- NULL
    # stopifnot(all(c("pxl_col_in_fullres", "pxl_row_in_fullres", "key") %in% colnames(d)))
    
    if (clustervar %in% c(
      "layer_guess",
      "layer_guess_reordered",
      "layer_guess_reordered_short",
      "spatialLIBD"
    )) {
      title <- gsub(clustervar, "LIBD Layers", title)
    }
    img <- SpatialExperiment::imgRaster(spe, sample_id = sampleid, image_id = image_id)
    
    p <- ggplot(
      d,
      aes(
        x = pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id),
        y = pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id),
        fill = factor(!!sym(clustervar)),
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
    p <- p +
      geom_point(
        shape = 21,
        size = point_size,
        stroke = 0,
        colour = "transparent",
        alpha = alpha
      ) +
      coord_cartesian(expand = FALSE) +
      scale_fill_manual(values = colors) +
      xlim(0, ncol(img)) +
      ylim(nrow(img), 0) +
      xlab("") + ylab("") +
      labs(fill = NULL) +
      guides(fill = guide_legend(override.aes = list(size = 3))) +
      ggtitle(title) +
      theme_set(theme_bw(base_size = 20),theme(legend.key=element_blank())) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
    return(p)
  }

vis_clus_AS <- function(spe,
                     sampleid,
                     clustervar,
                     colors = c(
                       "#b2df8a",
                       "#e41a1c",
                       "#377eb8",
                       "#4daf4a",
                       "#ff7f00",
                       "gold",
                       "#a65628",
                       "#999999",
                       "black",
                       "grey",
                       "white",
                       "purple"
                     ),
                     spatial = TRUE,
                     image_id = "lowres",
                     alpha = 1,
                     point_size = 1.25,
                     ...) {
  spe_sub <- spe[, spe$sample_id == sampleid]
  d <- as.data.frame(cbind(colData(spe_sub), SpatialExperiment::spatialCoords(spe_sub)), optional = TRUE)
  
  vis_clus_p_AS(
    spe = spe_sub,
    d = d,
    clustervar = clustervar,
    sampleid = sampleid,
    spatial = spatial,
    title = paste0(sampleid, ...),
    colors = get_colors(colors, d[, clustervar]),
    image_id = image_id,
    alpha = alpha,
    point_size = point_size
  )
}

pdf(file = here::here("plots","03_BayesSpace",paste0("polychrome_vis_clus_bayesSpace_harmony_square",k,".pdf")))
for (i in seq_along(sample_ids)){
  my_plot <- vis_clus(
    spe = spe,
    clustervar = bayesSpace_name,
    sampleid = sample_ids[i],
    colors = mycolors
  )
  print(my_plot)
}
dev.off()