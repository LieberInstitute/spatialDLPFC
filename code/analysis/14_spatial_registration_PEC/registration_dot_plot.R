

#' Plot a dotplot of spatail registation annotations
#'
#' @param annotation_df data.frame from annotate_registered_clusters
#' @param color_by name of column to color dots by
#'
#' @return
#' @export
#'
#' @examples
#' 
#' test_anno_a <- data.frame(cluster = c("Astro", "Excit", "Oligo"),
#' layer_confidence = "good",
#' layer_label = paste0("SpD", 1:3),
#' Dataset = "A"
#' )
#' 
#' test_anno_b <- data.frame(cluster = c("Astro", "Excit", "Oligo"),
#' layer_confidence = "good",
#' layer_label = paste0("SpD", c(1, 1, 3)),
#' Dataset = "B"
#' )
#' 
#' ## Which layer lables match previous annotations? 
#' ct_anno <- data.frame(cluster = c("Astro", "Excit", "Oligo"),
#' layer_label = paste0("SpD", 1:3),
#' Match = TRUE
#' )
#' 
#' ## combine datasets
#' test_anno <- dplyr::bind_rows(test_anno_a, test_anno_b) |>
#' dplyr::mutate(layer_label = factor(layer_label, levels = paste0("SpD", 1:4)))
#' 
#' registration_dot_plot(test_anno_a)
#' registration_dot_plot(test_anno_a, ct_anno = ct_anno)
#' registration_dot_plot(test_anno)
#' registration_dot_plot(test_anno, ct_anno = ct_anno)
#' 
registration_dot_plot <- function(annotation_df, 
                                  color_by = "Dataset",
                                  grid_fill = list(`TRUE` = "grey80", `FALSE` = "white"),
                                  ct_anno = NULL){
  
  if(is.factor(annotation_df$layer_label)){
    layer_lable_levels <- levels(annotation_df$layer_label)
  } else { 
    layer_lable_levels <- unique(annotation_df$layer_label)
  }
  
  tile_data <- tidyr::expand_grid(cluster = annotation_df$cluster, layer_label = layer_lable_levels) 
  
  if(is.null(ct_anno)){
    tile_data$Match <- FALSE
  } else {
    tile_data <- tile_data |>
      dplyr::left_join(ct_anno) |>
      tidyr::replace_na(list(Match = FALSE))
  }
  
  ex_dotplot <- ggplot2::ggplot(tile_data, 
                                ggplot2::aes(x = cluster, y = layer_label)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Match), color = "grey10") + 
    ggplot2::geom_point(data = annotation_df,
                        ggplot2::aes(color = .data[[color_by]]), 
               position = ggplot2::position_dodge(width = .8)) +
    ggplot2::scale_fill_manual(values = grid_fill, guide="none") +
    ggplot2::theme_bw()
  
  return(ex_dotplot)
}
