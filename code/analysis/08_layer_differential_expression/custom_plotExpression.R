
#' Plot SCE Expression of marker genes over clusters
#'
#' @param sce A SingleCellExperiment object containing expression values
#' @param genes list of genes to plot
#' @param assay name of assay to plot
#' @param cat name of category
#' @param fill_colors optional color pallet for cat
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' custom_plotExpression(example_sce, c("Gene_0001", "Gene_0002", "Gene_0003", "Gene_0004"), cat = "Mutation_Status")
#' custom_plotExpression(example_sce,
#'     c("Gene_0001", "Gene_0002", "Gene_0003", "Gene_0004"),
#'     cat = "Mutation_Status",
#'     fill_colors = c(negative = "green", positive = "pink")
#' )
custom_plotExpression <- function(sce, genes, assay = "logcounts", cat, highlight_sample = "None", fill_colors = NULL, title = NULL) {
    cat_df <- as.data.frame(colData(sce))[, c("sample_id", cat), drop = FALSE]
    
    expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes, ,drop = FALSE]))
    cat <- cat_df[expression_long$Var2, ]
    colnames(cat) <- c("sample_id", "cat") 
    expression_long <- cbind(expression_long, cat)
    
    # return(expression_long)

    expression_violin <- ggplot(data = expression_long, aes(x = cat, y = value, fill = cat)) +
        geom_violin(scale = "width") +
        geom_jitter(aes(color = sample_id == highlight_sample), size = 1) +
        facet_wrap(~Var1) +
        labs(
            y = paste0("Expression (", assay, ")"),
            title = title
        ) +
        theme_bw() +
        theme(
            legend.position = "None", 
            # axis.title.x = element_blank(),
            # axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(face = "italic"),
            text = element_text(size = 15)
            
        ) +
        stat_summary(
            fun = median,
            # fun.min = median,
            # fun.max = median,
            geom = "crossbar",
            width = 0.3
        ) +
      scale_color_manual(values = c(`TRUE` = "yellow", `FALSE` = "black"))

    if (!is.null(fill_colors)) expression_violin <- expression_violin + scale_fill_manual(values = fill_colors)

    # expression_violin
    return(expression_violin)
}



