
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
#' my_plotExpression(example_sce, c("Gene_0001", "Gene_0002", "Gene_0003", "Gene_0004"), cat = "Mutation_Status")
#' my_plotExpression(example_sce,
#'     c("Gene_0001", "Gene_0002", "Gene_0003", "Gene_0004"),
#'     cat = "Mutation_Status",
#'     fill_colors = c(negative = "green", positive = "pink")
#' )
my_plotExpression <- function(sce, genes, assay = "logcounts", cat = "cellType", fill_colors = NULL, title = NULL) {
    cat_df <- as.data.frame(colData(sce))[, cat, drop = FALSE]
    expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes, ]))

    cat <- cat_df[expression_long$Var2, ]
    expression_long <- cbind(expression_long, cat)

    expression_violin <- ggplot(data = expression_long, aes(x = cat, y = value, fill = cat)) +
        geom_violin(scale = "width") +
        facet_wrap(~Var1, ncol = 2) +
        labs(
            y = paste0("Expression (", assay, ")"),
            title = title
        ) +
        theme_bw() +
        theme(
            legend.position = "None", axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(face = "italic")
        ) +
        stat_summary(
            fun = median,
            # fun.min = median,
            # fun.max = median,
            geom = "crossbar",
            width = 0.3
        )

    if (!is.null(fill_colors)) expression_violin <- expression_violin + scale_fill_manual(values = fill_colors)

    # expression_violin
    return(expression_violin)
}


my_plotMarkers <- function(sce, marker_list, assay = "logcounts", cat = "cellType", fill_colors = NULL, pdf_fn) {
    message("plotting: ", pdf_fn)
    pdf(pdf_fn, height = 6, width = 8)

    for (i in 1:length(marker_list)) {
        message(names(marker_list)[[i]])

        m <- marker_list[[i]]
        markers <- m[m %in% rownames(sce)]

        # if(m != markers) message("Missing...",paste(m[!m %in% markers], collapse = ", "))
        if(!identical(m, markers)) message("Missing markers...")
        print(
            my_plotExpression(sce,
                genes = markers,
                title = names(marker_list)[[i]],
                cat = cat,
                fill_colors = fill_colors
            )
        )
    }

    dev.off()
}


#' Plot SCE Expression of clusters over marker genes
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
#' example_sce <- mockSCE()
#' example_sce <- logNormCounts(example_sce)
#' marker_list <- list("CT1" = c("Gene_0001", "Gene_0002"), "CT2" = c("Gene_0003", "Gene_0004"))
#'
#' my_plotExpression_flip(example_sce, marker_list)
#'
#' my_plotExpression(example_sce,
#'     c("Gene_0001", "Gene_0002", "Gene_0003", "Gene_0004"),
#'     cat = "Mutation_Status",
#'     fill_colors = c(negative = "green", positive = "pink")
#' )
my_plotExpression_flip <- function(sce, marker_list, assay = "logcounts", fill_colors = NULL, title = NULL) {
    message("N nuc: ", ncol(sce))
    genes <- unlist(marker_list)
    # cat_df <- as.data.frame(colData(sce))[,cat, drop = FALSE]
    genes <- genes %in% rownames(sce)
    message("N genes: ", length(genes))
    expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes, ]))

    marker_df <- stack(marker_list)
    names(marker_df) <- c("Var1", "marker")

    expression_long <- dplyr::left_join(expression_long, marker_df, by = "Var1")

    message("creating plot")
    expression_violin <- ggplot(data = expression_long, aes(x = Var1, y = value, fill = marker)) +
        geom_violin(scale = "width") +
        facet_wrap(~marker, scales = "free_x") +
        labs(
            y = paste0("Expression (", assay, ")"),
            title = title
        ) +
        theme_bw() +
        theme(
            legend.position = "None", axis.title.x = element_blank(),
            # axis.text.x=element_text(angle=45,hjust=1),
            strip.text.x = element_text(face = "italic")
        )
    # +
    #   stat_summary(fun = median,
    #                geom = "crossbar",
    #                width = 0.3)
    #
    if (!is.null(fill_colors)) expression_violin <- expression_violin + scale_fill_manual(values = fill_colors)

    # expression_violin
    return(expression_violin)
}

#'
#' my_plotClusterMarkers(example_sce, marker_list, cat = "Mutation_Status", pdf_fn = "cluster_marker_test.pdf")
#'
my_plotClusterMarkers <- function(sce, marker_list, assay = "logcounts", cat = "cellType", fill_colors = NULL, pdf_fn) {
    clusIndexes <- rafalib::splitit(sce[[cat]])
    message(length(clusIndexes))

    message("plotting: ", pdf_fn)
    pdf(pdf_fn, height = 6, width = 8)

    for (i in 1:length(clusIndexes)) {
        message(names(clusIndexes)[[i]])
        sce_temp <- sce[, clusIndexes[[i]]]

        print(
            my_plotExpression_flip(sce_temp,
                marker_list,
                title = paste(names(clusIndexes)[[i]], "n=", ncol(sce_temp)),
                fill_colors = fill_colors
            )
        )
    }

    dev.off()
}

#' my_plotClusterMarkers(example_sce, marker_list, cat = "Mutation_Status", pdf_fn = "cluster_marker_test.pdf")
