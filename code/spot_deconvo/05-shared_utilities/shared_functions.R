################################################################################
#   Functions related to spot plots
################################################################################

#   Wrapper around 'vis_gene' and 'vis_clus' to ensure exactly consistent
#   aspect ratio and point sizes, for use in producing manuscript figures
#
#   spe:            passed to 'spe' in either 'vis_gene' or 'vis_clus'
#   sample_id:      passed to 'sampleid'
#   title:          title for the plot, expected to be one line (avoid use of
#                   "\n")
#   var_name:       passed to 'geneid' for 'vis_gene' and 'clustervar' for
#                   'vis_clus'
#   include_legend: (logical) if FALSE, remove the legend
#   is_discrete:    (logical) if TRUE, use 'vis_clus'; otherwise, use 'vis_gene'
#   colors:         passed to 'colors' for 'vis_gene' if not [is_discrete]
#   assayname:      passed to 'assayname' for 'vis_gene' if not [is_discrete]
#   minCount:       passed to 'minCount' for 'vis_gene' if not [is_discrete]
#
#   Returns a ggplot object
spot_plot <- function(spe, sample_id, title, var_name, include_legend, is_discrete,
    colors = NULL, assayname = "logcounts", minCount = 0.5) {
    POINT_SIZE <- 2.3

    #   If the quantity to plot is discrete, use 'vis_clus'. Otherwise use
    #   'vis_gene'.
    if (is_discrete) {
        #   For 'vis_clus' only, supply a color scale if 'color' is not NULL
        if (is.null(colors)) {
            p <- vis_clus(
                spe,
                sampleid = sample_id, clustervar = var_name,
                return_plots = TRUE, spatial = FALSE, point_size = POINT_SIZE
            )
        } else {
            p <- vis_clus(
                spe,
                sampleid = sample_id, clustervar = var_name,
                return_plots = TRUE, spatial = FALSE, colors = colors,
                point_size = POINT_SIZE
            )
        }
    } else {
        p <- vis_gene(
            spe,
            sampleid = sample_id, geneid = var_name, return_plots = TRUE,
            spatial = FALSE, point_size = POINT_SIZE, assayname = assayname,
            cont_colors = viridisLite::plasma(21), alpha = 1,
            minCount = minCount
        )
    }

    #   Remove the legend if requested
    if (!include_legend) {
        p <- p + theme(legend.position = "none")
    }

    #   Overwrite the title
    p <- p + labs(title = title)

    return(p)
}

#   Set the largest value for the color/fill scale to [upper_limit] and return
#   a copy of the plot 'p' with that modification. 'min_count' should be the
#   value of 'minCount' passed to 'spot_plot', used to create 'p'.
overwrite_scale <- function(p, upper_limit, min_count) {
    p +
        scale_fill_gradientn(
            colors = viridisLite::plasma(21),
            limits = c(0, upper_limit), na.value = c("#CCCCCC40")
        ) +
        scale_color_gradientn(
            colors = viridisLite::plasma(21),
            limits = c(0, upper_limit), na.value = c("#CCCCCC40")
        ) +
        labs(fill = paste(" min >", min_count))
}

#   Given a list of ggplot objects 'plot_list', write up to 3 PDFs with
#   variations of this list. Write the default version as a grid with [ncol]
#   columns, expected to have a legend and a 1-line title; this is designed for
#   best internal viewing. Write a second version with legends but no title, for
#   use in grid-based plots in the manuscript. If [include_individual], write a
#   third version with one plot per page, the existing titles, and no legends.
#
#   plot_list:          list of ggplot objects to plot to multiple PDFs
#   n_col:              (integer) number of columns in grid versions of the
#                       plots
#   plot_dir:           (character) path to base directory for writing plots
#   file_prefix:        (character) start of filename (without extension)
#   include_individual: (logical) if TRUE, write a 3rd version of plots with one
#                       plot per page
#
#   Returns NULL.
write_spot_plots <- function(plot_list, n_col, plot_dir, file_prefix, include_individual) {
    #   Scaling factor for 'plot_grid'. 1.03 seems to reduce whitespace without
    #   introducing overlap, but larger values quickly introduce overlap
    SCALE <- 1.03

    #   Create two alternate versions of the default plot for use in the
    #   manuscript. One version has one plot per page, with consistent
    #   point size and aspect ratio with other individual spot plots. The
    #   other is a grid version with a legend but no title, consistent with
    #   grid plots in the manuscript
    plot_list_individual <- list()
    plot_list_no_title <- list()
    for (i in 1:length(plot_list)) {
        plot_list_individual[[i]] <- plot_list[[i]] +
            theme(legend.position = "none")
        plot_list_no_title[[i]] <- plot_list[[i]] + labs(title = NULL)
    }

    #   Internal-viewing version
    pdf(
        file.path(plot_dir, paste0(file_prefix, ".pdf")),
        width = 7 * n_col,
        height = 7 * length(plot_list) / n_col
    )
    print(plot_grid(plotlist = plot_list, ncol = n_col), scale = SCALE)
    dev.off()

    #   Individual version, if requested
    if (include_individual) {
        pdf(file.path(plot_dir, paste0(file_prefix, "_individual.pdf")))
        print(plot_list_individual)
        dev.off()
    }

    #   No-title version
    pdf(
        file.path(
            plot_dir, paste0(file_prefix, "_no_title.pdf")
        ),
        width = 7 * n_col,
        height = 7 * length(plot_list) / n_col
    )
    print(plot_grid(plotlist = plot_list_no_title, ncol = n_col, scale = SCALE))
    dev.off()
}

################################################################################
#   Other functions
################################################################################

#   Make a list of which layers we expect each cell type to be most highly
#   expressed in. Necessary for the below function 'layer_dist_barplot'
corresponding_layers <- list(
    "Astro" = "L1",
    "EndoMural" = "L1",
    "Excit" = paste0("L", 2:6),
    "Excit_L2_3" = c("L2", "L3"),
    "Excit_L3" = "L3",
    "Excit_L3_4_5" = c("L3", "L4", "L5"),
    "Excit_L4" = "L4",
    "Excit_L5" = "L5",
    "Excit_L5_6" = c("L5", "L6"),
    "Excit_L6" = "L6",
    "Inhib" = paste0("L", 2:6),
    "Micro" = c("L1", "WM"),
    "Oligo" = "WM",
    "OPC" = c("L1", "WM")
)

#   Given a tibble with columns 'label' (manual layer label), 'deconvo_tool',
#   'cell_type', and 'count', write a set of barplots to PDF under [plot_dir]
#   with name [filename]. 'ylab' give the y-axis label; 'x_var' is the x-axis
#   variable (as a string); 'fill_var' is the fill variable as a string;
#   'fill_scale' is passed to 'scale_fill_manual(values = [fill_scale])';
#   'fill_lab' is the fill label; 'xlab' is the x-axis label
#
#   The barplots are faceted by deconvo_tool, with x-axis including each
#   manually annotated layer. Each barplot includes counts for each cell type
#   in each layer. Each cell type is expected to have a maximal value (across
#   all bars in the facet) at a particular layer; an "O" is placed at the layer
#   with maximal value for each cell type if the layer is "correct" (e.g.
#   'Excit_L3' has maximal value in layer 3), and an "X" is placed if the layer
#   is incorrect. Total counts of "O"s for each facet across all cell types is
#   tallied and reported in the facet titles.
layer_dist_barplot <- function(counts_df, filename, ylab, x_var, fill_var, fill_scale, fill_lab, xlab) {
    #   Add a column 'layer_match' to indicate rows where each cell type
    #   has a maximal value across layers. We'll mark these with an "X"
    #   on the barplots
    counts_df <- counts_df |>
        group_by(deconvo_tool, cell_type) |>
        mutate(layer_match = count == max(count)) |>
        ungroup()

    #   Add a column 'correct_layer' indicating whether for a cell type
    #   and deconvo tool, the cell_type has maximal value in the correct/
    #   expected layer
    counts_df$correct_layer <- sapply(
        1:nrow(counts_df),
        function(i) {
            counts_df$layer_match[i] &&
                (counts_df$label[i] %in%
                    corresponding_layers[[
                        as.character(counts_df$cell_type)[i]
                    ]]
                )
        }
    )

    #   For each deconvo tool, add up how many times cell types have maximal
    #   value in the correct layers
    correct_df <- counts_df |>
        group_by(deconvo_tool) |>
        summarize(num_matches = sum(correct_layer)) |>
        ungroup()
    print("Number of times cell types have maximal value in the correct layer:")
    print(correct_df)

    print("Full list of which cell types matched the expected layer, by method:")
    counts_df |>
        group_by(deconvo_tool) |>
        filter(correct_layer) |>
        select(cell_type) |>
        ungroup() |>
        print(n = nrow(counts_df))

    #   Add the "layer accuracy" in the facet titles in the upcoming plot
    correct_labeller <- paste0(
        correct_df$deconvo_tool, ": ", correct_df$num_matches, "/",
        length(cell_types)
    )
    names(correct_labeller) <- correct_df |> pull(deconvo_tool)
    correct_labeller <- labeller(deconvo_tool = correct_labeller)

    p <- ggplot(
        counts_df,
        aes_string(x = x_var, y = "count", fill = fill_var)
    ) +
        facet_wrap(~deconvo_tool, labeller = correct_labeller) +
        geom_bar(stat = "identity") +
        labs(x = xlab, y = ylab, fill = fill_lab) +
        scale_fill_manual(values = fill_scale) +
        geom_text(
            aes(label = ifelse(correct_layer, "O", ifelse(layer_match, "X", ""))),
            position = position_stack(vjust = 0.5)
        ) +
        theme_bw(base_size = 16) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

    pdf(file.path(plot_dir, filename), width = 10, height = 5)
    print(p)
    dev.off()
}
