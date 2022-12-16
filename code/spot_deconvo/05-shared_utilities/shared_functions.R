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
#
#   Returns a ggplot object
spot_plot = function(
        spe, sample_id, title, var_name, include_legend, is_discrete,
        colors = NULL, assayname = "logcounts"
        ) {
    
    POINT_SIZE = 1.85
    
    #   If the quantity to plot is discrete, use 'vis_clus'. Otherwise use
    #   'vis_gene'.
    if (is_discrete) {
        #   For 'vis_clus' only, supply a color scale if 'color' is not NULL
        if (is.null(colors)) {
            p <- vis_clus(
                spe, sampleid = sample_id, clustervar = var_name,
                return_plots = TRUE, spatial = FALSE, point_size = POINT_SIZE
            )
        } else {
            p <- vis_clus(
                spe, sampleid = sample_id, clustervar = var_name,
                return_plots = TRUE, spatial = FALSE, colors = colors,
                point_size = POINT_SIZE
            )
        }
    } else {
        p <- vis_gene(
            spe, sampleid = sample_id, geneid = var_name, return_plots = TRUE,
            spatial = FALSE, point_size = POINT_SIZE, assayname = assayname
        )
    }
    
    #   Remove the legend if requested
    if (! include_legend) {
        p = p + theme(legend.position = "none")
    }
    
    #   Use fixed coordinate scale and remove the caption
    p = p + labs(title = title, caption = NULL)
    
    return(p)
}
