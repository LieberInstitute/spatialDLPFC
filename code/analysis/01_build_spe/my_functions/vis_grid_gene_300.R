vis_grid_gene_300 <-
    function(
        spe,
        geneid = "SCGB2A2; ENSG00000110484",
        pdf_file,
        assayname = "logcounts",
        minCount = 0,
        return_plots = FALSE,
        spatial = TRUE,
        viridis = TRUE,
        height = 24,
        width = 36,
        image_id = "lowres",
        alpha = 1,
        cont_colors = if (viridis) viridisLite::viridis(21) else c("aquamarine4", "springgreen", "goldenrod", "red"),
        sample_order = unique(spe$sample_id),
        point_size = 1.25,
        ...) {
        stopifnot(all(sample_order %in% unique(spe$sample_id)))

        plots <- lapply(sample_order, function(sampleid) {
            vis_gene_300(
                spe,
                sampleid,
                geneid,
                spatial,
                assayname,
                minCount,
                viridis,
                image_id = image_id,
                alpha = alpha,
                cont_colors = cont_colors,
                point_size = point_size,
                ...
            )
        })
        names(plots) <- sample_order

        if (!return_plots) {
            pdf(pdf_file, height = height, width = width)
            print(cowplot::plot_grid(plotlist = plots))
            dev.off()
            return(pdf_file)
        } else {
            return(plots)
        }
    }
