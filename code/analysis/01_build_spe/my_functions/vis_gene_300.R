vis_gene_300 <-
    function(spe,
    sampleid,
    geneid = "SCGB2A2; ENSG00000110484",
    spatial = TRUE,
    assayname = "logcounts",
    minCount = 0,
    viridis = TRUE,
    image_id = "lowres",
    alpha = 1,
    cont_colors = if (viridis) viridisLite::viridis(21) else c("aquamarine4", "springgreen", "goldenrod", "red"),
    point_size = 1.25,
    ...) {
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
        p <- vis_gene_p_300(
            spe = spe_sub,
            d = d,
            sampleid = sampleid,
            spatial = spatial,
            title = paste(
                sampleid,
                geneid,
                ...
            ),
            viridis = viridis,
            image_id = image_id,
            alpha = alpha,
            cont_colors = cont_colors,
            point_size = point_size
        )
        p + labs(caption = paste(
            if (!geneid %in% colnames(colData(spe_sub))) {
                assayname
            } else {
                NULL
            },
            "min >", minCount
        ))
    }
