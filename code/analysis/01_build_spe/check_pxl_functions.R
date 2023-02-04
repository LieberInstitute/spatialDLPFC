## Write a function for checking the array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres data
check_pxl <- function(spe_source, spe_target) {
    # stopifnot(all(unique(spe_source$sample_id) %in% unique(spe_target$sample_id)))
    # spe_target$sample_id <- gsub("_2|_all|_extra_reads|DLPFC_|_manual_alignment", "", spe_target$sample_id)
    spe_target <- add_key(spe_target)
    spe_source <- add_key(spe_source)

    m <- match(spe_source$key, spe_target$key)
    if (any(is.na(m))) {
        message("Subsetting to those spots that match in spe_source")
        print(table(is.na(m)))
        spe_source <- spe_source[, !is.na(m)]
        m <- m[!is.na(m)]
    }

    if (identical(spe_source$array_row, spe_target$array_row[m]) &
        identical(spe_source$array_col, spe_target$array_col[m]) &
        identical(spatialCoords(spe_source)[, "pxl_row_in_fullres"], spatialCoords(spe_target)[m, "pxl_row_in_fullres"]) &
        identical(spatialCoords(spe_source)[, "pxl_col_in_fullres"], spatialCoords(spe_target)[m, "pxl_col_in_fullres"])) {
        return("All good!")
    }

    list(
        "array_row" = table(spe_source$array_row == spe_target$array_row[m]),
        "array_col" = table(spe_source$array_col == spe_target$array_col[m]),
        "pxl_row_in_fullres_diff" = summary(spatialCoords(spe_source)[, "pxl_row_in_fullres"] - spatialCoords(spe_target)[m, "pxl_row_in_fullres"]),
        "pxl_col_in_fullres_diff" = summary(spatialCoords(spe_source)[, "pxl_col_in_fullres"] - spatialCoords(spe_target)[m, "pxl_col_in_fullres"]),
        "swap_row_col_diff" = summary(spatialCoords(spe_source)[, "pxl_row_in_fullres"] - spatialCoords(spe_target)[m, "pxl_col_in_fullres"]),
        "swap_col_row_diff" = summary(spatialCoords(spe_source)[, "pxl_col_in_fullres"] - spatialCoords(spe_target)[m, "pxl_row_in_fullres"]),
        "pxl_row_in_fullres_tab" = table(spatialCoords(spe_source)[, "pxl_row_in_fullres"] == spatialCoords(spe_target)[m, "pxl_col_in_fullres"]),
        "pxl_row_in_fullres_tab" = table(spatialCoords(spe_source)[, "pxl_col_in_fullres"] == spatialCoords(spe_target)[m, "pxl_row_in_fullres"])
    )
}

## Check frame limits
check_limits <- function(spe_source, spe_target) {
    x <- list(
        "source" = frame_limits(spe_source, sampleid = "Br6522_ant", image_id = "lowres"),
        "target" = frame_limits(spe_target, sampleid = "Br6522_ant", image_id = "lowres")
    )
    do.call(rbind, lapply(x, as.data.frame))
}

## Make a plot for checking how things are looking
check_plot <- function(spe_source, spe_target, auto_crop = FALSE) {
    spe_target <- add_key(spe_target)
    spe_source <- add_key(spe_source)

    m <- match(spe_source$key, spe_target$key)

    p_source <- vis_gene(
        spe = spe_source,
        sampleid = "Br6522_ant",
        geneid = "expr_chrM_ratio",
        ... = " spe_source",
        auto_crop = auto_crop
    )
    p_target <- vis_gene(
        spe = spe_target[, m],
        sampleid = "Br6522_ant",
        geneid = "expr_chrM_ratio",
        ... = " spe_target",
        auto_crop = auto_crop
    )
    print(cowplot::plot_grid(plotlist = list(p_source, p_target)))
    return(NULL)
}
