# Function to calculate the gene_name1, gene_name2 co-localization
# of an SPE object
# Return a SPE object with co-localization status as one variable in the colData
calc_coloc <- function(spe, gene_name1, gene_name2, sample_id = NULL){

    if(!is.null(sample_id)){
        if(!sample_id %in% unique(spe$sample_id))
            stop("Can't find sample_id")
        spe <- spe[,spe$sample_id == sample_id]
    }

    if(!all(c(gene_name1, gene_name2) %in% rowData(spe)$gene_name))
        stop("Can't find gene names")

    gene1_idx <- which(rowData(spe)$gene_name == gene_name1)
    gene2_idx <- which(rowData(spe)$gene_name == gene_name2)

    cnt_df <- counts(spe)[c(gene1_idx, gene2_idx),] |>
        as.matrix() |> t() |> data.frame()

    stopifnot(setequal(names(cnt_df),rowData(spe)[names(cnt_df),"gene_id"]))

    old_names <- names(cnt_df)
    names(cnt_df) <- paste0("gene", 1:2)
    # Check if the correct rows are selected

    coloc_cat <- cnt_df |>
        pmap_chr(.f = function(gene1, gene2){
            if(gene1!=0 && gene2!=0) return("co-localize")
            else if(gene1 != 0) return(gene_name1)
            else if(gene2 != 0) return(gene_name2)
            else return("Neither")
        })

    colData(spe)$coloc <- coloc_cat |> factor()

    return(spe)
}


# Function to plot spatial co-localization derived from `calc_coloc` in the
# Spatial domain context, specifically `Sp9D7 ~ L6`
# Return a spot plot oject if `save.path=NULL`, otherwise, save pdf to `save.path`
vis_coloc <- function(spe, gene_name1, gene_name2, sample_id,
                      save.path = NULL){

    ret_plot <- spatialLIBD::vis_clus(
        spe, clustervar = "coloc", point_size = 2.5,
        sampleid = sample_id, spatial = FALSE,
        colors = c("#d7191c", "#abd9e9", "#2c7bb6", "#CCCCCC40") |>
            set_names(c("co-localize", gene_name1, gene_name2, "Neither"))
    ) + labs(title = "") +
        scale_fill_manual(
            # breaks = c(TRUE),
            # labels = c("Sp09D07~L6"),
            breaks = c("co-localize", gene_name1, gene_name2),
            values = c("#d7191c", "#abd9e9", "#2c7bb6") |>
                set_names(c("co-localize", gene_name1, gene_name2)),
            labels = c(paste0(gene_1, " & ", gene_2),
                       paste0(gene_1, " only"),
                       paste0(gene_2, " only")),
            na.value = "#CCCCCC40"

        )

    tmp <- ret_plot +
        geom_point(aes(color = spe$BayesSpace_harmony_09==7),
                   shape = 21,
                   fill = "transparent",
                   size = 2.5,
                   stroke = 0.15) +
        scale_color_manual(
            breaks = c(TRUE),
            labels = c("Sp09D07~L6"),
            values = c("#757575", "transparent"),
            na.value = "transparent"

        ) +
        labs(color = NULL,
             fill = NULL)


    if(!is.null(save.path)){
        ggsave(
            paste0(
                save.path, "/",
                paste0(
                    paste(sample_id, gene_name1, gene_name2, sep  = "-"),
                    ".pdf")
            ),
            plot = tmp,
            height = 8,
            width = 9
        )

    }

    return(ret_plot)
}


# Plot all SPD ------------------------------------------------------------
vars_spd <- c(
    "BayesSpace_harmony_09",
    "BayesSpace_harmony_16"
) |>
    set_names()

# Factor the Spd Labels ----
bayes_layers <- here(
    "processed-data", "rdata", "spe",
    "08_spatial_registration",
    "bayesSpace_layer_annotations.Rdata"
) |>
    load() |>
    get() |>
    select(Annotation = bayesSpace, layer_long = cluster, layer_combo) |>
    filter(Annotation %in% c("k09", "k16"))


# Function to plot spatial co-localization derived from `calc_coloc` in the
# all spatial domains of `Sp9D7`
# Return a spot plot object if `save.path=NULL`, otherwise, save pdf to `save.path`
vis_coloc_spd <- function(spe, gene_name1, gene_name2, sample_id,
                          save.path = NULL){

    ret_plot <- spatialLIBD::vis_clus(
        spe, clustervar = "coloc", point_size = 2.5,
        sampleid = sample_id, spatial = FALSE,
        colors = c("#fc8d62", "#66c2a5", "#8da0cb", "#CCCCCC40") |>
            set_names(c("co-localize", gene_name1, gene_name2, "Neither"))
    ) + labs(title = "")

    layer_df <- bayes_layers |> filter(Annotation == "k09")
    SpD_fct <- factor(spe$BayesSpace_harmony_09,
                      levels = str_sub(layer_df$layer_long, -2, -1) |>
                          as.integer(),
                      labels = layer_df$layer_combo)

    tmp <- ret_plot +
        geom_point(aes(shape = spe$coloc), size = 1) +
        geom_point(aes(color = SpD_fct),
                   shape = 21,
                   fill = "transparent",
                   size = 2.5,
                   stroke = 0.1) +
        scale_color_manual(
            values = set_names(
                Polychrome::palette36.colors(9)[str_sub(layer_df$layer_long, -2, -1) |>
                                                    as.integer()],
                layer_df$layer_combo
            ),
            na.value = "transparent"
        ) +
        scale_shape_manual(
            # breaks = c(TRUE),
            # labels = c("Sp09D07~L6"),
            breaks = c("co-localize", gene_name1, gene_name2
            ),
            values = c(8,NA,NA # 8,3,4
            )|>
                set_names(c("co-localize", gene_name1, gene_name2
                )),
            na.value = NA

        ) +
        scale_fill_manual(
            breaks = c("co-localize", gene_name1, gene_name2),
            values = c("#fc8d62", "#66c2a5", "#8da0cb") |>
                set_names(c("co-localize", gene_name1, gene_name2)),
            na.value = "#CCCCCC40"

        ) +
        labs(color = NULL,
             shape = NULL)



    if(!is.null(save.path)){
        ggsave(
            paste0(
                save.path, "/",
                paste0(
                    paste(sample_id, gene_name1, gene_name2, "all_Spd",sep  = "-"),
                    ".pdf")
            ),
            plot = tmp,
            height = 8,
            width = 9
        )

    }

    return(ret_plot)
}
