# Load Packages
library(tidyverse)
library(here)
library(SpatialExperiment)
library(ggpubr)
library(SingleCellExperiment)

# Load Spe Object
spe_dat <- readRDS(
    here(
        "processed-data/rdata/spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
    )
)

# Config Color ------------------------------------------------------------
sce <- readRDS("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC_annotated/se.rds")
cell_type_colors_layer <- metadata(sce)$cell_type_colors_layer[levels(sce$cellType_layer)]
cell_type_colors_broad <- metadata(sce)$cell_type_colors_broad[levels(sce$cellType_broad_hc)]


factor_cell_type_broad <- function(vec){
    factor(vec,
           levels = c("astro",  "endomural", "micro", "oligo", "opc", "excit", "inhib"
           ),
           labels = c("Astro",  "EndoMural", "Micro", "Oligo", "OPC", "Excit", "Inhib" )
    )
}

factor_cell_type_layer <- function(vec){
    factor(vec,
           # levels = c("Astro",  "EndoMural", "Micro", "Oligo", "OPC",
           #            "Excit_L2_3", "Excit_L3", "Excit_L3_4_5", "Excit_L4","Excit_L5",
           #            "Excit_L5_6","Excit_L6","Inhib"),
           levels = c("astro",  "endomural", "micro", "oligo", "opc",
                      "excit_l2_3", "excit_l3", "excit_l3_4_5", "excit_l4","excit_l5",
                      "excit_l5_6","excit_l6","inhib"),
           labels = c("Astro",  "EndoMural", "Micro", "Oligo", "OPC",
                      "Excit_L2/3", "Excit_L3", "Excit_L3/4/5", "Excit_L4","Excit_L5",
                      "Excit_L5/6","Excit_L6","Inhib")

    )
}


spe <- spe_dat[,spe_dat$sample_id == "Br8667_mid"]

deconv_comb <- expand_grid(
    res = c("broad", "layer"),
    deconv = c("tangram", "cell2location","spotlight")
)

res <- "layer"
deconv <- "cell2location"

spatial <- FALSE
auto_crop <- TRUE
image_id <- "lowres"
sampleid <- "Br8667_mid"


factor_cell_type <- factor_cell_type_layer
cell_type_colors <- cell_type_colors_layer
fill_vals <- c("astro",  "endomural", "micro", "oligo", "opc",
               "excit_l2_3", "excit_l3", "excit_l3_4_5", "excit_l4","excit_l5",
               "excit_l5_6","excit_l6","inhib")
fill_labs <-c("Astro",  "EndoMural", "Micro", "Oligo", "OPC",
              "Excit_L2/3", "Excit_L3", "Excit_L3/4/5", "Excit_L4","Excit_L5",
              "Excit_L5/6","Excit_L6","Inhib")
cell_type_colors <- set_names(cell_type_colors[fill_labs],
                              fill_vals)

if(res == "broad"){
    factor_cell_type <- factor_cell_type_broad
    cell_type_colors <- cell_type_colors_broad
    fill_vals <- c("astro",  "endomural", "micro", "oligo", "opc", "excit", "inhib")
    fill_labs <-c("Astro",  "EndoMural", "Micro", "Oligo", "OPC", "Excit", "Inhib" )
    cell_type_colors <- set_names(cell_type_colors[fill_labs],
                                  fill_vals)
}

plot_df <-  colData(spe) |> data.frame() |>
    dplyr::select(
        starts_with(paste0(res, "_", deconv)), #cell_names
    ) |>
        rename_with(.fn = str_remove, pattern = paste0(res, "_", deconv, "_"))

cell_names <- names(plot_df)

top_three_mat <- plot_df |>
    apply(MARGIN = 1,
          FUN = function(x) head(order(x, decreasing = TRUE),3),
          simplify = FALSE) |>
    map_dfr(
        .f  = function(top_3_ind){
            ret_vec <- rep(0, length(cell_names)) |>
                               set_names(cell_names)
            ret_vec[top_3_ind] <- 1

            return(ret_vec)
        }
    )

stopifnot(all(rowSums(top_three_mat) == 3)) # Error prevention



top_six_mat <- plot_df |>
    apply(MARGIN = 1,
          FUN = function(x) head(order(x, decreasing = TRUE),6),
          simplify = FALSE) |>
    map_dfr(
        .f  = function(top_ind){
            ret_vec <- rep(0, length(cell_names)) |>
                set_names(cell_names)
            ret_vec[top_ind] <- 1

            return(ret_vec)
        }
    )

stopifnot(all(rowSums(top_six_mat) == 6)) # Error prevention



# add spatial coordinate

plot_df_fnl <- top_three_mat |>
    cbind(
        colData(spe) |> data.frame() |>
            dplyr::select(
                row, col  # Spatial Coordinates
            )
    )



all_cell_names <- grep(
    pattern = paste0(res, "_", deconv),
    x = names(colData(spe)),
    value = TRUE) |>
    str_remove(pattern = paste0(res, "_", deconv, "_"))




# plot_df <- plot_df |> select(-colap_cells)

# browser()
ret_plot <- vis_scatter_pie(spe,
                            plot_df = plot_df_fnl,
                            sampleid = sampleid,
                            image_id = image_id,
                            cell_names = all_cell_names,
                            spatial = FALSE)


ggsave(here("plots/15_cell_composition/",
            "scatter_pie",
            paste0("scatter_pie_",
                   # res, "_",deconv,"_Br8667_mid_spatial.pdf")),
                   res, "_",deconv,"_Br8667_mid_top3.pdf")),
       plot = ret_plot,
       height = 7.5,
       width = 9
)



# Top 6 -------------------------------------------------------------------

top6_plot_df_fnl <- top_six_mat |>
    cbind(
        colData(spe) |> data.frame() |>
            dplyr::select(
                row, col  # Spatial Coordinates
            )
    )


# plot_df <- plot_df |> select(-colap_cells)

# browser()
top_6_ret_plot <- vis_scatter_pie(spe,
                            plot_df = top6_plot_df_fnl,
                            sampleid = sampleid,
                            image_id = image_id,
                            cell_names = all_cell_names,
                            spatial = FALSE)


ggsave(here("plots/15_cell_composition/",
            "scatter_pie",
            paste0("scatter_pie_",
                   # res, "_",deconv,"_Br8667_mid_spatial.pdf")),
                   res, "_",deconv,"_Br8667_mid_top6.pdf")),
       plot = top_6_ret_plot,
       height = 7.5,
       width = 9
)


# Coloc Spots -------------------------------------------------------------
coloc_spe <- calc_coloc(spe, "EFNA5", "EPHA5", sample_id = NULL)

coloc_spe$coloc=="co-localize"

coloc_ret_plot <- vis_scatter_pie(spe[,coloc_spe$coloc=="co-localize"],
                            plot_df = plot_df_fnl[coloc_spe$coloc=="co-localize",],
                            sampleid = sampleid,
                            image_id = image_id,
                            cell_names = all_cell_names,
                            auto_crop = FALSE,
                            spatial = FALSE)

ggsave(here("plots/15_cell_composition/",
            "scatter_pie",
            paste0("scatter_pie_",
                   # res, "_",deconv,"_Br8667_mid_spatial.pdf")),
                   res, "_",deconv,"_Br8667_mid_top3_coloc.pdf")),
       plot = coloc_ret_plot,
       height = 7.5,
       width = 9
)


