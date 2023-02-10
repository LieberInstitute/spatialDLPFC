
# Plot --------------------------------------------------------------------

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

# deconv_comb |>
#     head(1) |>
#     pwalk(.f=function(res, deconv){

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
        row, col  # Spatial Coordinates
    ) |>
    rename_with(.fn = str_remove, pattern = paste0(res, "_", deconv, "_"))
all_cell_names <- grep(
    pattern = paste0(res, "_", deconv),
    x = names(colData(spe)),
    value = TRUE) |>
    str_remove(pattern = paste0(res, "_", deconv, "_"))




# plot_df <- plot_df |> select(-colap_cells)

# browser()
ret_plot <- vis_scatter_pie(spe,
                            plot_df = plot_df,
                            sampleid = sampleid,
                            image_id = image_id,
                            cell_names = all_cell_names,
                            spatial = FALSE)


ggsave(here("plots/15_cell_composition/",
            "scatter_pie",
            paste0("scatter_pie_",
                   # res, "_",deconv,"_Br8667_mid_spatial.pdf")),
                   res, "_",deconv,"_Br8667_mid.pdf")),
       plot = ret_plot,
       height = 7.5,
       width = 9
)
