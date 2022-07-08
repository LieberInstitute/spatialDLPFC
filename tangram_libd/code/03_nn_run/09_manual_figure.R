#   The "main tangram deconvolution figure" (Fig 4a from the paper) is hard to
#   generate outside of squidpy, and requires cell centroid info, which doesn't
#   seem critically important to the value of the figure. Here I explore if we
#   can generate a similar figure that only requires cell counts per spot.

library('here')
library('ggplot2')
library('png')
library('spatialLIBD')
library('jaffelab')
library('rjson')
library('grid')

sample_id = '151508'

img_path = here(
    "tangram_libd", "raw-data", "03_nn_run", "hires_histology",
    paste0(sample_id, ".png")
)

clusters_path = here(
    "tangram_libd", "processed-data", "03_nn_run", "tangram_out_DLPFC",
    sample_id, "clusters.csv"
)

scale_path = file.path(
    '/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X',
    sample_id, 'scalefactors_json.json'
)

plot_path = here(
    "tangram_libd", "plots", "03_nn_run", "DLPFC",
    paste0("my_spot_deconvo_fig_", sample_id, ".pdf")
)

spatial_coords_names = c('pxl_row_in_fullres', 'pxl_col_in_fullres')

#   Fetch and subset SPE to this sample
spe = fetch_data(type = 'spe')
spe = spe[,spe$sample_id == sample_id]

#   Read in image, cell counts, and image scale factors
img = readPNG(img_path)
clusters = read.csv(clusters_path)
scale_json = fromJSON(file = scale_path)

#   Add spatial coordinates to 'clusters', the data frame of cell counts per
#   spot
clusters$barcode = ss(clusters$key, '_')
stopifnot(all(clusters$barcode == rownames(spatialCoords(spe))))
clusters = cbind(clusters, spatialCoords(spe))

#   Infer the cell types used
cell_types = colnames(clusters)[
    ! (colnames(clusters) %in% c(
        'key', 'cell_count', 'barcode', spatial_coords_names
        )
    )
]

#   Make a long data frame, where each row is a cell 
df_list = list()
i = 1
for (barcode in rownames(clusters)) {
    for (cell_type in cell_types) {
        for (j in seq_len(clusters[barcode, cell_type])) {
            df_list[[i]] = c(
                barcode,
                as.numeric(clusters[barcode, spatial_coords_names]),
                cell_type
            )
            i = i + 1
        }
    }
}

df_long = data.frame(do.call(rbind, df_list))
colnames(df_long) = c('barcode', spatial_coords_names, 'cell_type')

#   Make sure spatialCoords are numeric, and scaled to represent
#   high-resolution pixels
for (colname in spatial_coords_names) {
    df_long[, colname] =  scale_json$tissue_hires_scalef * 
        as.numeric(df_long[, colname])
}

#   Reverse y coord of spatialCoords(spe) to agree with the coordinate system
#   ggplot2 is using
df_long[[spatial_coords_names[2]]] = dim(img)[2] -
    df_long[[spatial_coords_names[2]]]

#   Improve plotting speed and file size by keeping a copy of just the spots.
#   Get spot radius in high-res pixels
spots = df_long[match(unique(df_long$barcode), df_long$barcode),]
spots$spot_radius = scale_json$spot_diameter_fullres * 
    scale_json$tissue_hires_scalef / 2

#   Create the deconvolution figure
p = ggplot() +
    #   The high-res image as the plot's background
    annotation_custom(
        rasterGrob(
            img, width = unit(1,"npc"), height = unit(1,"npc")
        ),
        -Inf, Inf, -Inf, Inf
    ) +
    #   The axes should exactly span the pixel range of the image, and pixels
    #   should be squares
    scale_x_continuous(limits = c(0, dim(img)[1]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, dim(img)[2]), expand = c(0, 0)) +
    coord_fixed() +
    #   Plot the individual cells, jittered within a spot
    geom_jitter(
        data = df_long,
        aes_string(
            x = spatial_coords_names[1], y = spatial_coords_names[2],
            fill = 'cell_type'
        ),
        size = 0.45, width = 4, height = 4, color = "black", shape = 21,
        stroke = 0.05
    ) +
    #   Plot empty circles for each spot (to show the spot grid)
    geom_circle(
        data = spots,
        mapping = aes_string(
            x0 = spatial_coords_names[1], y0 = spatial_coords_names[2],
            r = 'spot_radius'
        ),
        size = 0.1
    ) +
    #   Brighten colors and improve legend
    scale_fill_hue(l = 80, name = "Cell type") +
    guides(fill = guide_legend(override.aes = list(size = 5)))

pdf(plot_path)
print(p)
dev.off()
