#   The "main tangram deconvolution figure" (Fig 4a from the paper) is hard to
#   generate outside of squidpy, and requires cell centroid info, which doesn't
#   seem critically important to the value of the figure. Here I generate a
#   similar figure that only requires cell counts per spot.

library('here')
library('ggplot2')
library('ggforce')
library('grid')
library('png')
library('jaffelab')
library('rjson')
library('SpatialExperiment')

sample_id_path = here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", "sample_ids.txt"
)
sample_id = readLines(sample_id_path)[as.numeric(Sys.getenv('SGE_TASK_ID'))]
print(paste0('Plotting sample ', sample_id, '.'))

img_path = here(
    'processed-data', '01_spaceranger_IF', sample_id, 'outs', 'spatial',
    'tissue_hires_image.png'
)

clusters_path = here(
    "processed-data", "spot_deconvo", "01-tangram", "IF", sample_id,
    'clusters.csv'
)

scale_path = here(
    'processed-data', '01_spaceranger_IF', sample_id, 'outs', 'spatial',
    'scalefactors_json.json'
)

spe_in = here(
    "processed-data", "rdata", "spe_IF", "01_build_spe_IF", "spe.rds"
)

plot_path = here(
    "plots", "spot_deconvo", "01-tangram", "IF", sample_id,
    "my_spot_deconvo_fig.pdf"
)

spatial_coords_names = c('pxl_row_in_fullres', 'pxl_col_in_fullres')

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)

#   Load and subset SpatialExperiment
spe = readRDS(spe_in)
spe = spe[,spe$sample_id == sample_id]

#   Read in image, cell counts, and image scale factors
img = readPNG(img_path)
img = aperm(img, c(2, 1, 3))
clusters = read.csv(clusters_path)
scale_json = fromJSON(file = scale_path)

clusters$barcode = ss(clusters$key, '_', 1)
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
df_long[[spatial_coords_names[2]]] = dim(img)[1] -
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
    scale_x_continuous(limits = c(0, dim(img)[2]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, dim(img)[1]), expand = c(0, 0)) +
    coord_fixed() +
    #   Plot the individual cells, jittered within a spot
    geom_jitter(
        data = df_long,
        aes_string(
            x = spatial_coords_names[1], y = spatial_coords_names[2],
            fill = 'cell_type'
        ),
        size = 0.25, width = 4, height = 4, color = "black", shape = 21,
        stroke = 0.05
    ) +
    #   Plot empty circles for each spot (to show the spot grid)
    geom_circle(
        data = spots,
        mapping = aes_string(
            x0 = spatial_coords_names[1], y0 = spatial_coords_names[2],
            r = 'spot_radius'
        ),
        size = 0.05, color = "white"
    ) +
    #   Brighten colors and improve legend
    scale_fill_hue(l = 80, name = "Cell type") +
    guides(fill = guide_legend(override.aes = list(size = 5)))

pdf(plot_path)
print(p)
dev.off()
