library("spatialLIBD")
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <-
  here(
    "plots",
    "99_spatial_plotting",
    "02_test_plots"
  )
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Load full data
load(
  here(
    "processed-data",
    "rdata",
    "spe",
    "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata"
  )
)

rowData(spe)$gene_search <- rowData(spe)$gene_name
spe$bayesSpace_harmony_9 <- as.factor(spe$bayesSpace_harmony_9)
spe$bayesSpace_harmony_16 <- as.factor(spe$bayesSpace_harmony_16)

#### Explore Frame Lims ####
frame_lims <- read.csv(here("processed-data","rdata", "spe","99_spatial_plotting","frame_limits.csv"))
head(frame_lims)

frame_rect_raw <- frame_lims |>
  ggplot(aes(xmin = x_min, xmax = x_max, ymin= y_min, ymax=y_max)) +
  geom_rect(fill = NA, color = "black", alpha = .5) + 
  coord_fixed() 

ggsave(frame_rect_raw, filename = here(plot_dir, "frame_rect_raw.png"))

all_scale_factors <- map_dbl(frame_lims$sample_id, 
                             ~SpatialExperiment::scaleFactors(spe, sample_id = .x, image_id = "lowres"))

frame_lims2 <- frame_lims |>
  add_column(scale_factor = all_scale_factors) |>
  mutate(x_diff = x_max - x_min,
         y_diff = y_max - y_min,
         ratio = x_diff/y_diff,
         area = x_diff * y_diff,
         area_scale = area*scale_factor,
         x_diff_scale = x_diff*scale_factor,
         y_diff_scale = y_diff*scale_factor)

frame_lims2 |> select(sample_id, area, area_scale) |> arrange(area)

frame_lims2 |> select(sample_id, area_scale) |> arrange(area_scale)
#      sample_id area_scale
# 1   Br8325_mid    4872040
# 2   Br2720_ant    4896830
# ...
# 29  Br6432_mid    7668611
# 30 Br6432_post    7681723

frame_rect_diff <- frame_lims2 |>
  ggplot(aes(xmin = 0, xmax = x_diff, ymin= 0, ymax=y_diff, color = ratio)) +
  geom_rect(fill = NA, alpha = .5) + 
  coord_fixed() 

ggsave(frame_rect_diff, filename = here(plot_dir, "frame_rect_diff.png"))

frame_diff_density <- frame_lims2 |>
  select(sample_id, x_diff, y_diff) |>
  pivot_longer(!sample_id, names_to = "dim", values_to ="pixle") |>
  ggplot(aes(x = pixle)) +
  geom_histogram(binwidth=100) +
  facet_wrap(~dim, ncol = 1)

ggsave(frame_diff_density, filename = here(plot_dir, "frame_diff_density.png"))

frame_rect_diff_scale <- frame_lims2 |>
  ggplot(aes(xmin = 0, xmax = x_diff_scale, ymin= 0, ymax=y_diff_scale, color = ratio)) +
  geom_rect(fill = NA, alpha = .5) + 
  coord_fixed() 

ggsave(frame_rect_diff_scale, filename = here(plot_dir, "frame_rect_diff_scale.png"))


## two groups of frames 
# small ~ 18500 x 17600
# large ~ 19200 x 18250

# frame_adj <- c(x_left = -40,
#                x_right = 30, # good  
#                y_up = -35, 
#                y_down = 45)

frame_adj <- c(x_left = -2670,
               x_right = 2070, # good  
               y_down = -2360, 
               y_up = 3022)

frame_edge_lims <- 

frame_area_scatter <- frame_lims2 |>
  ggplot(aes(x = area, y = area_scale, color = x_diff < 19000)) +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 2) +
  geom_point()

ggsave(frame_area_scatter, filename = here(plot_dir, "frame_area_scatter.png"), width = 10)


#### Test limits in plotting regular data ####
samples <- c("Br8325_ant","Br2720_ant","Br8667_mid","Br6522_ant","Br2743_mid","Br6432_post","Br8667_post")
# samples <- c("Br8325_mid","Br6432_post") # min and max scale areas

map(samples,~SpatialExperiment::scaleFactors(spe, sample_id = .x, image_id = "lowres"))

frame_lims2 |> filter(sample_id %in% samples)

walk(samples, function(samp){
  fl <- frame_lims[match(samp, frame_lims$sample_id),-1]
  scale_factor <- SpatialExperiment::scaleFactors(spe, sample_id = samp, image_id = "lowres")
  fls <- fl*scale_factor
  vis_gene_test <- vis_gene(
    spe = spe,
    point_size = 1.2,
    sampleid = samp,
    geneid = "PCP4"
  ) + 
    coord_fixed() 
  
  fls_plus <- fls + (frame_adj*scale_factor)
  
  vis_gene_test_lims <- vis_gene_test+
    geom_vline(xintercept = c(fls$x_min, fls$x_max), color = "red", linetype = "dashed") +
    geom_hline(yintercept = c(fls$y_min, fls$y_max), color = "red", linetype = "dashed") +
    ## outside frame
    geom_vline(xintercept = fls_plus$x_min, color = "blue", linetype = "dashed") +
    geom_vline(xintercept = fls_plus$x_max, color = "green", linetype = "dashed") +
    geom_hline(yintercept = fls_plus$y_min, color = "purple", linetype = "dashed") +
    geom_hline(yintercept = fls_plus$y_max , color = "goldenrod", linetype = "dashed")
  
  ggsave(vis_gene_test_lims, filename = here(plot_dir, paste0("vis_gene_ggsave_",samp,".png")))
  
})

# ggsave(vis_gene_test, filename = here(plot_dir, "vis_gene_ggsave.pdf"))

# geom_rect(xmin = fls$x_min, xmax = fls$x_max,  ## Doesn't work??
#           ymin = fls$y_min, ymax = fls$y_max,
#           color = "red", linetype = "dashed"
#           ) 



img <- SpatialExperiment::imgRaster(spe, sample_id = "Br2720_ant", image_id = "lowres")
dim(img)

# sample_id x_min x_max y_min y_max
# 1   Br2720_ant 15782 33439  2876 21502
sample <- "Br2720_ant"
frame_lims[match(sample,frame_lims$sample_id),]


#### Test vis crop functions ####
source("vis_gene_crop.R")

vis_gene_custom_test <- vis_gene_crop(
  spe = spe,
  point_size = 2.2,
  frame_lim_df = frame_lims,
  sampleid = "Br6432_ant",
  geneid = "PCP4"
) 

ggsave(vis_gene_custom_test, filename = here(plot_dir, "vis_gene_crop.png"))

walk(samples, function(samp){
  
  vis_gene_crop_plot <- vis_gene_crop(
    spe = spe,
    point_size = 2.2,
    frame_lim_df = frame_lims,
    sampleid = samp,
    geneid = "PCP4"
  )  
  ggsave(vis_gene_crop_plot, filename = here(plot_dir, "vis_gen_crop", paste0("vis_gene_crop_",samp,".png")))
  
})

#### Crop to fudcial frame ####
frame_adj <- list(x_left = -2350,
               x_right = 1950, # good  
               y_down = -2500, # Add space for legend 
               y_up = 2650)


frame_edge_lims <- frame_lims |>
  mutate(x_min = x_min + frame_adj$x_left,
         x_max = x_max + frame_adj$x_right,
         y_min = y_min + frame_adj$y_down,
         y_max = y_max + frame_adj$y_up)

vis_gene_crop_frame <- vis_gene_crop(
  spe = spe,
  point_size = 2,
  frame_lim_df = frame_edge_lims,
  sampleid = "Br6432_ant",
  geneid = "PCP4",
  legend_overlap = TRUE
) 

ggsave(vis_gene_crop_frame, filename = here(plot_dir, "vis_gene_crop_frame.png"))
ggsave(vis_gene_crop_frame, filename = here(plot_dir, "vis_gene_crop_frame.pdf"))

walk(frame_edge_lims$sample_id, function(samp){
  # message(samp)
  vis_gene_crop_plot <- vis_gene_crop(
    spe = spe,
    point_size = 2.2,
    frame_lim_df = frame_edge_lims,
    sampleid = samp,
    geneid = "PCP4",
    legend_overlap = TRUE
  )  
  ggsave(vis_gene_crop_plot, filename = here(plot_dir, "vis_gen_crop_frame", paste0("vis_gene_crop_frame_",samp,".png")))
  
})


