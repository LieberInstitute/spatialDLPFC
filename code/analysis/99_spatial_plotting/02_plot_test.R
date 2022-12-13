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

frame_lims2 <- frame_lims |>
  mutate(x_diff = x_max - x_min,
         y_diff = y_max - y_min,
         ratio = x_diff/y_diff,
         area = x_diff * y_diff
  )

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


frame_area <- 
  
  ## two groups of frames 
  # small ~ 18500 x 17600
  # large ~ 19200 x 18250
  
  # frame_adj <- c(x_left = -40,
  #                x_right = 30, # good  
  #                y_up = -35, 
  #                y_down = 45)
  
  frame_adj <- c(x_left = -2670,
                 x_right = 2070, # good  
                 y_up = -2360, 
                 y_down = 3022)

#### Test limits in plotting regular data ####
samples <- c("Br8325_ant","Br2720_ant","Br8667_mid","Br6522_ant","Br2743_mid")

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

vis_gene_test_lim <- vis_gene(
  spe = spe,
  point_size = 1.2,
  sampleid = "Br2720_ant",
  geneid = "CLDN5"
) + 
  coord_fixed() +
  xlim(sl$x_min*0.0148894, sl$x_max*0.0148894) +
  ylim(sl$y_min*0.0148894, sl$y_max*0.0148894)

ggsave(vis_gene_test_lim, filename = here(plot_dir, "vis_gene_lim_ggsave.png"))

vis_gene_test_lim_nohist <- vis_gene(
  spe = spe,
  point_size = 1.2,
  sampleid = "Br2720_ant",
  geneid = "CLDN5",
  spatial = FALSE
) + 
  coord_fixed()

ggsave(vis_gene_test_lim_nohist, filename = here(plot_dir, "vis_gene_test_lim_nohist_ggsave.png"))



img <- SpatialExperiment::imgRaster(spe, sample_id = "Br2720_ant", image_id = "lowres")
dim(img)

# sample_id x_min x_max y_min y_max
# 1   Br2720_ant 15782 33439  2876 21502
sample <- "Br2720_ant"
frame_lims[match(sample,frame_lims$sample_id),]


## try custom functions
source("vis_gene_crop.R")

vis_gene_custom_test <- vis_gene_crop(
  spe = spe,
  point_size = 2,
  frame_lim_df = frame_lims,
  sampleid = "Br6432_ant",
  geneid = "CLDN5"
) 

ggsave(vis_gene_custom_test, filename = here(plot_dir, "vis_gene_custom.png"))

