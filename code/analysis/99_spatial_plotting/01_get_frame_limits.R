
library("tidyverse")
library("here")

data_dir <-
  here(
    "processed-data", "rdata", "spe","99_spatial_plotting"
  )
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

spaceranger_dirs <- list.files(here("processed-data", "rerun_spaceranger"))
names(spaceranger_dirs) <- map_chr(strsplit(spaceranger_dirs,"_"),~paste0(.x[[2]],"_",.x[[3]]))

# "outs/spatial/tissue_positions_list.csv

sr_csv <- map_chr(spaceranger_dirs, ~here("processed-data", "rerun_spaceranger", .x, "outs", "spatial", "tissue_positions_list.csv"))
all(map_lgl(sr_csv, file.exists))
# TRUE

sr_positions <- map(sr_csv, ~read.csv(.x, header = FALSE))

head(sr_positions[[1]])

frame_lims <- map_dfr(sr_positions, function(p){
  list(x_min = min(p$V5),
       x_max = max(p$V5),
       y_min = min(p$V6),
       y_max = max(p$V6)
  )
})

frame_lims <- frame_lims %>% add_column(sample_id = names(sr_positions), .before = "x_min")

frame_lims
# A tibble: 30 × 5
# sample_id   x_min x_max y_min y_max
# <chr>       <int> <int> <int> <int>
# 1 Br2720_ant  15782 33439  2876 21502
# 2 Br2720_mid   8766 26397  3331 21929
# 3 Br2720_post  8729 26340  3062 21639
# 4 Br2743_ant   6302 24526  7171 26375
# 5 Br2743_mid   8039 26386  4398 23722
# 6 Br2743_post 11158 29455  5242 24520
# 7 Br3942_ant   6108 24308  4812 23991
# 8 Br3942_mid   7818 26073  4009 23246
# 9 Br3942_post 11069 29197  2759 21870
# 10 Br6423_ant   5747 23981  4871 24084

## check x & y are correct
# spe_temp <- spe[,spe$sample_id == "Br2720_ant"]
# # Y columns
# min(spatialCoords(spe_temp)[,"pxl_col_in_fullres"]) #4578
# max(spatialCoords(spe_temp)[,"pxl_col_in_fullres"]) #21245
# # X rows
# min(spatialCoords(spe_temp)[,"pxl_row_in_fullres"]) #15783
# max(spatialCoords(spe_temp)[,"pxl_row_in_fullres"]) #31063

write.csv(frame_lims, file = here(data_dir, "frame_limits.csv"), row.names = FALSE)

## dimensions are all about the same
frame_lims |>
  transmute(x_diff = x_max - x_min,
         y_diff = y_max - y_min) |>
  summary()

# x_diff          y_diff     
# Min.   :17584   Min.   :18550  
# 1st Qu.:17624   1st Qu.:18590  
# Median :17690   Median :18658  
# Mean   :17886   Mean   :18858  
# 3rd Qu.:18232   3rd Qu.:19211  
# Max.   :18347   Max.   :19324  
