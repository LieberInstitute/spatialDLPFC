
library("tidyverse")
library("here")

data_dir <-
  here(
    "processed-data",
    "rdata",
    "spe",
    "14_spatial_registration_PEC"
  )

list.files(data_dir, pattern = "PEC_correlation_annotation")
# [1] "PEC_correlation_annotation_Dx_V5.Rdata" "PEC_correlation_annotation_Dx.Rdata"   
# [3] "PEC_correlation_annotation_V5.Rdata"    "PEC_correlation_annotation.Rdata"

## Dx data
load(here(data_dir, "PEC_correlation_annotation_Dx_V5.Rdata"), verbose = TRUE)
pe_correlation_annotation_dx_v5 <- pe_correlation_annotation

load(here(data_dir, "PEC_correlation_annotation_Dx.Rdata"), verbose = TRUE)
pe_correlation_annotation_dx <- pe_correlation_annotation
  
## All data
load(here(data_dir, "PEC_correlation_annotation_V5.Rdata"), verbose = TRUE)
pe_correlation_annotation_v5 <- pe_correlation_annotation

load(here(data_dir, "PEC_correlation_annotation.Rdata"), verbose = TRUE)
names(pe_correlation_annotation_v5$CMC)

create_anno_long <- function(pe_correlation_annotation){
  layer_anno <- transpose(pe_correlation_annotation)$layer_anno
  names(layer_anno)
  
  layer_anno <- do.call("rbind", layer_anno)
  
  if(!"Dataset" %in% colnames(layer_anno)){
    layer_anno <- layer_anno |>
      rownames_to_column("Dataset") |>
      mutate(Dataset = gsub("\\.[0-9]+", "", Dataset)) 
  }
  
  layer_anno_long <- layer_anno |>
    select(Dataset, cluster, -layer_label, ends_with("label")) |>
    pivot_longer(!c(Dataset, cluster),
                 names_to = "Annotation",
                 values_to = "label"
    ) |>
    mutate(
      confidence = !grepl("\\*", label),
      # layers = str_split(gsub("\\*", "", label), "/"),
      Annotation = gsub("_label", "", Annotation)
    ) 
  
  return(layer_anno_long)
}

## All
anno_long_v5 <- create_anno_long(pe_correlation_annotation_v5) |> 
  rename(label_v5 = label, confidence_v5 = confidence)

anno_long <- create_anno_long(pe_correlation_annotation)

anno_compare_all <- anno_long_v5 |>
  left_join(anno_long) |>
  mutate(match_v5 = label == label_v5)

anno_compare_all |> count(Dataset, match_v5)

anno_compare_all |> filter(!match_v5) |> count(Dataset, Annotation)

anno_diff_all <- anno_compare_all |> 
  filter(!match_v5) |> 
  select(-starts_with("confidence")) |>
  arrange(Dataset, Annotation) 


write.csv(anno_diff_all, file = here(data_dir, "PEC_correlation_annotation_V5_diff.csv"), row.names = FALSE)

## Dx
anno_long_dx_v5 <- create_anno_long(pe_correlation_annotation_dx_v5) |> 
  rename(label_v5 = label, confidence_v5 = confidence)
  
anno_long_dx <- create_anno_long(pe_correlation_annotation_dx)

anno_compare_dx <- anno_long_dx_v5 |>
  left_join(anno_long_dx) |>
  mutate(match_v5 = label == label_v5)
  
anno_compare_dx |> count(Dataset, match_v5)

