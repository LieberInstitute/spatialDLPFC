
library("spatialLIBD")
library("tidyverse")
library("xlsx")
library("jaffelab")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

## plot dir
plot_dir <- here("plots", "14_spatial_registration_PEC", "06_PEC_correlation_annotation_Dx")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <-
    here(
        "processed-data",
        "rdata",
        "spe",
        "14_spatial_registration_PEC"
    )

#### Load Layer and k Data  ####
layer_modeling_results <- fetch_data(type = "modeling_results")

paths <-
    list(k09 = "modeling_results_BayesSpace_k09.Rdata", k16 = "modeling_results_BayesSpace_k16.Rdata")

modeling_results <- lapply(paths, function(x) {
    get(load(
        here(
            "processed-data",
            "rdata",
            "spe",
            "07_layer_differential_expression",
            x
        )
    ))
})

modeling_results <-
    c(list(layer = layer_modeling_results), modeling_results)
names(modeling_results)

pec_cell_types <- c(
  # Non-neuronal cells (8)
  "Astro",
  "Endo",
  "Immune",
  "Micro",
  "OPC",
  "Oligo",
  "PC",
  "SMC",      
  #Excit (9)
  "L2/3 IT",
  "L4 IT",
  "L5 ET",
  "L5 IT",
  "L5/6 NP",
  "L6 CT",
  "L6 IT",
  "L6 IT Car3",
  "L6b",
  # Inhib (10)
  "Chandelier",
  "Lamp5",
  "Lamp5 Lhx6",
  "Pax6",
  "Pvalb",
  "Sncg",
  "Sst",
  "Sst Chodl",
  "VLMC",
  "Vip"
)

pec_cell_type_tb <- tibble(cell_type = factor(pec_cell_types, levels = pec_cell_types), 
                          cluster =  make.names(cell_type),
                          ct_cat = unlist(map2(c("Non-neuronal", "Excit", "Inhib"),c(8, 9, 10), rep)))

# layer_anno[[1]] |> left_join(pec_cell_type_tb) |> filter(is.na(cell_type))

correlate_and_annotate <- function(dataset, make_cor_plot = FALSE) {
    #### Load Registration Stats ####
    message(Sys.time(), " - ", dataset)
    registration_stats <- readRDS(here(
        "processed-data",
        "rdata",
        "spe",
        "14_spatial_registration_PEC",
        paste0("registration_stats_Dx_", dataset, ".rds")
    ))
    
    cor_top100 <- map(registration_stats, function(dx_stats){
      
      registration_t_stats <-
        dx_stats[, grep("^t_stat", colnames(dx_stats))]
      
      colnames(registration_t_stats) <-
        gsub("^t_stat_", "", colnames(registration_t_stats))
      
      # Fix Cell Types

      #### Correlate with modeling results ####
      cor_top100 <-
        map(
          modeling_results,
          ~ layer_stat_cor(
            registration_t_stats,
            .x,
            model_type = "enrichment",
            reverse = FALSE,
            top_n = 100
          )
        )
      
      # rownames(cor_top100) <-
      #   cell_types[rownames(cor_top100)]
      
      return(cor_top100)
    })


    #### Annotate Layers ####
    layer_anno <- map2(cor_top100, names(cor_top100),function(dx_cor, dx){
      
      layer_anno <- map2(dx_cor, names(dx_cor), function(cor, name) {
        anno <- annotate_registered_clusters(
          cor_stats_layer = cor,
          confidence_threshold = 0.25,
          cutoff_merge_ratio = 0.10
        )
        colnames(anno) <- gsub("layer", name, colnames(anno))
        return(anno)
      })
      
      layer_anno <- reduce(layer_anno, left_join, by = "cluster")
      layer_anno$PrimaryDx <- dx
      
      return(layer_anno)
    })
    
    layer_anno2 <- do.call("bind_rows", layer_anno)
    # layer_anno$Dataset <- dataset
    
    return(list(cor_top100 = cor_top100, layer_anno = layer_anno2))
}
        
# datasets <- c("CMC", "DevBrain-snRNAseq", "IsoHuB", "SZBDMulti", "UCLA-ASD", "Urban-DLPFC")
datasets <-
  c( # No isohub not. case-control
    "CMC",
    "DevBrain-snRNAseq",
    "UCLA-ASD",
    "MultiomeBrain-DLPFC",
    "SZBDMulti-Seq"
  )
names(datasets) <- datasets

## Calculate correlations and annotations for each dataset
pe_correlation_annotation <- map(datasets, correlate_and_annotate)

save(pe_correlation_annotation,
    file = here(data_dir, "PEC_correlation_annotation_Dx.Rdata")
)

# load(here(data_dir, "PEC_correlation_annotation_Dx.Rdata"), verbose = TRUE)

#### Save Output to XLSX sheet ####
key <- data.frame(
    data = c("annotation", paste0("cor_", names(modeling_results))),
    description = c(
        "Annotations of spatial registration",
        "Correlation values with manual layer annotations",
        "Correlation values with k09 domains",
        "Correlation values with k16 domains"
    )
)

## Clear file and write key
annotation_xlsx <- here(data_dir, "PEC_spatial_annotations.xlsx")
write.xlsx(
    key,
    file = annotation_xlsx,
    sheetName = "Key",
    append = FALSE,
    row.names = FALSE
)

## write annotations
walk2(
    pe_correlation_annotation,
    names(pe_correlation_annotation),
    ~ write.xlsx(
        .x$layer_anno,
        file = annotation_xlsx,
        sheetName = paste0("annotation_", .y),
        append = TRUE,
        row.names = FALSE
    )
)

## write correlations
walk2(pe_correlation_annotation, names(pe_correlation_annotation), function(data, name) {
    name <- paste0("cor_", name)
    # message(name)
    walk2(
        data$cor_top100,
        names(data$cor_top100),
        ~ write.xlsx(
            t(.x),
            file = annotation_xlsx,
            sheetName = paste0(name, "_", .y),
            append = TRUE,
            row.names = TRUE
        )
    )
})

#### Compare annotations for each cell type ####
## prep data
source(
    here("code", "analysis", "12_spatial_registration_sn", "utils.R"),
    echo = TRUE,
    max.deparse.length = 500
)

## Match with bayesSpace spatial annotations for k09 k16 exploration
bayes_layers <- get(load(here("processed-data", "rdata", "spe", "08_spatial_registration", "bayesSpace_layer_annotations.Rdata"))) |>
  select(Annotation = bayesSpace, layers = cluster, layer_combo) |>
  filter(Annotation %in% c("k09", "k16"))

## Add Layer levels + update factor
spatial_layer_anno <- tibble(
  Annotation = "layer",
  layers = c(paste0("Layer", 1:6), "WM"),
  layer_combo = layers
) |>
  bind_rows(bayes_layers) |>
  mutate(layer_combo = factor(layer_combo, levels = c(paste0("Layer", 1:6), "WM", levels(bayes_layers$layer_combo))),
         Annotation = factor(Annotation, levels = c("k16", "k09", "layer")))

## Extract Layer anno, pivot each longer 
layer_anno <- transpose(pe_correlation_annotation)$layer_anno

layer_anno_long <- map(layer_anno,  ~.x|> 
    left_join(pec_cell_type_tb,  by = "cluster") |>
    select(PrimaryDx, cell_type, ct_cat, ends_with("label")) |>
    pivot_longer(!c(PrimaryDx, cell_type, ct_cat),
        names_to = "Annotation",
        values_to = "label"
    ) |>
    mutate(
        confidence = !grepl("\\*", label),
        layers = str_split(gsub("\\*", "", label), "/"),
        Annotation = gsub("_label", "", Annotation)
    ) |>
    unnest_longer("layers") |>
    # mutate(layers = ifelse(Annotation == "layer" & grepl("^[0-9]",layers), paste0("L",layers), layers))
    mutate(
      layers = ifelse(Annotation == "layer",
            ifelse(
                grepl("^[0-9]", layers),
                paste0("Layer", layers),
                gsub("L","Layer",layers)
            ),
            layers
        ),
        # layers_long = gsub("Layers","L",layers),
      layer_confidence = ifelse(confidence, "good", "poor"),
      Annotation = factor(Annotation, levels = c("k16", "k09", "layer"))) |> 
      left_join(spatial_layer_anno)
    )


sort(unique(unlist(map(layer_anno_long, ~.x$PrimaryDx))))
# [1] "Autism Spectrum Disorder" "Bipolar"                  "Bipolar Disorder"         "Control"                 
# [5] "Schizophrenia"            "Williams Syndrome"     

dx_colors = c(
  Control = "#2e282a",
  `Autism Spectrum Disorder` = "#76b041",
  `Bipolar` = "#17bebb",
  `Bipolar Disorder` = "#17bebb",
  Schizophrenia = "#e4572e",
  `Williams Syndrome` = "#CC9C00"
)

#### Plot layer annotation ####
## prep cell type annotation
spatial_layer_anno_long <- spatial_layer_anno |> 
  mutate(ml = str_split(gsub("S.* ~ |ayer","", layer_combo), "/")) |> ## ml = "manual layer"
  unnest_longer("ml") |> 
  mutate(ml = ifelse(grepl("^[0-9]", ml), paste0("L", ml), ml))

cell_type_anno <- pec_cell_type_tb |>
  mutate(layer_label = ifelse(grepl("^L[0-9]", cell_type), sub(" .*", "", cell_type), NA)) |>
  filter(!is.na(layer_label)) |>
  mutate(ml = str_split(gsub("\\*", "", layer_label), "/")) |> ## ml = "manual layer"
  unnest_longer("ml") |>
  mutate(
    ml = gsub("b", "", ifelse(
      grepl("^[0-9]", ml), paste0("L", ml), ml)),
    Match = TRUE) |>
  left_join(spatial_layer_anno_long) |>
  select(cell_type, layer_combo, Match)


spatial_layer_anno |> 
  mutate(ml = str_split(gsub("S.* ~ |ayer","", layer_combo), "/")) |> ## ml = "manual layer"
  unnest_longer("ml")

source("registration_dot_plot.R")

## Plot Dx
walk2(layer_anno_long, names(layer_anno_long), function(data, name){
  
  temp_dx_colors <- dx_colors[unique(data$PrimaryDx)]
  
  dotplot <- registration_dot_plot2(data, ct_anno = cell_type_anno) +
    labs(title = name) +
    scale_color_manual(values = temp_dx_colors)
  
    ggsave(dotplot, filename = here(plot_dir, paste0("registration_anno_dotplot_dx_",name,".png")), width = 10)

    # dotplot_conf <- registration_dot_plot2(data, ct_anno = cell_type_anno, conf_only = TRUE) + 
    #   labs(title = name) +
    #   scale_color_manual(values = temp_dx_colors)
    # ggsave(dotplot_conf, filename = here(plot_dir, paste0("registration_anno_dotplot_dx_conf_",name,".png")), width = 10)
    # 
})

#### Combine Control, plot to compare Datasets ####

## IsoHuB??
control_anno_long <- do.call("rbind", map2(layer_anno_long, names(layer_anno_long),
                                           ~.x |> 
                                             mutate(Dataset = .y) |> 
                                             filter(PrimaryDx == "Control")
                                           )
                             )
   
dotplot_control_datasets <- registration_dot_plot2(control_anno_long, color_by = "Dataset", ct_anno = cell_type_anno) +
  theme(legend.position = "bottom")
ggsave(dotplot_control_datasets, filename = here(plot_dir, "registration_anno_dotplot_control.png"), width = 11)
ggsave(dotplot_control_datasets, filename = here(plot_dir, "registration_anno_dotplot_control.pdf"), width = 11)

#### ASD Astrocyte Subset ####

astro_anno <- do.call("rbind", map2(layer_anno_long[c("DevBrain-snRNAseq", "UCLA-ASD")], c("DevBrain", "UCLA-ASD"), 
     ~.x |> 
       filter(cell_type == "Astro") |>
       mutate(cell_type = paste0(cell_type, "\n", .y))))

asd_dx_colors <- dx_colors[unique(astro_anno$PrimaryDx)]

astro_asd_datasets <- registration_dot_plot2(astro_anno, color_by = "PrimaryDx") +
  scale_color_manual(values = asd_dx_colors) +
  # theme(axis.text.x = element_text(angle = 0))
  theme_bw()

ggsave(astro_asd_datasets, filename = here(plot_dir, "registration_anno_dotplot_astro_asd.png"), height = 5)
ggsave(astro_asd_datasets, filename = here(plot_dir, "registration_anno_dotplot_astro_asd.pdf"), height = 5)

#### ASD Astrocyte Heatmap ####




# sgejobs::job_single('03_correlate_spatial', create_shell = TRUE, memory = '25G', command = "Rscript 03_correlate_spatial.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
