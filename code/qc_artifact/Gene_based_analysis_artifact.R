rm(list = ls())
library(ggplot2)
library(Seurat)
library(SeuratObject)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(cowplot)

# Read manual annotationss
Br6522_ant_layers <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br6522_ant_layers.csv")
Br6522_mid_layers <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br6522_mid_layers.csv")
Br8667_post_layers <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br8667_post_layers.csv")

Br6522_ant_wrinkles <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br6522_ant_wrinkle.csv")
Br6522_mid_wrinkles <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br6522_mid_wrinkle.csv")
Br8667_post_wrinkles <- read.csv("/data/abattle4/prashanthi/dewrinkler/manual_annotations/spatialLIBD_ManualAnnotation_Br8667_post_wrinkle.csv")

scDir <- "/home/pravich2/scratch16-abattle4/prashanthi/dewrinkler/processed-data/NextSeq/Round3/"
Br6522_ant_spot_counts <- read.csv(paste0(scDir, "DLPFC_Br6522_ant_manual_alignment_all/outs/spatial/tissue_spot_counts_new.csv"))
Br6522_mid_spot_counts <- read.csv(paste0(scDir, "DLPFC_Br6522_mid_manual_alignment_all/outs/spatial/tissue_spot_counts_new.csv"))
Br8667_post_spot_counts <- read.csv(paste0(scDir, "DLPFC_Br8667_post_manual_alignment_all/outs/spatial/tissue_spot_counts_new.csv"))

datDir <- "/data/abattle4/prashanthi/dewrinkler/data/wrinkles_manual_thresholding/"
Br6522_ant <- readRDS(paste0(datDir, "Br6522_ant.rds"))
Br6522_mid <- readRDS(paste0(datDir, "Br6522_mid.rds"))
Br8667_post <- readRDS(paste0(datDir, "Br8667_post.rds"))

# Add manual metadata
manually_annotate <- function(sobj, layers_df, wrinkles_df, sc_df){
  spots <- rownames(sobj@meta.data)
  layers <- c()
  wrinkle <- c()
  nCells <- c()
  for(i in c(1:length(spots))){
    if(spots[i] %in% layers_df$spot_name){
      layers[i] <- layers_df$ManualAnnotation[layers_df$spot_name == spots[i]]
    }else{
      layers[i] <- "Unknown"
    }
    if(spots[i] %in% wrinkles_df$spot_name){
      wrinkle[i] <- wrinkles_df$ManualAnnotation[wrinkles_df$spot_name == spots[i]]
    }else{
      wrinkle[i] <- "None"
    }
    nCells[i] <- sc_df$Nmask_dark_blue[sc_df$barcode == spots[i]]
  }
  sobj@meta.data["Layers"] <- layers
  sobj@meta.data["Wrinkles"] <- wrinkle
  sobj@meta.data["nCells"] <- nCells
  sobj
}

Br6522_ant <- manually_annotate(Br6522_ant, Br6522_ant_layers, Br6522_ant_wrinkles, Br6522_ant_spot_counts)
Br6522_mid <- manually_annotate(Br6522_mid, Br6522_mid_layers, Br6522_mid_wrinkles, Br6522_mid_spot_counts)
Br8667_post <- manually_annotate(Br8667_post, Br8667_post_layers, Br8667_post_wrinkles, Br8667_post_spot_counts)

Br6522_ant <- subset(Br6522_ant, Layers!="Unknown")
Br6522_mid <- subset(Br6522_mid, Layers!="Unknown")
Br8667_post <- subset(Br8667_post, Layers!="Unknown")

Br6522_ant$Layers <- factor(Br6522_ant$Layers, levels = c("Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "WM"))
Br6522_mid$Layers <- factor(Br6522_mid$Layers, levels = c("Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "WM"))
Br8667_post$Layers <- factor(Br8667_post$Layers, levels = c("Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6"))

Br6522_ant$Wrinkles <- gsub("_", " ", Br6522_ant$Wrinkles)
Br6522_mid$Wrinkles <- gsub("_", " ", Br6522_mid$Wrinkles)
Br8667_post$Wrinkles <- gsub("_", " ", Br8667_post$Wrinkles)

Br6522_ant$Wrinkles <- factor(Br6522_ant$Wrinkles, levels = c("None", "Fold 1", "Shear 1", "Shear 2", "Shear 3", "Wrinkle 1", "Wrinkle 2", 
                                                            "Wrinkle 3", "Wrinkle 4", "Wrinkle 5", "Wrinkle 6", "Wrinkle 7", "Wrinkle 8", 
                                                            "Wrinkle 9"))
Br6522_mid$Wrinkles <- factor(Br6522_mid$Wrinkles, levels = c("None", "Fold 1", "Shear 1", "Shear 2", "Wrinkle 1", "Wrinkle 2", 
                                                            "Wrinkle 3", "Wrinkle 4", "Wrinkle 5", "Wrinkle 6", "Wrinkle 7", "Wrinkle 8", 
                                                            "Wrinkle 9", "Wrinkle 10", "Wrinkle 11", "Wrinkle 12", "Wrinkle 13"))
Br8667_post$Wrinkles <- factor(Br8667_post$Wrinkles, levels = c("None", "Shear 1", "Wrinkle 1", "Wrinkle 2", 
                                                              "Wrinkle 3", "Wrinkle 4", "Wrinkle 5", "Wrinkle 6", "Wrinkle 7", "Wrinkle 8", 
                                                              "Wrinkle 9", "Wrinkle 11", "Wrinkle 12", "Wrinkle 13"))

Br6522_ant_1 <- Br6522_ant
Br6522_mid_1 <- Br6522_mid
Br8667_post_1 <- Br8667_post

Br6522_ant_1$Wrinkles <- as.character(Br6522_ant_1$Wrinkles)
Br6522_ant_1$Wrinkles[Br6522_ant_1$Wrinkles == "None"] <- NA
Br6522_ant_1$Wrinkles <- factor(Br6522_ant_1$Wrinkles, c("Fold 1", "Shear 1", "Shear 2", "Shear 3", "Wrinkle 1", "Wrinkle 2", 
                                                         "Wrinkle 3", "Wrinkle 4", "Wrinkle 5", "Wrinkle 6", "Wrinkle 7", "Wrinkle 8", 
                                                         "Wrinkle 9"))

Br6522_mid_1$Wrinkles <- as.character(Br6522_mid_1$Wrinkles)
Br6522_mid_1$Wrinkles[Br6522_mid_1$Wrinkles == "None"] <- NA
Br6522_mid_1$Wrinkles <- factor(Br6522_mid_1$Wrinkles, c("Fold 1", "Shear 1", "Shear 2", "Wrinkle 1", "Wrinkle 2", 
                                                         "Wrinkle 3", "Wrinkle 4", "Wrinkle 5", "Wrinkle 6", "Wrinkle 7", "Wrinkle 8", 
                                                         "Wrinkle 9", "Wrinkle 10", "Wrinkle 11", "Wrinkle 12", "Wrinkle 13"))

Br8667_post_1$Wrinkles <- as.character(Br8667_post_1$Wrinkles)
Br8667_post_1$Wrinkles[Br8667_post_1$Wrinkles == "None"] <- NA
Br8667_post_1$Wrinkles <- factor(Br8667_post_1$Wrinkles, c("Shear 1", "Wrinkle 1", "Wrinkle 2", 
                                                           "Wrinkle 3", "Wrinkle 4", "Wrinkle 5", "Wrinkle 6", "Wrinkle 7", "Wrinkle 8", 
                                                           "Wrinkle 9", "Wrinkle 11", "Wrinkle 12", "Wrinkle 13"))

layer_palette <- c("Layer 1" = "#F0027F", "Layer 2" = "#377EB8", "Layer 3" = "#4DAF4A", 
                   "Layer 4" = "#984EA3", "Layer 5" = "#FFD700", "Layer 6" = "#FF7F00", 
                   "WM" = "#1A1A1A", "NA" = "transparent")
wrinkle_palette <- c("Fold 1" = "#e6194B", "Shear 1" = "#3cb44b", 
                     "Shear 2" = "#4363d8", "Shear 3" = "#f58231", 
                     "Wrinkle 1" = "#911eb4", "Wrinkle 2" = "#42d4f4", 
                     "Wrinkle 3" = "#f032e6", "Wrinkle 4" = "#bfef45",
                     "Wrinkle 5" = "#469990", "Wrinkle 6" = "#9A6324",
                     "Wrinkle 7" = "#800000", "Wrinkle 8" = "#aaffc3", 
                     "Wrinkle 9" = "#373e02", "Wrinkle 10" = "#000075",
                     "Wrinkle 11" = "#ffd8b1" , "Wrinkle 12" = "#fffac8", 
                     "Wrinkle 13" = "#ffe119", "NA" = "transparent")

make_spatial_plots <- function(sobj, title){
  Idents(sobj) <- sobj$Layers
  p1 <- SpatialDimPlot(sobj, cols = layer_palette) + theme_bw() + guides(fill=guide_legend(title="")) + xlab("") + ylab("") +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(face = "bold", size = 10), legend.text = element_text(size = 7))
  Idents(sobj) <- sobj$Wrinkles
  p2 <- SpatialDimPlot(sobj, cols = wrinkle_palette)+ theme_bw() + guides(fill=guide_legend(title="", ncol = 2)) + xlab("") + ylab("") +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 7))
  list(p1, p2)
}

sp_plot_1 <- make_spatial_plots(Br6522_ant_1, "Br6522 Anterior")
sp_plot_2 <- make_spatial_plots(Br6522_mid_1, "Br6522 Middle")
sp_plot_3 <- make_spatial_plots(Br8667_post_1, "Br8667 Posterior")

pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br6522_ant_layers.pdf")
sp_plot_1[[1]]
dev.off()

pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br6522_mid_layers.pdf")
sp_plot_2[[1]]
dev.off()

pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br8667_post_layers.pdf")
sp_plot_3[[1]]
dev.off()

pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br6522_ant_artifacts.pdf")
sp_plot_1[[2]]
dev.off()

pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br6522_mid_artifacts.pdf")
sp_plot_2[[2]]
dev.off()

pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br8667_post_artifacts.pdf")
sp_plot_3[[2]]
dev.off()

Br6522_ant$percent_reads_LG <- apply(Br6522_ant@assays$Spatial@counts, 2, max)/ apply(Br6522_ant@assays$Spatial@counts, 2, sum)
Br6522_mid$percent_reads_LG <- apply(Br6522_mid@assays$Spatial@counts, 2, max)/ apply(Br6522_mid@assays$Spatial@counts, 2, sum)
Br8667_post$percent_reads_LG <- apply(Br8667_post@assays$Spatial@counts, 2, max)/ apply(Br8667_post@assays$Spatial@counts, 2, sum)

Br6522_ant$percent_reads_LG <- Br6522_ant$percent_reads_LG*100
Br6522_mid$percent_reads_LG <- Br6522_mid$percent_reads_LG*100
Br8667_post$percent_reads_LG <- Br8667_post$percent_reads_LG*100

Br6522_ant$is_wrinkle <- !Br6522_ant$Wrinkles == "None"
Br6522_mid$is_wrinkle <- !Br6522_mid$Wrinkles == "None"
Br8667_post$is_wrinkle <- !Br8667_post$Wrinkles == "None"

make_metadata_layers_plots <- function(sobj, title){
  Layers <- levels(sobj$Layers)
  nSpots_normal <- c()
  nSpots_artifacts <- c()
  for(i in c(1:length(Layers))){
    if(sum(sobj$Layers == Layers[i]) == 0){
      nSpots_normal[i] <- 0
      nSpots_artifacts[i] <- 0
    }else{
      nSpots_normal[i] <- sum(!sobj$is_wrinkle[sobj$Layers == Layers[i]])
      nSpots_artifacts[i] <- sum(sobj$is_wrinkle[sobj$Layers == Layers[i]])
    }
  }
  count_df_artifact <- data.frame("Layers" = Layers, "nSpots" = nSpots_artifacts, 
                                "Artifact" = TRUE)
  count_df_normal <- data.frame("Layers" = Layers, "nSpots" = nSpots_normal, 
                                "Artifact" = FALSE)
  count_df <- rbind(count_df_artifact, count_df_normal)
  p0 <- ggplot(count_df, aes(x = Layers, y = nSpots, fill = Artifact)) + geom_bar(position="dodge", stat="identity") + theme_classic() +
        ylab("Number of spots") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title)
  
  stat.test1 <- sobj@meta.data %>%
    group_by(Layers) %>%
    wilcox_test(nCount_Spatial ~ is_wrinkle, alternative = "less") %>%
    add_significance("p")
  stat.test1 <- stat.test1 %>%
    add_xy_position(x = "Layers", dodge = 0.8)
  p1 <- ggplot(sobj@meta.data) + geom_boxplot( aes(x = Layers, y = nCount_Spatial, fill = is_wrinkle),outlier.size = 0.4) + theme_classic() + 
    ylab("Library Size") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title) + 
    stat_pvalue_manual(
      stat.test1,  label = "{p.signif}", 
      tip.length = 0, hide.ns = TRUE)
  
  stat.test2 <- sobj@meta.data %>%
    group_by(Layers) %>%
    wilcox_test(percent.mt ~ is_wrinkle, alternative = "greater") %>%
    add_significance("p")
  stat.test2 <- stat.test2 %>%
    add_xy_position(x = "Layers", dodge = 0.8)
  p2 <- ggplot(sobj@meta.data) + geom_boxplot(aes(x = Layers, y = percent.mt , fill = is_wrinkle), outlier.size = 0.4) + theme_classic() + 
    ylab("Percentage mitochondrial reads") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title) + 
    stat_pvalue_manual(
      stat.test2,  label = "{p.signif}", 
      tip.length = 0, hide.ns = TRUE)
  
  stat.test3 <- sobj@meta.data %>%
    group_by(Layers) %>%
    wilcox_test(nFeature_Spatial ~ is_wrinkle, alternative = "less") %>%
    add_significance("p")
  stat.test3 <- stat.test3 %>%
    add_xy_position(x = "Layers", dodge = 0.8)
  p3 <- ggplot(sobj@meta.data) + geom_boxplot(aes(x = Layers, y = nFeature_Spatial , fill = is_wrinkle), outlier.size = 0.4) + theme_classic() + 
    ylab("Detected genes") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title) +
  stat_pvalue_manual(
    stat.test3,  label = "{p.signif}", 
    tip.length = 0, hide.ns = TRUE)
  
  stat.test4 <- sobj@meta.data %>%
    group_by(Layers) %>%
    wilcox_test(nCells ~ is_wrinkle, alternative = "less") %>%
    add_significance("p")
  stat.test4 <- stat.test4 %>%
    add_xy_position(x = "Layers", dodge = 0.8)
  p4 <- ggplot(sobj@meta.data) + geom_boxplot(aes(x = Layers, y = nCells , fill = is_wrinkle), outlier.size = 0.4) + theme_classic() + 
    ylab("Nuclei detected (DAPI)") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title) +
    stat_pvalue_manual(
      stat.test4,  label = "{p.signif}", 
      tip.length = 0, hide.ns = TRUE)
  
  stat.test5 <- sobj@meta.data %>%
    group_by(Layers) %>%
    wilcox_test(percent_reads_LG ~ is_wrinkle, alternative = "greater") %>%
    add_significance("p")
  stat.test5 <- stat.test5 %>%
    add_xy_position(x = "Layers", dodge = 0.8)
  p5 <- ggplot(sobj@meta.data) + geom_boxplot(aes(x = Layers, y = percent_reads_LG , fill = is_wrinkle), outlier.size = 0.4) + theme_classic() + 
    ylab("Percent reads (Largest gene)") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title) +
    stat_pvalue_manual(
      stat.test5,  label = "{p.signif}", 
      tip.length = 0, hide.ns = TRUE)
  
  list(p0, p1, p2, p3, p4, p5)
}

md_plot_1 <- make_metadata_layers_plots(Br6522_ant, "Br6522 Anterior")
md_plot_2 <- make_metadata_layers_plots(Br6522_mid, "Br6522 Middle")
md_plot_3 <- make_metadata_layers_plots(Br8667_post, "Br8667 Posterior")

md_plot_1[[1]] + md_plot_2[[1]] + md_plot_3[[1]]
md_plot_1[[2]] + md_plot_2[[2]] + md_plot_3[[2]]
md_plot_1[[3]] + md_plot_2[[3]] + md_plot_3[[3]]
md_plot_1[[4]] + md_plot_2[[4]] + md_plot_3[[4]]
md_plot_1[[5]] + md_plot_2[[5]] + md_plot_3[[5]]
md_plot_1[[6]] + md_plot_2[[6]] + md_plot_3[[6]]

make_metadata_wrinkles_plots <- function(sobj, title, layer){
  Wrinkles <- levels(sobj$Wrinkles)
  nSpots_L1 <- c()
  nSpots_L2 <- c()
  nSpots_L3 <- c()
  nSpots_L4 <- c()
  nSpots_L5 <- c()
  nSpots_L6 <- c()
  nSpots_WM <- c()
  for(i in c(1:length(Wrinkles))){
      nSpots_L1[i] <- sum(sobj$Layers[sobj$Wrinkles == Wrinkles[i]] == "Layer 1")
      nSpots_L2[i] <- sum(sobj$Layers[sobj$Wrinkles == Wrinkles[i]] == "Layer 2")
      nSpots_L3[i] <- sum(sobj$Layers[sobj$Wrinkles == Wrinkles[i]] == "Layer 3")
      nSpots_L4[i] <- sum(sobj$Layers[sobj$Wrinkles == Wrinkles[i]] == "Layer 4")
      nSpots_L5[i] <- sum(sobj$Layers[sobj$Wrinkles == Wrinkles[i]] == "Layer 5")
      nSpots_L6[i] <- sum(sobj$Layers[sobj$Wrinkles == Wrinkles[i]] == "Layer 6")
      nSpots_WM[i] <- sum(sobj$Layers[sobj$Wrinkles == Wrinkles[i]] == "Layer WM")
  }
  count_df_L1 <- data.frame("Wrinkles" = Wrinkles, "nSpots" = nSpots_L1, 
                                  "Layer" = "Layer 1")
  count_df_L2 <- data.frame("Wrinkles" = Wrinkles, "nSpots" = nSpots_L2, 
                           "Layer" = "Layer 2")
  count_df_L3 <- data.frame("Wrinkles" = Wrinkles, "nSpots" = nSpots_L3, 
                           "Layer" = "Layer 3")
  count_df_L4 <- data.frame("Wrinkles" = Wrinkles, "nSpots" = nSpots_L4, 
                           "Layer" = "Layer 4")
  count_df_L5 <- data.frame("Wrinkles" = Wrinkles, "nSpots" = nSpots_L5, 
                           "Layer" = "Layer 5")
  count_df_L6 <- data.frame("Wrinkles" = Wrinkles, "nSpots" = nSpots_L1, 
                           "Layer" = "Layer 1")
  count_df_WM <- data.frame("Wrinkles" = Wrinkles, "nSpots" = nSpots_L6, 
                            "Layer" = "Layer 6")
  
  count_df <- rbind(count_df_L1, count_df_L2, count_df_L3,
                    count_df_L4, count_df_L5, count_df_L6, count_df_WM)
  count_df <- count_df[!count_df$Wrinkles == "None", ]
  p0 <- ggplot(count_df, aes(Wrinkles, nSpots, fill = Layer)) + geom_bar(position="dodge", stat="identity") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(title) + scale_fill_brewer(palette = "Dark2") + xlab("")
  
  stat.test1 <- sobj@meta.data %>%
    group_by(Layers) %>%
    wilcox_test(nCount_Spatial ~ Wrinkles, alternative = "less", ref.group = "None") %>%
    add_significance("p")
  stat.test1 <- stat.test1 %>%
    add_xy_position(x = "Layers", dodge = 0.8)
  p1 <- ggplot(sobj@meta.data) + geom_boxplot(aes(x = Layers, y = nCount_Spatial, fill = Wrinkles),outlier.size = 0.4) + theme_classic() + 
    ylab("Library Size") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title) + 
    stat_pvalue_manual(
      stat.test1,  label = "{p.signif}", 
      tip.length = 0, hide.ns = TRUE) 
  
  stat.test2 <- sobj@meta.data %>%
    group_by(Layers) %>%
    wilcox_test(percent.mt ~ Wrinkles, alternative = "greater", ref.group = "None") %>%
    add_significance("p")
  stat.test2 <- stat.test2 %>%
    add_xy_position(x = "Layers", dodge = 0.8)
  p2 <- ggplot(sobj@meta.data) + geom_boxplot(aes(x = Layers, y = percent.mt , fill = Wrinkles), outlier.size = 0.4) + theme_classic() + 
    ylab("Percentage mitochondrial reads") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title) + 
    stat_pvalue_manual(
      stat.test2,  label = "{p.signif}", 
      tip.length = 0, hide.ns = TRUE)
  
  stat.test3 <- sobj@meta.data %>%
    group_by(Layers) %>%
    wilcox_test(nFeature_Spatial ~ Wrinkles, alternative = "less", ref.group = "None") %>%
    add_significance("p")
  stat.test3 <- stat.test3 %>%
    add_xy_position(x = "Layers", dodge = 0.8)
  p3 <- ggplot(sobj@meta.data) + geom_boxplot(aes(x = Layers, y = nFeature_Spatial , fill = Wrinkles), outlier.size = 0.4) + theme_classic() + 
    ylab("Detected genes") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title) +
    stat_pvalue_manual(
      stat.test3,  label = "{p.signif}", 
      tip.length = 0, hide.ns = TRUE)
  
  stat.test4 <- sobj@meta.data %>%
    group_by(Layers) %>%
    wilcox_test(nCells ~ Wrinkles, alternative = "less", ref.group = "None") %>%
    add_significance("p")
  stat.test4 <- stat.test4 %>%
    add_xy_position(x = "Layers", dodge = 0.8)
  p4 <- ggplot(sobj@meta.data) + geom_boxplot(aes(x = Layers, y = nCells , fill = Wrinkles), outlier.size = 0.4) + theme_classic() + 
    ylab("Nuclei detected (DAPI)") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title) +
    stat_pvalue_manual(
      stat.test4,  label = "{p.signif}", 
      tip.length = 0, hide.ns = TRUE)
  
  stat.test5 <- sobj@meta.data %>%
    group_by(Layers) %>%
    wilcox_test(percent_reads_LG ~ Wrinkles, alternative = "greater", ref.group = "None") %>%
    add_significance("p")
  stat.test5 <- stat.test5 %>%
    add_xy_position(x = "Layers", dodge = 0.8)
  p5 <- ggplot(sobj@meta.data) + geom_boxplot(aes(x = Layers, y = percent_reads_LG , fill = Wrinkles), outlier.size = 0.4) + theme_classic() + 
    ylab("Percent reads (Largest gene)") + guides(fill=guide_legend(title="Artifact")) + ggtitle(title) +
    stat_pvalue_manual(
      stat.test5,  label = "{p.signif}", 
      tip.length = 0, hide.ns = TRUE)
  
  list(p0, p1, p2, p3, p4, p5)
}

mdw_plot_1 <- make_metadata_wrinkles_plots(Br6522_ant, "Br6522 Anterior")
mdw_plot_2 <- make_metadata_wrinkles_plots(Br6522_mid, "Br6522 Middle")
mdw_plot_3 <- make_metadata_wrinkles_plots(Br8667_post, "Br8667 Posterior")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mdw_legend_1 <- get_legend(mdw_plot_1[[2]])
mdw_legend_2 <- get_legend(mdw_plot_2[[2]])
mdw_legend_3 <- get_legend(mdw_plot_3[[2]])

mdw_plot_1[[2]] <- mdw_plot_1[[2]] + theme(legend.position="none")
mdw_plot_1[[3]] <- mdw_plot_1[[3]] + theme(legend.position="none")
mdw_plot_1[[4]] <- mdw_plot_1[[4]] + theme(legend.position="none")
mdw_plot_1[[5]] <- mdw_plot_1[[5]] + theme(legend.position="none")
mdw_plot_1[[6]] <- mdw_plot_1[[6]] + theme(legend.position="none")

mdw_plot_2[[2]] <- mdw_plot_2[[2]] + theme(legend.position="none")
mdw_plot_2[[3]] <- mdw_plot_2[[3]] + theme(legend.position="none")
mdw_plot_2[[4]] <- mdw_plot_2[[4]] + theme(legend.position="none")
mdw_plot_2[[5]] <- mdw_plot_2[[5]] + theme(legend.position="none")
mdw_plot_2[[6]] <- mdw_plot_2[[6]] + theme(legend.position="none")

mdw_plot_3[[2]] <- mdw_plot_3[[2]] + theme(legend.position="none")
mdw_plot_3[[3]] <- mdw_plot_3[[3]] + theme(legend.position="none")
mdw_plot_3[[4]] <- mdw_plot_3[[4]] + theme(legend.position="none")
mdw_plot_3[[5]] <- mdw_plot_3[[5]] + theme(legend.position="none")
mdw_plot_3[[6]] <- mdw_plot_3[[6]] + theme(legend.position="none")

Br6522_ant$subject <- "Br6522 anterior"
Br6522_mid$subject <- "Br6522 middle"
Br8667_post$subject <- "Br8667 posterior"

combined_metadata <- rbind(Br6522_ant@meta.data, Br6522_mid@meta.data, Br8667_post@meta.data)
combined_metadata$sample <- combined_metadata$subject
combined_metadata$sample[combined_metadata$is_wrinkle] <- paste(combined_metadata$sample[combined_metadata$is_wrinkle], "(artifact)")
combined_metadata$sample <- factor(combined_metadata$sample, levels = c("Br8667 posterior", "Br8667 posterior (artifact)", 
                                                                           "Br6522 middle", "Br6522 middle (artifact)", 
                                                                           "Br6522 anterior", "Br6522 anterior (artifact)"))

stat.test1 <- combined_metadata %>%
  group_by(Layers) %>%
  wilcox_test(nCount_Spatial ~ sample, alternative = "less") %>%
  add_significance("p")
stat.test1 <- stat.test1 %>% 
  add_xy_position(x = "Layers", dodge = 0.8, fun = "median_mad")
stat.test1$group_combined <- paste(stat.test1$group1, stat.test1$group2, sep = "-")
stat.test1 <- stat.test1[stat.test1$group_combined %in% c("Br6522 anterior-Br6522 anterior (artifact)", "Br6522 middle-Br6522 middle (artifact)", "Br8667 posterior-Br8667 posterior (artifact)"), ]
p1 <- ggplot(combined_metadata) + geom_boxplot(aes(x = Layers, y = nCount_Spatial, fill = sample),outlier.shape = NA) + 
  theme_classic() + guides(fill=guide_legend(title="Spot category")) + theme(legend.direction="horizontal", axis.title.y = element_text(size = 10))+
  ylab("Library Size") + xlab("") + scale_fill_brewer(palette = "Paired") +  
  stat_pvalue_manual(stat.test1,  label = "{p.adj.signif}", tip.length = 0, hide.ns = TRUE) 

stat.test2 <- combined_metadata %>%
  group_by(Layers) %>%
  wilcox_test(nFeature_Spatial ~ sample, alternative = "less") %>%
  add_significance("p")
stat.test2 <- stat.test2 %>% 
  add_xy_position(x = "Layers", dodge = 0.8, fun = "median_mad")
stat.test2$group_combined <- paste(stat.test2$group1, stat.test2$group2, sep = "-")
stat.test2 <- stat.test2[stat.test2$group_combined %in% c("Br6522 anterior-Br6522 anterior (artifact)", "Br6522 middle-Br6522 middle (artifact)", "Br8667 posterior-Br8667 posterior (artifact)"), ]
p2 <- ggplot(combined_metadata) + geom_boxplot(aes(x = Layers, y = nFeature_Spatial, fill = sample),outlier.shape = NA) + 
  theme_classic() +  guides(fill=guide_legend(title="Spot category")) + theme(legend.direction="horizontal", axis.title.y = element_text(size = 10)) +
  ylab("Number of detected genes") + xlab("") + scale_fill_brewer(palette = "Paired") +  
  stat_pvalue_manual(stat.test2,  label = "{p.adj.signif}",tip.length = 0, hide.ns = TRUE)

stat.test3 <- combined_metadata %>%
  group_by(Layers) %>%
  wilcox_test(percent.mt ~ sample, alternative = "greater") %>%
  add_significance("p")
stat.test3 <- stat.test3 %>% 
  add_xy_position(x = "Layers", dodge = 0.8, fun = "median_mad")
stat.test3$group_combined <- paste(stat.test3$group1, stat.test3$group2, sep = "-")
stat.test3 <- stat.test3[stat.test3$group_combined %in% c("Br6522 anterior-Br6522 anterior (artifact)", "Br6522 middle-Br6522 middle (artifact)", "Br8667 posterior-Br8667 posterior (artifact)"), ]
p3 <- ggplot(combined_metadata) + geom_boxplot(aes(x = Layers, y = percent.mt, fill = sample),outlier.shape = NA) + 
  theme_classic() +  guides(fill=guide_legend(title="Spot category")) + theme(legend.direction="horizontal", axis.title.y = element_text(size = 10)) +
  ylab("Percent mito") + xlab("") + scale_fill_brewer(palette = "Paired") +  
  stat_pvalue_manual(stat.test3,  label = "{p.adj.signif}",tip.length = 0, hide.ns = TRUE)


stat.test4 <- combined_metadata %>%
  group_by(Layers) %>%
  wilcox_test(nCells ~ sample, alternative = "less") %>%
  add_significance("p")
stat.test4 <- stat.test4 %>% 
  add_xy_position(x = "Layers", dodge = 0.8, fun = "median_mad")
stat.test4$group_combined <- paste(stat.test4$group1, stat.test4$group2, sep = "-")
stat.test4 <- stat.test4[stat.test4$group_combined %in% c("Br6522 anterior-Br6522 anterior (artifact)", "Br6522 middle-Br6522 middle (artifact)", "Br8667 posterior-Br8667 posterior (artifact)"), ]
stat.test4$y.position <- stat.test4$y.position + 25
p4 <- ggplot(combined_metadata) + geom_boxplot(aes(x = Layers, y = nCells, fill = sample),outlier.shape = NA) + 
  theme_classic() + guides(fill=guide_legend(title="Spot category")) + theme(legend.direction="horizontal", axis.title.y = element_text(size = 10)) +
  ylab("#Cells (Vistoseg)") + xlab("") + scale_fill_brewer(palette = "Paired") +  
  stat_pvalue_manual(stat.test4,  label = "{p.adj.signif}",tip.length = 0, hide.ns = TRUE)

qc_legend <- get_legend(p1)

pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Library_size.pdf", width = 9, height = 4)
p1 + theme(legend.position = "none")
dev.off()
pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Number_of_detected_genes.pdf", width = 9, height = 4)
p2 + theme(legend.position = "none")
dev.off()
pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Percent_mito.pdf", width = 9, height = 4)
p3 + theme(legend.position = "none")
dev.off()
pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Number_of_cells.pdf", width = 9, height = 4)
p4 + theme(legend.position = "none")
dev.off()
pdf(file = "/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/legend_QC.pdf", width = 9, height = 0.5)
plot_grid(qc_legend)
dev.off()

ggplotRegression <- function (fit, x, y) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) + theme_classic() + 
                        xlab(x) + ylab(y) + theme(plot.title = element_text(size = 10), axis.title = element_text(size = 10))
}

fit_wrinkle <- lm(nCount_Spatial ~ nCells, Br6522_ant@meta.data[Br6522_ant$is_wrinkle, ])
fit_normal <- lm(nCount_Spatial ~ nCells, Br6522_ant@meta.data[!Br6522_ant$is_wrinkle, ])
plot_grid(ggplotRegression(fit_normal, "Number of cells", "Library size"), ggplotRegression(fit_wrinkle, "Number of cells", "Library size"), nrow = 2)

fit_wrinkle <- lm(nCount_Spatial ~ nCells, Br6522_mid@meta.data[Br6522_mid$is_wrinkle, ])
fit_normal <- lm(nCount_Spatial ~ nCells, Br6522_mid@meta.data[!Br6522_mid$is_wrinkle, ])
plot_grid(ggplotRegression(fit_normal, "Number of cells", "Library size"), ggplotRegression(fit_wrinkle, "Number of cells", "Library size"), nrow = 2)

fit_wrinkle <- lm(nCount_Spatial ~ nCells, Br8667_post@meta.data[Br8667_post$is_wrinkle, ])
fit_normal <- lm(nCount_Spatial ~ nCells, Br8667_post@meta.data[!Br8667_post$is_wrinkle, ])
plot_grid(ggplotRegression(fit_normal, "Number of cells", "Library size"), ggplotRegression(fit_wrinkle, "Number of cells", "Library size"), nrow = 2)

set.seed(100)
make_plots <- function(data){
  gene_attr <- data.frame(mean = rowMeans(data), detection_rate = rowMeans(data > 0), var = apply(data, 1, var))
  gene_attr$log_mean <- log10(gene_attr$mean)
  gene_attr$log_var <- log10(gene_attr$var)
  rownames(gene_attr) <- rownames(data)
  cell_attr <- data.frame(n_umi = colSums(data), n_gene = colSums(data > 0))
  rownames(cell_attr) <- colnames(data)
  p1 <- ggplot(gene_attr, aes(log_mean, log_var)) + geom_point(alpha = 0.3, shape = 16) + 
    geom_density_2d(size = 0.3) + geom_abline(intercept = 0, slope = 1, color = "red")+ theme_classic()
  x = seq(from = -3, to = 2, length.out = 1000)
  poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
  p2 <- ggplot(gene_attr, aes(log_mean, detection_rate)) + geom_point(alpha = 0.3, shape = 16) + 
    geom_line(data = poisson_model, color = "red") + theme_gray(base_size = 8) + theme_classic()
  p3 <- ggplot(cell_attr, aes(n_umi, n_gene)) + geom_point(alpha = 0.3, shape = 16) + geom_density_2d(size = 0.3) + theme_classic()
  p1 + p2 + p3
}

geneData <- readRDS("/data/abattle4/prashanthi/dewrinkler/data/gene_df.rds")
Br6522_ant@assays$Spatial@meta.features$nCells <- rowSums(Br6522_ant@assays$Spatial@counts > 0)
Br6522_mid@assays$Spatial@meta.features$nCells <- rowSums(Br6522_mid@assays$Spatial@counts > 0)
Br8667_post@assays$Spatial@meta.features$nCells <- rowSums(Br8667_post@assays$Spatial@counts > 0)
Br6522_ant@assays$Spatial@meta.features$geneType <- geneData$gene_type
Br6522_mid@assays$Spatial@meta.features$geneType <- geneData$gene_type
Br8667_post@assays$Spatial@meta.features$geneType <- geneData$gene_type

Br6522_ant@assays$Spatial@meta.features$pCells <- Br6522_ant@assays$Spatial@meta.features$nCells/dim(Br6522_ant)[2]
Br6522_mid@assays$Spatial@meta.features$pCells <- Br6522_mid@assays$Spatial@meta.features$nCells/dim(Br6522_mid)[2]
Br8667_post@assays$Spatial@meta.features$pCells <- Br8667_post@assays$Spatial@meta.features$nCells/dim(Br8667_post)[2]

# Apply QC filters
# Gene level QC
Br6522_ant <- Br6522_ant[Br6522_ant@assays$Spatial@meta.features$geneType == "protein_coding", ]
Br6522_mid <- Br6522_mid[Br6522_mid@assays$Spatial@meta.features$geneType == "protein_coding", ]
Br8667_post <- Br8667_post[Br8667_post@assays$Spatial@meta.features$geneType == "protein_coding", ]

Br6522_ant <- Br6522_ant[!grepl("^MT-", rownames(Br6522_ant)), ]
Br6522_mid <- Br6522_mid[!grepl("^MT-", rownames(Br6522_mid)), ]
Br8667_post <- Br8667_post[!grepl("^MT-", rownames(Br8667_post)), ]

Br6522_ant <- Br6522_ant[Br6522_ant@assays$Spatial@meta.features$pCells >= 0.05, ]
Br6522_mid <- Br6522_mid[Br6522_mid@assays$Spatial@meta.features$pCells >= 0.05, ]
Br8667_post <- Br8667_post[Br8667_post@assays$Spatial@meta.features$pCells >= 0.05, ]

Br6522_ant_normal <- subset(Br6522_ant, Wrinkles == "None")
Br6522_mid_normal <- subset(Br6522_mid, Wrinkles == "None")
Br8667_post_normal <- subset(Br8667_post, Wrinkles == "None")

Br6522_ant_normal@assays$Spatial@meta.features$nCells <- rowSums(Br6522_ant_normal@assays$Spatial@counts > 0)
Br6522_mid_normal@assays$Spatial@meta.features$nCells <- rowSums(Br6522_mid_normal@assays$Spatial@counts > 0)
Br8667_post_normal@assays$Spatial@meta.features$nCells <- rowSums(Br8667_post_normal@assays$Spatial@counts > 0)

Br6522_ant_normal@assays$Spatial@meta.features$pCells <- Br6522_ant_normal@assays$Spatial@meta.features$nCells/dim(Br6522_ant_normal)[2]
Br6522_mid_normal@assays$Spatial@meta.features$pCells <- Br6522_mid_normal@assays$Spatial@meta.features$nCells/dim(Br6522_mid_normal)[2]
Br8667_post_normal@assays$Spatial@meta.features$pCells <- Br8667_post_normal@assays$Spatial@meta.features$nCells/dim(Br8667_post_normal)[2]

Br6522_ant_normal <- Br6522_ant_normal[Br6522_ant_normal@assays$Spatial@meta.features$pCells >= 0.05, ]
Br6522_mid_normal <- Br6522_mid_normal[Br6522_mid_normal@assays$Spatial@meta.features$pCells >= 0.05, ]
Br8667_post_normal <- Br8667_post_normal[Br8667_post_normal@assays$Spatial@meta.features$pCells >= 0.05, ]

Br6522_ant_all_data <- Br6522_ant@assays$Spatial@counts
Br6522_mid_all_data <- Br6522_mid@assays$Spatial@counts
Br8667_post_all_data <- Br8667_post@assays$Spatial@counts

Br6522_ant_normal_data <- Br6522_ant_normal@assays$Spatial@counts
Br6522_mid_normal_data <- Br6522_mid_normal@assays$Spatial@counts
Br8667_post_normal_data <- Br8667_post_normal@assays$Spatial@counts

make_plots(Br6522_ant_all_data)
make_plots(Br6522_mid_all_data)
make_plots(Br8667_post_all_data)

make_plots(Br6522_ant_normal_data)
make_plots(Br6522_mid_normal_data)
make_plots(Br8667_post_normal_data)

Br6522_ant_normal_vst_out <- sctransform::vst(Br6522_ant_normal_data,latent_var = c("log_umi"), return_gene_attr = TRUE, 
                            return_cell_attr = TRUE, n_genes = NULL,  verbosity = 1)
Br6522_mid_normal_vst_out <- sctransform::vst(Br6522_mid_normal_data, latent_var = c("log_umi"), return_gene_attr = TRUE, 
                                       return_cell_attr = TRUE,n_genes = NULL,  verbosity = 1)
Br8667_post_normal_vst_out <- sctransform::vst(Br8667_post_normal_data, latent_var = c("log_umi"), return_gene_attr = TRUE, 
                                       return_cell_attr = TRUE, n_genes = NULL, verbosity = 1)


Br6522_ant_additional_params <- data.frame(Br6522_ant$is_wrinkle)
colnames(Br6522_ant_additional_params) <- "is_wrinkle"
Br6522_ant_all_vst_out <- sctransform::vst(Br6522_ant_all_data, cell_attr = Br6522_ant_additional_params,latent_var = c("log_umi", "is_wrinkle"), return_gene_attr = TRUE, 
                                              return_cell_attr = TRUE, n_genes = NULL,  verbosity = 1)
Br6522_ant_all_vst_out_libSize_only <- sctransform::vst(Br6522_ant_all_data, cell_attr = Br6522_ant_additional_params,latent_var = c("log_umi"), return_gene_attr = TRUE, 
                                           return_cell_attr = TRUE, n_genes = NULL,  verbosity = 1)
Br6522_ant_all_vst_out_wrinkle_only <- sctransform::vst(Br6522_ant_all_data, cell_attr = Br6522_ant_additional_params,latent_var = c("is_wrinkle"), return_gene_attr = TRUE, 
                                                        return_cell_attr = TRUE, n_genes = NULL,  verbosity = 1)


Br6522_mid_additional_params <- data.frame(Br6522_mid$is_wrinkle)
colnames(Br6522_mid_additional_params) <- "is_wrinkle"
Br6522_mid_all_vst_out <- sctransform::vst(Br6522_mid_all_data, cell_attr = Br6522_mid_additional_params, latent_var = c("log_umi", "is_wrinkle"), return_gene_attr = TRUE, 
                                              return_cell_attr = TRUE,n_genes = NULL,  verbosity = 1)
Br6522_mid_all_vst_out_libSize_only <- sctransform::vst(Br6522_mid_all_data, cell_attr = Br6522_mid_additional_params, latent_var = c("log_umi"), return_gene_attr = TRUE, 
                                           return_cell_attr = TRUE,n_genes = NULL,  verbosity = 1)
Br6522_mid_all_vst_out_wrinkle_only <- sctransform::vst(Br6522_mid_all_data, cell_attr = Br6522_mid_additional_params, latent_var = c("is_wrinkle"), return_gene_attr = TRUE, 
                                                        return_cell_attr = TRUE,n_genes = NULL,  verbosity = 1)

Br8667_post_additional_params <- data.frame(Br8667_post$is_wrinkle)
colnames(Br8667_post_additional_params) <- "is_wrinkle"
Br8667_post_all_vst_out <- sctransform::vst(Br8667_post_all_data, cell_attr = Br8667_post_additional_params,latent_var = c("log_umi", "is_wrinkle"), return_gene_attr = TRUE, 
                                               return_cell_attr = TRUE, n_genes = NULL, verbosity = 1)
Br8667_post_all_vst_out_libSize_only <- sctransform::vst(Br8667_post_all_data, cell_attr = Br8667_post_additional_params, latent_var = c("log_umi"), return_gene_attr = TRUE, 
                                                        return_cell_attr = TRUE,n_genes = NULL,  verbosity = 1)
Br8667_post_all_vst_out_wrinkle_only <- sctransform::vst(Br8667_post_all_data, cell_attr = Br8667_post_additional_params, latent_var = c("is_wrinkle"), return_gene_attr = TRUE, 
                                                        return_cell_attr = TRUE,n_genes = NULL,  verbosity = 1)


sctransform::plot_model_pars(Br6522_ant_normal_vst_out, show_theta = TRUE)
sctransform::plot_model_pars(Br6522_mid_normal_vst_out, show_theta = TRUE)
sctransform::plot_model_pars(Br8667_post_normal_vst_out, show_theta = TRUE)

make_compare_plots <- function(vst_out, vst_out_wrinkle_only, sample_label){
  wrinkle_only_params_df <- data.frame(rbind(cbind(vst_out_wrinkle_only$gene_attr$gmean, vst_out_wrinkle_only$model_pars, "Single gene estimate"), 
                                             cbind(vst_out_wrinkle_only$gene_attr$gmean, vst_out_wrinkle_only$model_pars_fit, "Regularized")))
  colnames(wrinkle_only_params_df)[1] <- "gmean"
  colnames(wrinkle_only_params_df)[5] <- "estimate_type"
  wrinkle_only_params_df$gmean <- as.numeric(wrinkle_only_params_df$gmean)
  wrinkle_only_params_df$theta <- as.numeric(wrinkle_only_params_df$theta)
  wrinkle_only_params_df$X.Intercept. <- as.numeric(wrinkle_only_params_df$X.Intercept.)
  wrinkle_only_params_df$is_wrinkleTRUE <- as.numeric(wrinkle_only_params_df$is_wrinkleTRUE)
  wrinkle_only_params_df$estimate_type <- factor(wrinkle_only_params_df$estimate_type, levels = c("Regularized", "Single gene estimate"))
  wrinkle_only_params_df$gmean <- log10(wrinkle_only_params_df$gmean)
  p0 <- ggplot(wrinkle_only_params_df, aes(x = gmean, y = is_wrinkleTRUE, colour = estimate_type)) + 
    geom_point(alpha = 0.5, size = 0.7) + theme_classic() +
    xlab("Geometric mean of gene [log10]") + ylab(expression(beta[2])) + geom_hline(yintercept = 0, lty = 2, color = "black") +
    guides(colour=guide_legend(title= "Estimate type")) + theme(axis.title.x = element_text(size = 10))
  
  all_params_df <- data.frame(rbind(cbind(vst_out$gene_attr$gmean, vst_out$model_pars, "Single gene estimate"), cbind(vst_out$gene_attr$gmean, vst_out$model_pars_fit, "Regularized")))
  colnames(all_params_df)[1] <- "gmean"
  colnames(all_params_df)[6] <- "estimate_type"
  all_params_df$gmean <- as.numeric(all_params_df$gmean)
  all_params_df$theta <- as.numeric(all_params_df$theta)
  all_params_df$X.Intercept. <- as.numeric(all_params_df$X.Intercept.)
  all_params_df$log_umi <- as.numeric(all_params_df$log_umi)
  all_params_df$is_wrinkleTRUE <- as.numeric(all_params_df$is_wrinkleTRUE)
  all_params_df$estimate_type <- factor(all_params_df$estimate_type, levels = c("Regularized", "Single gene estimate"))
  all_params_df$gmean <- log10(all_params_df$gmean)
  p1 <- ggplot(all_params_df, aes(x = gmean, y = log_umi, colour = estimate_type)) + 
    geom_point(alpha = 0.5, size = 0.7) + theme_classic() +
    xlab("Geometric mean of gene [log10]") + ylab(expression(beta[1])) + guides(colour=guide_legend(title= "Estimate type")) +
    theme(axis.title.x = element_text(size = 10)) + geom_hline(yintercept = 0, lty = 2, color = "black")
  p2 <- ggplot(all_params_df, aes(x = gmean, y = is_wrinkleTRUE, colour = estimate_type)) + 
    geom_point(alpha = 0.5, size = 0.7) + theme_classic() +
    xlab("Geometric mean of gene [log10]") + ylab(expression(beta[2])) + guides(colour=guide_legend(title= "Estimate type")) +
    theme(axis.title.x = element_text(size = 10)) + geom_hline(yintercept = 0, lty = 2, color = "black")
  estimate_type_legend <- get_legend(p0)
  p0 <- p0 + theme(legend.position = "none")
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  
  title_model_1 <- ggdraw() + draw_label(expression("Model 1: log(E["~x[i]~"])" == ~ beta[0] ~ + ~ beta[2] ~ w), size = 12)
  title_model_2 <- ggdraw() + draw_label(expression("Model 2: log(E["~x[i]~"])" == ~ beta[0] ~ + ~ beta[1] ~ m ~ + ~ beta[2] ~ w), size = 12)
  sample_name <- ggdraw() + draw_label(sample_label, size = 12)
  
  plot_grid(sample_name, 
            plot_grid(title_model_1, title_model_2, p0, plot_grid(p2, p1), rel_heights = c(1, 5), rel_widths = c(1, 2)), 
            estimate_type_legend, rel_widths = c(1, 6, 1), nrow = 1)
}

p1 <- make_compare_plots(Br6522_ant_all_vst_out, Br6522_ant_all_vst_out_wrinkle_only, "Br6522 Anterior")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br6522_ant_NB_regression.pdf", width = 12, height = 3)
p1
dev.off()
p2 <- make_compare_plots(Br6522_mid_all_vst_out, Br6522_mid_all_vst_out_wrinkle_only, "Br6522 Middle")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br6522_mid_NB_regression.pdf", width = 12, height = 3)
p2
dev.off()
p3 <- make_compare_plots(Br8667_post_all_vst_out, Br8667_post_all_vst_out_wrinkle_only, "Br8667 Posterior")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br8667_post_NB_regression.pdf", width = 12, height = 3)
p3
dev.off()

plot_params <- function(normal_vst_out, all_vst_out, title){
  common.genes <-  intersect(rownames(normal_vst_out$model_pars), rownames(all_vst_out$model_pars))
  normal_params <- normal_vst_out$model_pars
  normal_params_fit <- normal_vst_out$model_pars_fit
  all_params <- all_vst_out$model_pars
  all_params_fit <- all_vst_out$model_pars_fit
  normal_params <- normal_params[match(common.genes, rownames(normal_params)), ]
  normal_params_fit <- normal_params_fit[match(common.genes, rownames(normal_params_fit)), ]
  all_params <- all_params[match(common.genes, rownames(all_params)), ]
  all_params_fit <- all_params_fit[match(common.genes, rownames(all_params_fit)), ]
  plot_df <- data.frame("Normal_intercept" = normal_params[ ,2], 
                        "Normal_intercept_regularized" = normal_params_fit[ ,2], 
                        "All_intercept" = all_params[ ,2], 
                        "All_intercept_regularized" = all_params_fit[ ,2])

  ggplot(plot_df)  +  geom_point(aes(x = Normal_intercept, y = All_intercept, colour = "Single gene estimate"), alpha = 0.5) + 
    geom_point(aes(x = Normal_intercept_regularized, y = All_intercept_regularized, colour = "Regularized"), alpha = 0.5) + geom_abline(slope = 1, intercept = 0, lty = 2) + 
    theme_classic() + xlab("Intercept (excl. artifacts)") + ylab("Intercept") + ggtitle(title) + 
    scale_colour_manual(name = "Estimate type", values = c("Regularized" = "#F8766D", "Single gene estimate" = "#00BFC4")) +
    theme(plot.title = element_text(face = "bold", size = 14), axis.title = element_text(size = 12), legend.title = element_text(size = 12), 
          legend.text = element_text(size = 10))
}

p <-  plot_params(Br6522_ant_normal_vst_out, Br6522_ant_all_vst_out_libSize_only, "Br6522 Anterior")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br6522_ant_intercept.pdf", width = 7, height = 5)
p
dev.off()

p <-  plot_params(Br6522_mid_normal_vst_out, Br6522_mid_all_vst_out_libSize_only, "Br6522 Middle")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br6522_mid_intercept.pdf", width = 7, height = 5)
p
dev.off()

p <-  plot_params(Br8667_post_normal_vst_out, Br8667_post_all_vst_out_libSize_only, "Br8667 Posterior")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S4/Br8667_post_intercept.pdf", width = 7, height = 5)
p
dev.off()

