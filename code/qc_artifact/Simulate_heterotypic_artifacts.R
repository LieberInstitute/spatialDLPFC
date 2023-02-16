rm(list = ls())
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(rstatix)
library(cowplot)
library(DoubletFinder)
set.seed(1001)

manually_annotate <- function(sobj, layers_df, wrinkles_df, sc_df) {
    spots <- rownames(sobj@meta.data)
    layers <- c()
    wrinkle <- c()
    nCells <- c()
    for (i in c(1:length(spots))) {
        if (spots[i] %in% layers_df$spot_name) {
            layers[i] <- layers_df$ManualAnnotation[layers_df$spot_name == spots[i]]
        } else {
            layers[i] <- "Unknown"
        }
        if (spots[i] %in% wrinkles_df$spot_name) {
            wrinkle[i] <- wrinkles_df$ManualAnnotation[wrinkles_df$spot_name == spots[i]]
        } else {
            wrinkle[i] <- "None"
        }
        nCells[i] <- sc_df$Nmask_dark_blue[sc_df$barcode == spots[i]]
    }
    sobj@meta.data["Layers"] <- layers
    sobj@meta.data["Wrinkles"] <- wrinkle
    sobj@meta.data["nCells"] <- nCells
    sobj
}

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

Br6522_ant <- manually_annotate(Br6522_ant, Br6522_ant_layers, Br6522_ant_wrinkles, Br6522_ant_spot_counts)
Br6522_mid <- manually_annotate(Br6522_mid, Br6522_mid_layers, Br6522_mid_wrinkles, Br6522_mid_spot_counts)
Br8667_post <- manually_annotate(Br8667_post, Br8667_post_layers, Br8667_post_wrinkles, Br8667_post_spot_counts)

# Exclude points that do not have a definitive layer assignment
Br6522_ant <- Br6522_ant[, !Br6522_ant$Layers == "Unknown"]
Br6522_mid <- Br6522_mid[, !Br6522_mid$Layers == "Unknown"]
Br8667_post <- Br8667_post[, !Br8667_post$Layers == "Unknown"]

# Format the metadata
Br6522_ant$Layers <- factor(Br6522_ant$Layers, levels = c("Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "WM"))
Br6522_mid$Layers <- factor(Br6522_mid$Layers, levels = c("Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "WM"))
Br8667_post$Layers <- factor(Br8667_post$Layers, levels = c("Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6"))

Br6522_ant$Wrinkles <- factor(Br6522_ant$Wrinkles, levels = c(
    "None", "Fold_1", "Shear_1", "Shear_2", "Shear_3", "Wrinkle_1", "Wrinkle_2",
    "Wrinkle_3", "Wrinkle_4", "Wrinkle_5", "Wrinkle_6", "Wrinkle_7", "Wrinkle_8",
    "Wrinkle_9"
))
Br6522_mid$Wrinkles <- factor(Br6522_mid$Wrinkles, levels = c(
    "None", "Fold_1", "Shear_1", "Shear_2", "Wrinkle_1", "Wrinkle_2",
    "Wrinkle_3", "Wrinkle_4", "Wrinkle_5", "Wrinkle_6", "Wrinkle_7", "Wrinkle_8",
    "Wrinkle_9", "Wrinkle_10", "Wrinkle_11", "Wrinkle_12", "Wrinkle_13"
))
Br8667_post$Wrinkles <- factor(Br8667_post$Wrinkles, levels = c(
    "None", "Shear_1", "Wrinkle_1", "Wrinkle_2",
    "Wrinkle_3", "Wrinkle_4", "Wrinkle_5", "Wrinkle_6", "Wrinkle_7", "Wrinkle_8",
    "Wrinkle_9", "Wrinkle_11", "Wrinkle_12", "Wrinkle_13"
))
Br6522_ant$is_wrinkle <- !Br6522_ant$Wrinkles == "None"
Br6522_mid$is_wrinkle <- !Br6522_mid$Wrinkles == "None"
Br8667_post$is_wrinkle <- !Br8667_post$Wrinkles == "None"

simulate_doublets <- function(sobj) {
    sobj <- subset(sobj, Wrinkles == "None")
    nSim <- round(0.05 * dim(sobj)[2])
    layer_combinations <- combn(levels(sobj$Layers), 2)
    simulated_doublets <- list()
    for (i in c(1:dim(layer_combinations)[2])) {
        print(i)
        L1 <- layer_combinations[1, i]
        L2 <- layer_combinations[2, i]
        L1_sobj <- subset(sobj, Layers == L1)
        L2_sobj <- subset(sobj, Layers == L2)
        prop <- runif(nSim, min = 0.3, max = 0.7)
        # prop <- rep(0.5, nSim)
        sim_Dat <- matrix(NA, nrow = dim(L1_sobj)[1], ncol = nSim)
        dat1 <- L1_sobj@assays$Spatial@counts
        dat2 <- L2_sobj@assays$Spatial@counts
        for (j in c(1:nSim)) {
            sim_Dat[, j] <- (prop[j] * dat1[, sample(c(1:dim(dat1)[2]), 1)]) + ((1 - prop[j]) * dat2[, sample(c(1:dim(dat2)[2]), 1)])
        }
        sim_Dat <- data.frame(t(sim_Dat))
        colnames(sim_Dat) <- rownames(sobj)
        rownames(sim_Dat) <- paste(gsub(" ", "", L1), gsub(" ", "", L2), c(1:dim(sim_Dat)[1]), sep = "_")
        sim_Dat$L1 <- L1
        sim_Dat$L2 <- L2
        sim_Dat$Layer <- paste(sim_Dat$L1, sim_Dat$L2, sep = "/")
        sim_Dat$prop_L1 <- prop
        sim_Dat$Type <- "Simulated"
        simulated_doublets[[i]] <- sim_Dat
        rm(sim_Dat)
        rm(prop)
    }
    simulated_doublets <- do.call(rbind, simulated_doublets)
    simulated_doublets
}

Br6522_ant_simulated <- simulate_doublets(Br6522_ant)
Br6522_mid_simulated <- simulate_doublets(Br6522_mid)
Br8667_post_simulated <- simulate_doublets(Br8667_post)

Br6522_ant$Type <- ifelse(Br6522_ant$Wrinkles == "None", "Non-artifact", "Artifact")
Br6522_mid$Type <- ifelse(Br6522_mid$Wrinkles == "None", "Non-artifact", "Artifact")
Br8667_post$Type <- ifelse(Br8667_post$Wrinkles == "None", "Non-artifact", "Artifact")


Br6522_ant_expr <- cbind(Br6522_ant@assays$Spatial@counts, t(Br6522_ant_simulated[, 1:36601]))
Br6522_mid_expr <- cbind(Br6522_mid@assays$Spatial@counts, t(Br6522_mid_simulated[, 1:36601]))
Br8667_post_expr <- cbind(Br8667_post@assays$Spatial@counts, t(Br8667_post_simulated[, 1:36601]))

Br6522_ant_all <- CreateSeuratObject(Br6522_ant_expr, project = "Br6522 Anterior")
Br6522_mid_all <- CreateSeuratObject(Br6522_mid_expr, project = "Br6522 Middle")
Br8667_post_all <- CreateSeuratObject(Br8667_post_expr, project = "Br8667 Posterior")

Br6522_ant_all$Type <- c(as.character(Br6522_ant$Type), rep("Simulated doublet", dim(Br6522_ant_simulated)[1]))
Br6522_ant_all$prop_L1 <- c(rep(1, dim(Br6522_ant@assays$Spatial@counts)[2]), as.numeric(Br6522_ant_simulated$prop_L1))
Br6522_ant_all$Wrinkles <- c(as.character(Br6522_ant$Wrinkles), rep("Doublet", dim(Br6522_ant_simulated)[1]))
Br6522_ant_all$Layers <- c(as.character(Br6522_ant$Layers), as.character(Br6522_ant_simulated$Layer))
Br6522_ant_all$Layers <- as.factor(Br6522_ant_all$Layers)
Br6522_ant_all$Type <- as.factor(Br6522_ant_all$Type)
Br6522_ant_all$Wrinkles <- as.factor(Br6522_ant_all$Wrinkles)

Br6522_mid_all$Type <- c(as.character(Br6522_mid$Type), rep("Simulated doublet", dim(Br6522_mid_simulated)[1]))
Br6522_mid_all$prop_L1 <- c(rep(1, dim(Br6522_mid@assays$Spatial@counts)[2]), as.numeric(Br6522_mid_simulated$prop_L1))
Br6522_mid_all$Wrinkles <- c(as.character(Br6522_mid$Wrinkles), rep("Doublet", dim(Br6522_mid_simulated)[1]))
Br6522_mid_all$Layers <- c(as.character(Br6522_mid$Layers), as.character(Br6522_mid_simulated$Layer))
Br6522_mid_all$Layers <- as.factor(Br6522_mid_all$Layers)
Br6522_mid_all$Type <- as.factor(Br6522_mid_all$Type)
Br6522_mid_all$Wrinkles <- as.factor(Br6522_mid_all$Wrinkles)

Br8667_post_all$Type <- c(as.character(Br8667_post$Type), rep("Simulated doublet", dim(Br8667_post_simulated)[1]))
Br8667_post_all$prop_L1 <- c(rep(1, dim(Br8667_post@assays$Spatial@counts)[2]), as.numeric(Br8667_post_simulated$prop_L1))
Br8667_post_all$Wrinkles <- c(as.character(Br8667_post$Wrinkles), rep("Doublet", dim(Br8667_post_simulated)[1]))
Br8667_post_all$Layers <- c(as.character(Br8667_post$Layers), as.character(Br8667_post_simulated$Layer))
Br8667_post_all$Layers <- as.factor(Br8667_post_all$Layers)
Br8667_post_all$Type <- as.factor(Br8667_post_all$Type)
Br8667_post_all$Wrinkles <- as.factor(Br8667_post_all$Wrinkles)

Br6522_ant_all <- PercentageFeatureSet(Br6522_ant_all, "^MT-", col.name = "percent_mito")
Br6522_ant_all <- PercentageFeatureSet(Br6522_ant_all, "^RP[SL]", col.name = "percent_ribo")
selected_cells <- WhichCells(Br6522_ant_all, expression = nFeature_RNA > 200)
selected_f <- rownames(selected_cells)[Matrix::rowSums(Br6522_ant_all) > 3]
Br6522_ant_all <- subset(Br6522_ant_all, features = selected_f, cells = selected_cells)
Br6522_ant_all <- Br6522_ant_all[!grepl("MALAT1", rownames(Br6522_ant_all)), ]
Br6522_ant_all <- Br6522_ant_all[!grepl("^MT-", rownames(Br6522_ant_all)), ]
Br6522_ant_all <- Br6522_ant_all[!grepl("^RP[SL]", rownames(Br6522_ant_all)), ]
Br6522_ant_all <- NormalizeData(Br6522_ant_all)

Br6522_mid_all <- PercentageFeatureSet(Br6522_mid_all, "^MT-", col.name = "percent_mito")
Br6522_mid_all <- PercentageFeatureSet(Br6522_mid_all, "^RP[SL]", col.name = "percent_ribo")
selected_cells <- WhichCells(Br6522_mid_all, expression = nFeature_RNA > 200)
selected_f <- rownames(selected_cells)[Matrix::rowSums(Br6522_mid_all) > 3]
Br6522_mid_all <- subset(Br6522_mid_all, features = selected_f, cells = selected_cells)
Br6522_mid_all <- Br6522_mid_all[!grepl("MALAT1", rownames(Br6522_mid_all)), ]
Br6522_mid_all <- Br6522_mid_all[!grepl("^MT-", rownames(Br6522_mid_all)), ]
Br6522_mid_all <- Br6522_mid_all[!grepl("^RP[SL]", rownames(Br6522_mid_all)), ]
Br6522_mid_all <- NormalizeData(Br6522_mid_all)

Br8667_post_all <- PercentageFeatureSet(Br8667_post_all, "^MT-", col.name = "percent_mito")
Br8667_post_all <- PercentageFeatureSet(Br8667_post_all, "^RP[SL]", col.name = "percent_ribo")
selected_cells <- WhichCells(Br8667_post_all, expression = nFeature_RNA > 200)
selected_f <- rownames(selected_cells)[Matrix::rowSums(Br8667_post_all) > 3]
Br8667_post_all <- subset(Br8667_post_all, features = selected_f, cells = selected_cells)
Br8667_post_all <- Br8667_post_all[!grepl("MALAT1", rownames(Br8667_post_all)), ]
Br8667_post_all <- Br8667_post_all[!grepl("^MT-", rownames(Br8667_post_all)), ]
Br8667_post_all <- Br8667_post_all[!grepl("^RP[SL]", rownames(Br8667_post_all)), ]
Br8667_post_all <- NormalizeData(Br8667_post_all)

Br6522_ant_all <- FindVariableFeatures(Br6522_ant_all, verbose = F)
Br6522_ant_all <- ScaleData(Br6522_ant_all,
    vars.to.regress = c("nFeature_RNA", "percent_mito"),
    verbose = F
)
Br6522_ant_all <- RunPCA(Br6522_ant_all, verbose = FALSE)
stdv <- Br6522_ant[["pca"]]@stdev
sum.stdv <- sum(Br6522_ant[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(
    which((percent.stdv[1:length(percent.stdv) - 1] -
        percent.stdv[2:length(percent.stdv)]) > 0.1),
    decreasing = T
)[1] + 1
min.pc <- min(co1, co2)
Br6522_ant_all <- FindNeighbors(Br6522_ant_all, reduction = "pca", dims = 1:8, return.neighbor = TRUE)

Br6522_mid_all <- FindVariableFeatures(Br6522_mid_all, verbose = F)
Br6522_mid_all <- ScaleData(Br6522_mid_all,
    vars.to.regress = c("nFeature_RNA", "percent_mito"),
    verbose = F
)
Br6522_mid_all <- RunPCA(Br6522_mid_all, verbose = FALSE)
stdv <- Br6522_mid[["pca"]]@stdev
sum.stdv <- sum(Br6522_mid[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(
    which((percent.stdv[1:length(percent.stdv) - 1] -
        percent.stdv[2:length(percent.stdv)]) > 0.1),
    decreasing = T
)[1] + 1
min.pc <- min(co1, co2)
Br6522_mid_all <- FindNeighbors(Br6522_mid_all, reduction = "pca", dims = 1:8, return.neighbor = TRUE)

Br8667_post_all <- FindVariableFeatures(Br8667_post_all, verbose = F)
Br8667_post_all <- ScaleData(Br8667_post_all,
    vars.to.regress = c("nFeature_RNA", "percent_mito"),
    verbose = F
)
Br8667_post_all <- RunPCA(Br8667_post_all, verbose = FALSE)
stdv <- Br8667_post[["pca"]]@stdev
sum.stdv <- sum(Br8667_post[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(
    which((percent.stdv[1:length(percent.stdv) - 1] -
        percent.stdv[2:length(percent.stdv)]) > 0.1),
    decreasing = T
)[1] + 1
min.pc <- min(co1, co2)
Br8667_post_all <- FindNeighbors(Br8667_post_all, reduction = "pca", dims = 1:8, return.neighbor = TRUE)

add_info_neighbors <- function(sobj) {
    frac_layers <- matrix(NA, nrow = dim(sobj)[2], ncol = length(levels(sobj$Layers)))
    for (i in c(1:dim(sobj)[2])) {
        layers <- sobj$Layers[sobj@neighbors$RNA.nn@nn.idx[i, 1:10]]
        frac_layers[i, ] <- table(layers)
    }
    colnames(frac_layers) <- levels(sobj$Layers)
    sobj@meta.data <- cbind(sobj@meta.data, frac_layers)
    frac_type <- matrix(NA, nrow = dim(sobj)[2], ncol = length(levels(sobj$Type)))
    for (i in c(1:dim(sobj)[2])) {
        type <- sobj$Type[sobj@neighbors$RNA.nn@nn.idx[i, 1:10]]
        frac_type[i, ] <- table(type)
    }
    colnames(frac_type) <- levels(sobj$Type)
    sobj@meta.data <- cbind(sobj@meta.data, frac_type)
    sobj
}

Br6522_ant_all <- add_info_neighbors(Br6522_ant_all)
Br6522_mid_all <- add_info_neighbors(Br6522_mid_all)
Br8667_post_all <- add_info_neighbors(Br8667_post_all)

df <- Br6522_ant_all@meta.data
df <- df[!df$Type == "Simulated doublet", ]
df$Layer_Wrinkle <- paste(df$Layers, as.numeric(!df$Wrinkles == "None"), sep = "-")
df2 <- df[, c(6, 38, 39, 40)] %>%
    group_by(Wrinkles) %>%
    summarise(across(everything(), mean),
        .groups = "drop"
    ) %>%
    as.data.frame()
df3 <- df[, c(6, 38, 39, 40)] %>%
    group_by(Wrinkles) %>%
    summarise(across(everything(), sd),
        .groups = "drop"
    ) %>%
    as.data.frame()
df2 <- reshape2::melt(df2)
df3 <- reshape2::melt(df3)
df2$sd <- df3$value
df2$Wrinkles <- gsub("_", " ", df2$Wrinkles)
df2$Wrinkles <- factor(df2$Wrinkles, levels = c(
    "None", "Fold 1", "Shear 1", "Shear 2",
    "Shear 3", "Wrinkle 1", "Wrinkle 2",
    "Wrinkle 3", "Wrinkle 4", "Wrinkle 5",
    "Wrinkle 6", "Wrinkle 7", "Wrinkle 8",
    "Wrinkle 9"
))
df2$variable <- factor(df2$variable, levels = c("Simulated doublet", "Artifact", "Non-artifact"))
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/Br6522_ant_SNN_artifacts.pdf", width = 13, height = 3)
ggplot(df2, aes(x = Wrinkles, y = value, fill = variable)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin = value - 0.5 * sd, ymax = value + 0.5 * sd), width = .2, position = position_dodge(.9)) +
    theme_classic() +
    xlab("") +
    ylab("Proportion of NN") +
    ggtitle("Br6522 Anterior") +
    theme(legend.title = element_blank(), plot.title = element_text(face = "bold"))
dev.off()

df <- Br6522_mid_all@meta.data
df <- df[!df$Type == "Simulated doublet", ]
df$Layer_Wrinkle <- paste(df$Layers, as.numeric(!df$Wrinkles == "None"), sep = "-")
df2 <- df[, c(6, 38, 39, 40)] %>%
    group_by(Wrinkles) %>%
    summarise(across(everything(), mean),
        .groups = "drop"
    ) %>%
    as.data.frame()
df3 <- df[, c(6, 38, 39, 40)] %>%
    group_by(Wrinkles) %>%
    summarise(across(everything(), sd),
        .groups = "drop"
    ) %>%
    as.data.frame()
df2 <- reshape2::melt(df2)
df3 <- reshape2::melt(df3)
df2$sd <- df3$value
df2$Wrinkles <- gsub("_", " ", df2$Wrinkles)
df2$Wrinkles <- factor(df2$Wrinkles, levels = c(
    "None", "Fold 1", "Shear 1", "Shear 2",
    "Wrinkle 1", "Wrinkle 2",
    "Wrinkle 3", "Wrinkle 4", "Wrinkle 5",
    "Wrinkle 6", "Wrinkle 7", "Wrinkle 8",
    "Wrinkle 9", "Wrinkle 10", "Wrinkle 11",
    "Wrinkle 12", "Wrinkle 13"
))
df2$variable <- factor(df2$variable, levels = c("Simulated doublet", "Artifact", "Non-artifact"))
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/Br6522_mid_SNN_artifacts.pdf", width = 13, height = 3)
ggplot(df2, aes(x = Wrinkles, y = value, fill = variable)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin = value - 0.5 * sd, ymax = value + 0.5 * sd), width = .2, position = position_dodge(.9)) +
    theme_classic() +
    xlab("") +
    ylab("Proportion of NN") +
    ggtitle("Br6522 Middle") +
    theme(legend.title = element_blank(), plot.title = element_text(face = "bold"))
dev.off()

df <- Br8667_post_all@meta.data
df <- df[!df$Type == "Simulated doublet", ]
df$Layer_Wrinkle <- paste(df$Layers, as.numeric(!df$Wrinkles == "None"), sep = "-")
df2 <- df[, c(6, 25, 26, 27)] %>%
    group_by(Wrinkles) %>%
    summarise(across(everything(), mean),
        .groups = "drop"
    ) %>%
    as.data.frame()
df3 <- df[, c(6, 25, 26, 27)] %>%
    group_by(Wrinkles) %>%
    summarise(across(everything(), sd),
        .groups = "drop"
    ) %>%
    as.data.frame()
df2 <- reshape2::melt(df2)
df3 <- reshape2::melt(df3)
df2$sd <- df3$value
df2$Wrinkles <- gsub("_", " ", df2$Wrinkles)
df2$Wrinkles <- factor(df2$Wrinkles, levels = c(
    "None", "Shear 1",
    "Wrinkle 1", "Wrinkle 2",
    "Wrinkle 3", "Wrinkle 4", "Wrinkle 5",
    "Wrinkle 6", "Wrinkle 7", "Wrinkle 8",
    "Wrinkle 9", "Wrinkle 11",
    "Wrinkle 12", "Wrinkle 13"
))
df2$variable <- factor(df2$variable, levels = c("Simulated doublet", "Artifact", "Non-artifact"))
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/Br8667_post_SNN_artifacts.pdf", width = 13, height = 3)
ggplot(df2, aes(x = Wrinkles, y = value, fill = variable)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin = value - 0.5 * sd, ymax = value + 0.5 * sd), width = .2, position = position_dodge(.9)) +
    theme_classic() +
    xlab("") +
    ylab("Proportion of NN") +
    ggtitle("Br8667 Posterior") +
    theme(legend.title = element_blank(), plot.title = element_text(face = "bold"))
dev.off()


Br6522_ant_all <- RunUMAP(Br6522_ant_all, dims = 1:30)
Idents(Br6522_ant_all) <- Br6522_ant_all$Layers
L1_cells <- WhichCells(Br6522_ant_all, idents = "Layer 1")
L5_cells <- WhichCells(Br6522_ant_all, idents = "Layer 5")
L_1_5_cells <- WhichCells(Br6522_ant_all, idents = "Layer 1/Layer 5")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/Br6522_ant_simulated.pdf", width = 8, height = 4)
DimPlot(Br6522_ant_all,
    split.by = "Type", cells.highlight = list("Layer 1" = L1_cells, "Layer 5" = L5_cells, "Layer 1/Layer 5" = L_1_5_cells),
    cols.highlight = c("blue", "purple", "red"), cols = "grey"
) + ggtitle("Br6522 Anterior") + xlab("UMAP 1") + ylab("UMAP 2")
dev.off()

Br6522_mid_all <- RunUMAP(Br6522_mid_all, dims = 1:30)
Idents(Br6522_mid_all) <- Br6522_mid_all$Layers
L1_cells <- WhichCells(Br6522_mid_all, idents = "Layer 1")
L5_cells <- WhichCells(Br6522_mid_all, idents = "Layer 5")
L_1_5_cells <- WhichCells(Br6522_mid_all, idents = "Layer 1/Layer 5")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/Br6522_mid_simulated.pdf", width = 8, height = 4)
DimPlot(Br6522_mid_all,
    split.by = "Type", cells.highlight = list("Layer 1" = L1_cells, "Layer 5" = L5_cells, "Layer 1/Layer 5" = L_1_5_cells),
    cols.highlight = c("blue", "purple", "red"), cols = "grey"
) + ggtitle("Br6522 Middle") + xlab("UMAP 1") + ylab("UMAP 2")
dev.off()

Br8667_post_all <- RunUMAP(Br8667_post_all, dims = 1:30)
Idents(Br8667_post_all) <- Br8667_post_all$Layers
L2_cells <- WhichCells(Br8667_post_all, idents = "Layer 2")
L6_cells <- WhichCells(Br8667_post_all, idents = "Layer 6")
L_2_6_cells <- WhichCells(Br8667_post_all, idents = "Layer 2/Layer 6")
pdf("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/Br8667_post_simulated.pdf", width = 8, height = 4)
DimPlot(Br8667_post_all,
    split.by = "Type", cells.highlight = list("Layer 2" = L2_cells, "Layer 6" = L6_cells, "Layer 2/Layer 6" = L_2_6_cells),
    cols.highlight = c("blue", "purple", "red"), cols = "grey"
) + ggtitle("Br8667 Posterior") + xlab("UMAP 1") + ylab("UMAP 2")
dev.off()
