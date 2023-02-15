rm(list = ls())
library(data.table)
library(corTest)
library(corrplot)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
# Read cell scores
sample <- "Br6522_ant"
datDir <- "/data/abattle4/prashanthi/dewrinkler/data/tangram_results/"
map <- fread(paste0(datDir, sample, "/weights.csv"))
spot_meta <- read.csv(paste0(datDir, sample, "/spot_meta.csv"))
cell_meta <- read.csv(paste0(datDir, sample, "/cell_meta.csv"))
cellTypes <- unique(cell_meta$cellType_broad_hc)
map <- as.matrix(map)
map_cellType <- matrix(NA, nrow = dim(spot_meta)[1], ncol = length(cellTypes))

for(i in c(1:dim(map)[2])){
  weights_df <- data.frame(map[ ,i], cell_meta$cellType_broad_hc)
  colnames(weights_df) <- c("spot_weight", "cell_type")
  weights_df$cell_type <- factor(weights_df$cell_type, levels = cellTypes)
  summary_df <- aggregate(weights_df$spot_weight, by=list(Category=weights_df$cell_type), FUN=sum)
  map_cellType[i, ] <- summary_df$x
}

colnames(map_cellType) <- cellTypes
map_cellType <- data.frame(map_cellType)
map_cellType$sumCells <- rowSums(map_cellType[ ,1:7])
map_cellType$Layer <- spot_meta$LayerAnnot
map_cellType$Wrinkle <- spot_meta$WrinkleAnnot
map_cellType <- map_cellType[!map_cellType$Layer == "Unknown", ]
if(sample == "Br6522_ant"){
map_cellType$Wrinkle <- factor(map_cellType$Wrinkle, levels = c("None", "Fold_1", "Shear_1", "Shear_2", 
                                                                "Shear_3", "Wrinkle_1", "Wrinkle_2", "Wrinkle_3", 
                                                                "Wrinkle_4", "Wrinkle_5", "Wrinkle_6", "Wrinkle_7", 
                                                                "Wrinkle_8", "Wrinkle_9"))}
if(sample == "Br6522_mid"){
  map_cellType$Wrinkle <- factor(map_cellType$Wrinkle, levels = c("None", "Fold_1", "Shear_1", "Shear_2", "Wrinkle_1", 
                                                                  "Wrinkle_2", "Wrinkle_3", "Wrinkle_4", "Wrinkle_5", 
                                                                  "Wrinkle_6", "Wrinkle_7", 
                                                                  "Wrinkle_8", "Wrinkle_9", 
                                                                  "Wrinkle_10", "Wrinkle_11", "Wrinkle_12", "Wrinkle_13"))
  
}
if(sample == "Br8667_post"){
  map_cellType$Wrinkle <- factor(map_cellType$Wrinkle, levels = c("None", "Shear_1", "Wrinkle_1", 
                                                                  "Wrinkle_2", "Wrinkle_3", "Wrinkle_4", "Wrinkle_5", 
                                                                  "Wrinkle_6", "Wrinkle_7", 
                                                                  "Wrinkle_8", "Wrinkle_9",
                                                                  "Wrinkle_11", "Wrinkle_12", "Wrinkle_13"))
  
}
ggplot(map_cellType, aes(x = Layer, y = sumCells, fill = Wrinkle)) + geom_boxplot() + theme_classic()

diff_corr <- function(none, all){
  p_value_mat <- matrix(NA, nrow = 7, ncol = 7)
  for(i in c(1:7)){
    for(j in c(1:7)){
      if(i == j){
        p_value_mat[i,j] <- NA
      }else{
        res <- st1(all[ ,i], all[ ,j], none[ ,i], none[ ,j])
        p_value_mat[i, j] <- res$pval
      }
    }
  }
  colnames(p_value_mat) <- colnames(none)[1:7]
  rownames(p_value_mat) <- colnames(p_value_mat)
  p_value_mat[lower.tri(p_value_mat)] <- NA
  p_values_df <- melt(p_value_mat)
  colnames(p_values_df) <- c("cell_type_1", "cell_type_2", "p_value")
  p_values_df <- p_values_df[!is.na(p_values_df$p_value), ]
  # sig_mat <- p_value_mat
  # sig_mat[p_value_mat > 0.1] <- "n.s"
  # sig_mat[p_value_mat < 0.1 & p_value_mat > 0.05] <- "*"
  # sig_mat[p_value_mat < 0.05 & p_value_mat > 0.01] <- "**"
  # sig_mat[p_value_mat < 0.01] <- "***"
  # sig_mat[lower.tri(sig_mat)] <- ""
  # diag(sig_mat) <- ""
  # p_value_mat[lower.tri(p_value_mat)] <- NA
  # pheatmap(p_value_mat, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = sig_mat, 
  #          color = colorRampPalette(rev(brewer.pal(n = 7, name ="Reds")))(100), breaks= seq(0, 1, 0.01),
  #          main = title, border_color = NA, na_col = "transparent")
  p_values_df
}

if(!sample == "Br8667_post"){
L1 <- map_cellType[map_cellType$Layer == "Layer 1", ]
L1_mat <- as.matrix(L1[ ,1:7])
L1_mat <- scale(L1_mat, center = TRUE, scale = FALSE)
L1_pca <- svd(L1_mat)
L1_u <- L1_pca$u
colnames(L1_u) <- paste0("U", c(1:dim(L1_u)[2]))
L1 <- cbind(L1, L1_u)
L1$Wrinkle <- as.factor(L1$Wrinkle)
ggplot(L1, aes(x = U1, y = U2, colour = Wrinkle)) + geom_point() + theme_classic()
ggplot(L1, aes(x = U1, y = U2, colour = sumCells)) + geom_point() + theme_classic()
L1_none <- L1[L1$Wrinkle == "None", ]
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample ,"_layer1.pdf"), width = 5, height = 5)
corrplot(cor(L1[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample ,"_layer1_excl_artifact.pdf"), width = 5, height = 5)
corrplot(cor(L1_none[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
L1_pvalue <- diff_corr(L1_none, L1)
L1_pvalue$Layers <- "Layer 1"
}

L2 <- map_cellType[map_cellType$Layer == "Layer 2", ]
L2_mat <- as.matrix(L2[ ,1:7])
L2_mat <- scale(L2_mat, center = TRUE, scale = FALSE)
L2_pca <- svd(L2_mat)
L2_u <- L2_pca$u
colnames(L2_u) <- paste0("U", c(1:dim(L2_u)[2]))
L2 <- cbind(L2, L2_u)
L2$Wrinkle <- as.factor(L2$Wrinkle)
ggplot(L2, aes(x = U1, y = U2, colour = Wrinkle)) + geom_point() + theme_classic()
ggplot(L2, aes(x = U1, y = U2, colour = sumCells)) + geom_point() + theme_classic()
L2_none <- L2[L2$Wrinkle == "None", ]
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample ,"_layer2.pdf"), width = 5, height = 5)
corrplot(cor(L2[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer2_excl_artifact.pdf"), width = 5, height = 5)
corrplot(cor(L2_none[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
L2_pvalue <- diff_corr(L2_none, L2)
L2_pvalue$Layers <- "Layer 2"

L3 <- map_cellType[map_cellType$Layer == "Layer 3", ]
L3_mat <- as.matrix(L3[ ,1:7])
L3_mat <- scale(L3_mat, center = TRUE, scale = FALSE)
L3_pca <- svd(L3_mat)
L3_u <- L3_pca$u
colnames(L3_u) <- paste0("U", c(1:dim(L3_u)[2]))
L3 <- cbind(L3, L3_u)
L3$Wrinkle <- as.factor(L3$Wrinkle)
ggplot(L3, aes(x = U1, y = U2, colour = Wrinkle)) + geom_point() + theme_classic()
ggplot(L3, aes(x = U1, y = U2, colour = sumCells)) + geom_point() + theme_classic()
L3_none <- L3[L3$Wrinkle == "None", ]
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer3.pdf"), width = 5, height = 5)
corrplot(cor(L3[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer3_excl_artifact.pdf"), width = 5, height = 5)
corrplot(cor(L3_none[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
L3_pvalue <- diff_corr(L3_none, L3)
L3_pvalue$Layers <- "Layer 3"

L4 <- map_cellType[map_cellType$Layer == "Layer 4", ]
L4_mat <- as.matrix(L4[ ,1:7])
L4_mat <- scale(L4_mat, center = TRUE, scale = FALSE)
L4_pca <- svd(L4_mat)
L4_u <- L4_pca$u
colnames(L4_u) <- paste0("U", c(1:dim(L4_u)[2]))
L4 <- cbind(L4, L4_u)
L4$Wrinkle <- as.factor(L4$Wrinkle)
ggplot(L4, aes(x = U1, y = U2, colour = Wrinkle)) + geom_point() + theme_classic()
ggplot(L4, aes(x = U1, y = U2, colour = sumCells)) + geom_point() + theme_classic()
L4_none <- L4[L4$Wrinkle == "None", ]
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer4.pdf"), width = 5, height = 5)
corrplot(cor(L4[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer4_excl_artifact.pdf"), width = 5, height = 5)
corrplot(cor(L4_none[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
L4_pvalue <- diff_corr(L4_none, L4)
L4_pvalue$Layers <- "Layer 4"

L5 <- map_cellType[map_cellType$Layer == "Layer 5", ]
L5_mat <- as.matrix(L5[ ,1:7])
L5_mat <- scale(L5_mat, center = TRUE, scale = FALSE)
L5_pca <- svd(L5_mat)
L5_u <- L5_pca$u
colnames(L5_u) <- paste0("U", c(1:dim(L5_u)[2]))
L5 <- cbind(L5, L5_u)
L5$Wrinkle <- as.factor(L5$Wrinkle)
ggplot(L5, aes(x = U1, y = U2, colour = Wrinkle)) + geom_point() + theme_classic()
ggplot(L5, aes(x = U1, y = U2, colour = sumCells)) + geom_point() + theme_classic()
L5_none <- L5[L5$Wrinkle == "None", ]
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer5.pdf"), width = 5, height = 5)
corrplot(cor(L5[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer5_excl_artifact.pdf"), width = 5, height = 5)
corrplot(cor(L5_none[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
L5_pvalue <- diff_corr(L5_none, L5)
L5_pvalue$Layers <- "Layer 5"

L6 <- map_cellType[map_cellType$Layer == "Layer 6", ]
L6_mat <- as.matrix(L6[ ,1:7])
L6_mat <- scale(L6_mat, center = TRUE, scale = FALSE)
L6_pca <- svd(L6_mat)
L6_u <- L6_pca$u
colnames(L6_u) <- paste0("U", c(1:dim(L6_u)[2]))
L6 <- cbind(L6, L6_u)
L6$Wrinkle <- as.factor(L6$Wrinkle)
ggplot(L6, aes(x = U1, y = U2, colour = Wrinkle)) + geom_point() + theme_classic()
ggplot(L6, aes(x = U1, y = U2, colour = sumCells)) + geom_point() + theme_classic()
L6_none <- L6[L6$Wrinkle == "None", ]
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer6.pdf"), width = 5, height = 5)
corrplot(cor(L6[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer6_excl_artifact.pdf"), width = 5, height = 5)
corrplot(cor(L6_none[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
L6_pvalue <- diff_corr(L6_none, L6)
L6_pvalue$Layers <- "Layer 6"

if(!sample == "Br8667_post"){
WM <- map_cellType[map_cellType$Layer == "WM", ]
WM_mat <- as.matrix(WM[ ,1:7])
WM_mat <- scale(WM_mat, center = TRUE, scale = FALSE)
WM_pca <- svd(WM_mat)
WM_u <- WM_pca$u
colnames(WM_u) <- paste0("U", c(1:dim(WM_u)[2]))
WM <- cbind(WM, WM_u)
WM$Wrinkle <- as.factor(WM$Wrinkle)
ggplot(WM, aes(x = U1, y = U2, colour = Wrinkle)) + geom_point() + theme_classic()
ggplot(WM, aes(x = U1, y = U2, colour = sumCells)) + geom_point() + theme_classic()
WM_none <- WM[WM$Wrinkle == "None", ]
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_WM.pdf"), width = 5, height = 5)
corrplot(cor(WM[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_WM_excl_artifact.pdf"), width = 5, height = 5)
corrplot(cor(WM_none[ ,1:7]), method = 'color', diag = TRUE, type = 'upper', addCoef.col = 'black',
         order = 'original', tl.col = "black", number.digits = 2, number.cex =0.8)
dev.off()
WM_pvalue <- diff_corr(WM_none, WM)
WM_pvalue$Layers <- "WM"}

if(!sample == "Br8667_post"){
  pvalue_df <- rbind(L1_pvalue, L2_pvalue, L3_pvalue, L4_pvalue, L5_pvalue, L6_pvalue, WM_pvalue)
}else{
  pvalue_df <- rbind(L2_pvalue, L3_pvalue, L4_pvalue, L5_pvalue, L6_pvalue)
}

pvalue_df$corr_p_values <- p.adjust(pvalue_df$p_value, method = "fdr")

plot_heatmap <- function(pvalue_df, layer){
  p_value_mat <-dcast(pvalue_df[pvalue_df$Layers == layer, c(1, 2, 5)], cell_type_1 ~ cell_type_2)
  rownames(p_value_mat) <- p_value_mat$cell_type_1
  p_value_mat$cell_type_1 <- NULL
  p_value_mat <- as.matrix(p_value_mat)
  sig_mat <- p_value_mat
  sig_mat[p_value_mat > 0.1] <- "n.s"
  sig_mat[p_value_mat < 0.1 & p_value_mat > 0.05] <- "*"
  sig_mat[p_value_mat < 0.05 & p_value_mat > 0.01] <- "**"
  sig_mat[p_value_mat < 0.01] <- "***"
  sig_mat[is.na(sig_mat)] <- ""
  pheatmap(p_value_mat, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = sig_mat,
            color = colorRampPalette(rev(brewer.pal(n = 7, name ="Reds")))(100), breaks= seq(0, 1, 0.01),
            main = layer, border_color = NA, na_col = "transparent", number_color = "black", fontsize = 14)
}

if(!sample == "Br8667_post"){
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer1_pval.pdf"), width = 4, height = 4)
plot_heatmap(pvalue_df, "Layer 1")
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_WM_pval.pdf"), width = 4, height = 4)
plot_heatmap(pvalue_df, "WM")
dev.off()}

pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer2_pval.pdf"), width = 4, height = 4)
plot_heatmap(pvalue_df, "Layer 2")
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer3_pval.pdf"), width = 4, height = 4)
plot_heatmap(pvalue_df, "Layer 3")
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer4_pval.pdf"), width = 4, height = 4)
plot_heatmap(pvalue_df, "Layer 4")
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer5_pval.pdf"), width = 4, height = 4)
plot_heatmap(pvalue_df, "Layer 5")
dev.off()
pdf(paste0("/data/abattle4/prashanthi/dewrinkler/figures/fig_S6/", sample, "_layer6_pval.pdf"), width = 4, height = 4)
plot_heatmap(pvalue_df, "Layer 6")
dev.off()

pvalue_df$sample <- sample
pvalue_df <- pvalue_df[ ,c("sample", "Layers", "cell_type_1", "cell_type_2", "p_value", "corr_p_values")]
saveRDS(pvalue_df, paste0("/data/abattle4/prashanthi/dewrinkler/tables/", sample, "_pvalue_df.rds"))

    
    