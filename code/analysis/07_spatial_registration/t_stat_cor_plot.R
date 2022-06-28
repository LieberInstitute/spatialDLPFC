library(spatialLIBD)
library(ggplot2)

k = 7
load(file = here::here("processed-data","rdata","spe","08_layer_differential_expression",paste0("parsed_modeling_results_k",k,".Rdata")))
stats <-modeling_results$enrichment
rownames(stats) = stats$ensembl
dim(stats)
#[1] 12327    23
# colnames(stats)
# [1] "t_stat_1"  "t_stat_2"  "t_stat_3"  "t_stat_4"  "t_stat_5"  "t_stat_6" 
# [7] "t_stat_7"  "p_value_1" "p_value_2" "p_value_3" "p_value_4" "p_value_5"
# [13] "p_value_6" "p_value_7" "fdr_1"     "fdr_2"     "fdr_3"     "fdr_4"    
# [19] "fdr_5"     "fdr_6"     "fdr_7"     "ensembl"   "gene"  

tstats <- spatialLIBD::fetch_data("modeling_results")
tstats <-tstats$enrichment
rownames(tstats) = tstats$ensembl
dim(tstats)
#[1] 22331    23
colnames(tstats)
# [1] "t_stat_WM"      "t_stat_Layer1"  "t_stat_Layer2"  "t_stat_Layer3" 
# [5] "t_stat_Layer4"  "t_stat_Layer5"  "t_stat_Layer6"  "p_value_WM"    
# [9] "p_value_Layer1" "p_value_Layer2" "p_value_Layer3" "p_value_Layer4"
# [13] "p_value_Layer5" "p_value_Layer6" "fdr_WM"         "fdr_Layer1"    
# [17] "fdr_Layer2"     "fdr_Layer3"     "fdr_Layer4"     "fdr_Layer5"    
# [21] "fdr_Layer6"     "ensembl"        "gene"  

tstats <-
  tstats[, grep("[f|t]_stat_", colnames(tstats))]
colnames(tstats) <-
  gsub("[f|t]_stat_", "", colnames(tstats))
colnames(tstats)
#[1] "WM"     "Layer1" "Layer2" "Layer3" "Layer4" "Layer5" "Layer6"

stats <-
  stats[, grep("[f|t]_stat_", colnames(stats))]
colnames(stats) <-
  gsub("[f|t]_stat_", "", colnames(stats))
colnames(stats)
# [1] "1" "2" "3" "4" "5" "6" "7"

#from layer_stat_cor code https://github.com/LieberInstitute/spatialLIBD/blob/d44eea67ffd876ace89718024bf855dbce268111/R/layer_stat_cor.R#L77-L79
# top_n_index <- unique(as.vector(apply(tstats, 2, function(t) {
#   order(t, decreasing = TRUE)[seq_len(100)]
# })))
# length(top_n_index)
# #[1] 692
# summary(top_n_index)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 45    4935   10384   10574   15708   22225 

top_n_index <- data.frame(index = integer(),layer = factor())
for (i in 1:k){
  idx <- as.vector(order(tstats[,i], decreasing = TRUE)[seq_len(100)])
  lyr <- as.vector(rep(colnames(tstats)[i],100))
  df <- data.frame(idx,lyr)
  top_n_index <- rbind(top_n_index,df)
}
dim(top_n_index)
# [1] 700   2

top_n_index <- top_n_index[!duplicated(top_n_index$idx),]
dim(top_n_index)
# [1] 692   2

tstats_small <- tstats[top_n_index$idx, , drop = FALSE]
dim(tstats_small)
#[1] 692   7

tstats_small <- cbind(tstats_small,top_n_index$lyr)

mm<-match(rownames(tstats_small),rownames(stats))
length(which(!is.na(mm)))
#[1] 584

tstats_small <- tstats_small[!is.na(mm), ]
dim(tstats_small)
# [1] 584   8

stats_small <-stats[mm[!is.na(mm)], ]
dim(stats_small)
#[1] 584   7

top_genes <- rownames(tstats_small)
length(top_genes)
#[1] 584

#make plot 
pdf(file = here::here("plots","07_spatial_registration","t_cor_k7_wm.pdf"))
plot(stats_small[,7],tstats_small$WM)
dev.off()


#make nicer plot 
dat_small <- cbind(stats_small,tstats_small)
colnames(dat_small)
colnames(dat_small)[15] <- "Layer"

#from Louise's code https://github.com/LieberInstitute/goesHyde_mdd_rnaseq/blob/2bb13a25fad8d1260ef38a8f073be4387c1f9ed0/differential_expression/code/utils.R#L51-L59
pdf(file = here::here("plots","07_spatial_registration","ggplot_t_cor_k7_wm_colored.pdf"))
ggplot(dat_small, aes(x = dat_small[,7], y = WM, color = Layer)) +
  geom_point(alpha = 0.7, size = 0.8) +
  scale_color_manual(values = c("WM" = "#1A1A1A",
                                "Layer1"="#F0027F",
                                "Layer2"="#377EB8",
                                "Layer3" = "#4DAF4A",
                                "Layer4" = "#984EA3",
                                "Layer5" = "#FFD700",
                                "Layer6" = "#FF7F00")) +
  labs(x = "t-stats BayesSpace Cluster 7", 
       y = "t-stats Manual Annotations WM",
       title = "t-stat Correlation") +
  theme_bw() +
  geom_smooth(method=lm, se=FALSE, colour = "#00BA38") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20), axis.text = element_text(size = 25))
dev.off()


#Make cor plot with all genes and color by just the 
common_genes <- intersect(rownames(stats),rownames(tstats))
length(common_genes)
# [1] 12150
stats <- stats[common_genes,]
tstats <- tstats[common_genes,]

pdf(file = here::here("plots","07_spatial_registration","t_cor_k7_wm_all_genes.pdf"))
plot(stats[,7],tstats$WM)
dev.off()

dat <- cbind(stats,tstats)
dat$TopGene<- "FALSE"
dat[top_genes,]$TopGene <-"TRUE"

pdf(file = here::here("plots","07_spatial_registration","ggplot_t_cor_k7_wm_all_genes.pdf"))
ggplot(dat, aes(x = dat[,7], y = WM,color = TopGene)) +
  geom_point(alpha = 0.7, size = 0.5) +
  labs(x = "t-stats BayesSpace Cluster 7", 
       y = "t-stats Ground Truth WM",
       title = "t-stat Correlation") +
  theme_bw() +
  geom_smooth(method=lm, se=FALSE, colour = "#00BA38") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c ("TRUE" = "red", "FALSE" = "gray"))+
  guides(fill=guide_legend(title="Top Gene"))
dev.off()

