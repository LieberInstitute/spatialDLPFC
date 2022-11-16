library(spatialLIBD)
library(pheatmap)
library(Polychrome)
library(tidyverse)

# Load modeling results
load(file = here::here("processed-data", "rdata", "spe", "07_layer_differential_expression", "parsed_modeling_results_k9.Rdata"))
modeling_results_9 <- modeling_results$enrichment
load(file = here::here("processed-data", "rdata", "spe", "07_layer_differential_expression", "parsed_modeling_results_k16.Rdata"))
modeling_results_16 <- modeling_results$enrichment

table(modeling_results_9$t_stat_1 > 0, modeling_results_9$fdr_1 < 0.05)

num_genes <- c(
    nrow(modeling_results_9[modeling_results_9$t_stat_1 > 0 & modeling_results_9$fdr_1 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_2 > 0 & modeling_results_9$fdr_2 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_3 > 0 & modeling_results_9$fdr_3 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_4 > 0 & modeling_results_9$fdr_4 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_5 > 0 & modeling_results_9$fdr_5 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_6 > 0 & modeling_results_9$fdr_6 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_7 > 0 & modeling_results_9$fdr_7 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_8 > 0 & modeling_results_9$fdr_8 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_9 > 0 & modeling_results_9$fdr_9 < 0.05, ])
)

# pivot long, filter with dypler, count
long <- modeling_results_9 %>%
    select(ensembl, starts_with("t_stat")) %>%
    pivot_longer(!ensembl)
head(long)

# do same for fdr and then left join them
num_genes <- as.data.frame(num_genes)
colors_bayesSpace <- Polychrome::palette36.colors(9)
num_genes$colors <- as.character(colors_bayesSpace)
num_genes$colors <- as.factor(num_genes$colors)
names(num_genes$colors) <- c(1:9)
num_genes$num_genes <- as.numeric(num_genes$num_genes)
num_genes$cluster <- as.factor(c(1:9))
df <- num_genes

pdf(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/07_layer_differential_expression/num_enriched_genes_per_cluster_bar.pdf")
ggplot(df, aes(x = cluster, y = num_genes)) +
    geom_bar(stat = "identity", color = df$colors, fill = df$colors) +
    scale_color_manual(values = df$colors) +
    theme_bw() +
    xlab("BayesSpace Cluster") +
    ylab("Number of Enriched Genes")
dev.off()

#####
num_genes <- c(
    nrow(modeling_results_9[modeling_results_9$t_stat_1 < 0 & modeling_results_9$fdr_1 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_2 < 0 & modeling_results_9$fdr_2 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_3 < 0 & modeling_results_9$fdr_3 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_4 < 0 & modeling_results_9$fdr_4 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_5 < 0 & modeling_results_9$fdr_5 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_6 < 0 & modeling_results_9$fdr_6 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_7 < 0 & modeling_results_9$fdr_7 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_8 < 0 & modeling_results_9$fdr_8 < 0.05, ]),
    nrow(modeling_results_9[modeling_results_9$t_stat_9 < 0 & modeling_results_9$fdr_9 < 0.05, ])
)

# pivot long, filter with dypler, count
long <- modeling_results_9 %>%
    select(ensembl, starts_with("t_stat")) %>%
    pivot_longer(!ensembl)
head(long)

# do same for fdr and then left join them
num_genes <- as.data.frame(num_genes)
colors_bayesSpace <- Polychrome::palette36.colors(9)
num_genes$colors <- as.character(colors_bayesSpace)
num_genes$colors <- as.factor(num_genes$colors)
names(num_genes$colors) <- c(1:9)
num_genes$num_genes <- as.numeric(num_genes$num_genes)
num_genes$cluster <- as.factor(c(1:9))
df <- num_genes

pdf(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/07_layer_differential_expression/num_depleted_genes_per_cluster_bar.pdf")
ggplot(df, aes(x = cluster, y = num_genes)) +
    geom_bar(stat = "identity", fill = df$colors) +
    scale_color_manual(values = df$colors) +
    theme_bw() +
    xlab("BayesSpace Cluster") +
    ylab("Number of Depleted Genes")
dev.off()
