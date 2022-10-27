library("tidyverse")
library("here")
library("sessioninfo")
library("patchwork")

## setup output dirs
# # data_dir <- here("processed-data", "rdata", "spe", "13_nnSVG", "02_compile_nnSVG_output")
# if(!dir.exists(data_dir)) dir.create(data_dir)

plot_dir <- here("plots", "13_nnSVG", "03_summarize_nnSVG_output")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

#### Explore Data ####

## Load all
load(here("processed-data", "rdata", "spe", "13_nnSVG", "02_compile_nnSVG_output", "nnSVG_all.Rdata"))

nnSVG_all  |> filter(FDR < 0.05) |> count()

nnSVG_all_summary <- nnSVG_all |>
  group_by(domains, gene_id, gene_name) |> 
  summarize(n = n(), 
            n_signif = sum(FDR < 0.05),
            max_FDR = max(FDR),
            min_FDR = min(FDR),
            min_rank = min(rank),
            max_rank = max(rank),
            mean_rank = mean(rank), 
            median_rank = median(rank),
            top100_rank = sum(rank <= 100)) |>
  ungroup()  |>
  group_by(domains)

nnSVG_all_summary |> filter(n == 30) |> count()
# domains     n
# <chr>   <int>
# 1 4v16      465
# 2 5v9       568
nnSVG_all_summary |> filter(max_FDR !=1) |> ungroup()  |> count(n ==1)

nnSVG_all_summary |> filter(max_FDR !=1, n !=1)

nnSVG_all_summary |> arrange(-top100_rank) |> slice(1:3)

## lowest max rank for each domains
low_max <- nnSVG_all_summary |> filter(n > 20) |> arrange(max_rank) |> slice(1)

n_gene_distribution <- nnSVG_all_summary |>
  ggplot(aes(n, color = domains)) +
  geom_density()

ggsave(n_gene_distribution, filename = here(plot_dir, "n_gene_distribution.png"), width = 10)

## What is the FDR distribution like for low max?

low_max <- low_max |> left_join(nnSVG_all)

low_max |> count()

low_max_FDR <- low_max |>
  ggplot(aes(x = domains , y = FDR, fill = domains)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25) +
  geom_jitter(aes(color = domains)) +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed")

ggsave(low_max_FDR, filename = here(plot_dir, "low_max_FDR.png"))

## Mean rank vs. n
mean_rank_n_box <- nnSVG_all_summary |>
  ggplot(aes(x = as.factor(n), y = mean_rank, fill = domains)) +
  geom_boxplot() +
  facet_wrap(~domains) +
  theme(legend.position = "None")


ggsave(mean_rank_n_box, filename = here(plot_dir, "mean_rank_n_box.png"), width = 12)

mean_rank_top_n_scatter <- nnSVG_all_summary |> 
  arrange(mean_rank) |>
  slice(1:500) |>
  # ggplot(aes(x = as.factor(n), y = mean_rank, fill = domains)) +
  ggplot(aes(x = n, y = mean_rank, color = top100_rank)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~domains) 

ggsave(mean_rank_top_n_scatter, filename = here(plot_dir, "mean_rank_top_n_scatter.png"), width = 12)

nnSVG_all_summary |> 
  arrange(mean_rank) |>
  slice(1:5)

## n rank > 100 vs. mean
library("ggrepel")

mean_rank_n_top100 <- nnSVG_all_summary |>
  # slice(1:100) |>
  ggplot(aes(x = top100_rank, y = mean_rank)) +
  geom_point(alpha = 0.2) +
  geom_text_repel(aes(label = ifelse(top100_rank > 10, gene_name, NA)), size = 2) +
  # geom_jitter(alpha = 0.2) +
  facet_wrap(~domains) +
  theme(legend.position = "None")

ggsave(mean_rank_n_top100, filename = here(plot_dir, "mean_rank_n_top100.png"), width = 12)

mean_rank_n_top100 <- nnSVG_all_summary |>
  arrange(-top100_rank) |>
  slice(1:100) |>
  ggplot(aes(x = top100_rank, y = mean_rank)) +
  geom_point(alpha = 0.2) +
  geom_text_repel(aes(label = gene_name), size = 2) +
  # geom_text(aes(label = gene_name), size = 2) +
  # geom_jitter(alpha = 0.2) +
  facet_wrap(~domains) +
  theme(legend.position = "None")

ggsave(mean_rank_n_top100, filename = here(plot_dir, "mean_rank_n_top100_filter.png"), width = 12)




