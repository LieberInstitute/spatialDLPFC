
library("tidyverse")
library("here")
library("sessioninfo")

## setup output dirs
data_dir <- here("processed-data", "rdata", "spe", "13_nnSVG", "02_combine_output")
if(!dir.exists(data_dir)) dir.create(data_dir)

plot_dir <- here("plots", "13_nnSVG", "02_combine_output")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

#### Progress ####
output_dir <- here("processed-data", "rdata", "spe", "13_nnSVG", "01_nnSVG_pairwise")
output_fn <- list.files(output_dir, pattern = "nnSVG_k16", full.names = TRUE)
names(output_fn) <- gsub("nnSVG_k16-|.RData","", output_fn)

length(output_fn)
# [1] 169

output_t <- tibble(file = basename(output_fn)) |>
            mutate(label = gsub("nnSVG_k16-|.RData","", file)) |>
  separate(label, into = c("domains","Sample"), sep = "-")

Sys.time()
(n_done <- output_t |> count(domains))
# [1] "2022-10-24 10:25:15 EDT"
# domains     n
# <chr>   <int>
# 1 12v13      27
# 2 12v16      26
# 3 4v16       30
# 4 5v9        30
# 5 7v12       28
# 6 7v13       28


#### Get info from Log files ####
## function to overcome incomplete logs 
get_log_info <- function(logs, pattern, add = 0){
  info <- map(logs, ~(.x[grep(pattern,.x) + add]))
  info <- modify_if(info, ~ length(.) == 0,  ~ NA_character_)
  return(info)
}

## Read log files 
log_dir <- here("code", "analysis", "13_nnSVG", "logs")
log_fn <- list.files(log_dir, pattern = "01_nnSVG_pairwise_loop_", full.names = TRUE)
names(log_fn) <- gsub("01_nnSVG_pairwise_loop_", "nnSVG_", gsub(".txt","",basename(log_fn)))

length(log_fn)
# [1] 178

logs <- map(log_fn, readLines)

## Job completed
jobs_done <- map_lgl(logs, ~any(grepl("Job ends", .x)))
table(jobs_done)
# FALSE  TRUE 
# 2   176 

n_spots <- get_log_info(logs, "ncol")
error_spot<- grepl("Error",n_spots)
n_spots <- ifelse(error_spot, NA, parse_number(unlist(n_spots)))

error <- grepl("Error", logs)

log_t <- tibble(log = names(log_fn),
                Sample = gsub(".*(Br\\d+_.*),.*","\\1",get_log_info(logs, "Running sample: ")),
                n_spots =  n_spots,
                n_genes = parse_number(unlist(get_log_info(logs, "nrow"))),
                error_spot = error_spot,
                error_anny = error,
                runtime = unlist(get_log_info(logs, "Reproducibility information:", 3)),
                done = jobs_done) |> 
  separate(runtime, into = c("user", "system", "elapsed"), sep = " +", convert = TRUE, extra = "drop") |>
  separate(log, into = c(NA, "domains", "sample_i"), convert = TRUE, remove = FALSE) |>
  left_join(output_t)

Sys.time()
log_t |> count(domains, done)

# [1] "2022-10-24 10:54:41 EDT"
# domains done      n
# <chr>   <lgl> <int>
# 1 12v13   TRUE     30
# 2 12v16   FALSE     1
# 3 12v16   TRUE     29
# 4 4v16    TRUE     30
# 5 5v9     TRUE     30
# 6 7v12    FALSE     1
# 7 7v12    TRUE     29
# 8 7v13    TRUE     28

log_t |> filter(error) |> count(domains)
# domains     n
# 1 12v13       3
# 2 12v16       3
# 3 7v12        1

## Checkout runtime
runtime_long <- log_t |>
  filter(done) |> 
  mutate(prod = n_spots * n_genes) |>
  select(domains, Sample, n_spots, n_genes, elapsed, user, prod) |>
  pivot_longer(!c(domains, Sample, elapsed, user), names_to = "metric", values_to = "n")

runtime_plot <- runtime_long |>
  ggplot(aes(x = n, y = elapsed, color = domains)) +
  geom_point() +
  facet_wrap(~metric, scales = "free_x")

ggsave(runtime_plot, filename = here(plot_dir, "runtime_scatter.png"), width = 10)

## cells vs. genes
spots_v_gene_scatter <- log_t |>
  ggplot(aes(x= n_spots, y = n_genes, color = domains)) +
  geom_point() +
  geom_vline(xintercept = 65, color = "red", linetype = "dashed")

ggsave(spots_v_gene_scatter, filename = here(plot_dir, "spots_v_gene_scatter.png"))

spots_box <- log_t |>
  ggplot(aes(x= domains, y = n_spots, fill = domains)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25) +
  geom_jitter(aes(color = domains)) +
  geom_hline(yintercept = 65, color = "red", linetype = "dashed")

ggsave(spots_box, filename = here(plot_dir, "spots_box.png"))

genes_box <- log_t |>
  ggplot(aes(x= domains, y = n_genes, fill = domains)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25) +
  geom_jitter(aes(color = domains)) 

ggsave(genes_box, filename = here(plot_dir, "genes_box.png"))


library(patchwork)

ggsave(spots_v_gene_scatter + theme(legend.position = "None") +
         genes_box /
         (spots_box + coord_flip() + theme(legend.position = "None")), 
       filename = here(plot_dir, "spots_genes_summary.png"), width = 10)

## low genes all from one sample
log_t |> filter(n_genes < 1000)
# log           domains sample_i Sample     n_spots n_genes error  user system elapsed done  file                       
# <chr>         <chr>      <int> <chr>        <dbl>   <dbl> <lgl> <dbl>  <dbl>   <dbl> <lgl> <chr>                      
# 1 nnSVG_4v16.13 4v16          13 Br2720_ant    1390     521 FALSE  574.   5.59    582. TRUE  nnSVG_k16-4v16-Br2720_ant.…
# 2 nnSVG_5v9.13  5v9           13 Br2720_ant    1401     710 FALSE 1127.  12.6    1143. TRUE  nnSVG_k16-5v9-Br2720_ant.R…
# 3 nnSVG_7v12.13 7v12          13 Br2720_ant     873     459 FALSE 1430.  74.4    1507. TRUE  nnSVG_k16-7v12-Br2720_ant.…
# 4 nnSVG_7v13.13 7v13          13 Br2720_ant     926     509 FALSE  581.  11.8     595. TRUE  nnSVG_k16-7v13-Br2720_ant.…

log_t |> filter(Sample == "Br2720_ant")

#### load data ####
nnSVG_output <- lapply(output_fn, function(x) get(load(x)))

nnSVG_all <- do.call("rbind", nnSVG_output) |> 
  group_by(domains) |> 
  mutate(FDR = p.adjust(pval, method = "BH"))

dim(nnSVG_all)
# [1] 446052     24

#### Br6522_ant quick export for check ####

nnSVG_quick_check <- nnSVG_all |>
  # filter(Sample == "Br6522_ant", domains == "5v9") |> 
  filter(Sample == "Br6522_ant", rank <= 100) |>
  select(Sample, domains, gene_id, gene_name, rank, pval, FDR) |>
  arrange(rank) |>
  arrange(domains)

nnSVG_quick_check |> filter(rank <= 3) |> print(n=25)


write_csv(nnSVG_quick_check, file = here(data_dir, "nnSVG_quick_check_Br6522_ant-5v9.csv"))


#### Explore Data ####
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




