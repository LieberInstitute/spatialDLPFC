
library("tidyverse")
library("here")
library("sessioninfo")

#### Progress ####
data_dir <- here("processed-data", "rdata", "spe", "13_nnSVG", "01_nnSVG_pairwise")
output_fn <- list.files(data_dir, pattern = "nnSVG_k16", full.names = TRUE)
names(output_fn) <- gsub("nnSVG_k16-|.RData","", output_fn)

length(output_fn)
# [1] 141

output_t <- tibble(file = basename(output_fn)) |>
            mutate(label = gsub("nnSVG_k16-|.RData","", file)) |>
  separate(label, into = c("domains","Sample"), sep = "-")

Sys.time()
(n_done <- output_t |> count(domains))
[1] "2022-10-20 17:19:06 EDT"
# # A tibble: 6 × 2
# domains     n
# <chr>   <int>
# 1 12v13      26
# 2 12v16      20
# 3 4v16       30
# 4 5v9        30
# 5 7v12       27
# 6 7v13       26


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
# [1] 150

logs <- map(log_fn, readLines)

## Job completed
map_lgl(logs, ~any(grepl("Job ends", .x)))

n_spots <- get_log_info(logs, "ncol")
error <- grepl("Error",ncol)
n_spots <- ifelse(error, NA, parse_number(unlist(ncol)))

log_t <- tibble(log = names(log_fn),
                Sample = gsub(".*(Br\\d+_.*),.*","\\1",get_log_info(logs, "Running sample: ")),
                n_spots =  get_log_info(logs, "ncol"),
                n_genes = parse_number(unlist(get_log_info(logs, "nrow"))),
                error = error,
                runtime = unlist(get_log_info(logs, "Reproducibility information:", 3))) |> 
  separate(runtime, into = c("user", "system", "elapsed"), sep = " +", convert = TRUE, extra = "drop") |>
  separate(log, into = c(NA, "domains", "sample_i"), convert = TRUE, remove = FALSE) |>
  mutate(done = !is.na(user),
         error = grepl("Error",ncol),
         n_spots = ifelse(error, NA, parse_number(unlist(ncol)))) |>
  left_join(output_t)

Sys.time()
log_t |> count(domains, done)

# [1] "2022-10-20 15:16:27 EDT"
# domains done      n
# <chr>   <lgl> <int>
# 1 12v13   FALSE    13
# 2 12v13   TRUE     16
# 3 12v16   FALSE     3
# 4 12v16   TRUE      9
# 5 4v16    TRUE     30
# 6 5v9     FALSE     1
# 7 5v9     TRUE     29
# 8 7v12    FALSE    11
# 9 7v12    TRUE     18
# 10 7v13    FALSE     7
# 11 7v13    TRUE     15

log_t |> filter(error) |> count(domains)

log_t |> filter(done, !error, is.na(file))
log_t |> filter(done, !error) |> count(domains, output = !is.na(file))


## Checkout runtime
plot_dir <- here("plots", "13_nnSVG")

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
# [1] 414846     24

nnSVG_all  |> filter(FDR < 0.05) |> count()

nnSVG_all_summary <- nnSVG_all |>
  group_by(domains, gene_id) |> 
  summarize(n = n(), 
            n_signif = sum(FDR < 0.05),
            max_FDR = max(FDR),
            min_FDR = min(FDR),
            min_rank = min(rank),
            max_rank = max(rank),
            mean_rank = mean(rank), 
            median_rank = median(rank))

nnSVG_all_summary |> filter(n == 30) |> count()
# domains     n
# <chr>   <int>
# 1 4v16      465
# 2 5v9       568

nnSVG_all_summary |> filter(n_signif != 0) |> arrange(-n_signif)

nnSVG_all_summary |> filter(max_FDR !=1) |> ungroup()  |> count(n ==1)

nnSVG_all_summary |> filter(max_FDR !=1, n !=1)

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
