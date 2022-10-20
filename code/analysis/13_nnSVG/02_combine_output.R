
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
# [1] "2022-10-20 16:05:36 EDT"
# domains     n
# <chr>   <int>
# 1 12v13      24
# 2 12v16      11
# 3 4v16       30
# 4 5v9        29
# 5 7v12       26
# 6 7v13       21

domains <- readLines("nnSVG_domains.txt")


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

#### load data ####
nnSVG_output <- lapply(output_fn, function(x) get(load(x)))

nnSVG_all <- do.call("rbind", nnSVG_output) |> 
  group_by(domains) |> 
  mutate(FDR = p.adjust(pval, method = "BH"))

dim(nnSVG_all)
# [1] 364546     23


nnSVG_all |> arrange(rank) |> head() 

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

nnSVG_all_summary |> arrange(min_rank)

nnSVG_all_summary |> filter(n == 30)
nnSVG_all_summary |> filter(n_signif != 0) |> arrange(-n_signif)

nnSVG_all_summary |> filter(max_FDR !=1) |> ungroup()  |> count(n ==1)

nnSVG_all_summary |> filter(n > 20) |> arrange(max_rank) |> slice(1)

n_gene_distribution <- nnSVG_all_summary |>
  ggplot(aes(n, color = domains)) +
  geom_density()

ggsave(n_gene_distribution, filename = here(plot_dir, "n_gene_distribution.png"), width = 10)


