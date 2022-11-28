
library("tidyverse")
library("here")
library("sessioninfo")
library("patchwork")

## setup output dirs
data_dir <- here("processed-data", "rdata", "spe", "13_nnSVG", "02_compile_nnSVG_output")
if (!dir.exists(data_dir)) dir.create(data_dir)

plot_dir <- here("plots", "13_nnSVG", "02_compile_nnSVG_output")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

#### Progress ####
output_dir <- here("processed-data", "rdata", "spe", "13_nnSVG", "01_nnSVG_pairwise")
output_fn <- list.files(output_dir, pattern = "nnSVG_k16", full.names = TRUE)
names(output_fn) <- gsub("nnSVG_k16-|.RData", "", output_fn)

length(output_fn)
# [1] 172

output_t <- tibble(file = basename(output_fn)) |>
    mutate(label = gsub("nnSVG_k16-|.RData", "", file)) |>
    separate(label, into = c("domains", "Sample"), sep = "-")

Sys.time()
(n_done <- output_t |> count(domains))
# [1] "2022-10-28 11:26:42 EDT"
# domains     n
# <chr>   <int>
# 1 12v13      28
# 2 12v16      26
# 3 4v16       30
# 4 5v9        30
# 5 7v12       28
# 6 7v13       30


#### Get info from Log files ####
## function to overcome incomplete logs
get_log_info <- function(logs, pattern, add = 0) {
    info <- map(logs, ~ (.x[grep(pattern, .x) + add]))
    info <- modify_if(info, ~ length(.) == 0, ~NA_character_)
    return(info)
}

## Read log files
log_dir <- here("code", "analysis", "13_nnSVG", "logs")
log_fn <- list.files(log_dir, pattern = "01_nnSVG_pairwise_loop_", full.names = TRUE)
names(log_fn) <- gsub("01_nnSVG_pairwise_loop_", "nnSVG_", gsub(".txt", "", basename(log_fn)))

length(log_fn)
# [1] 180

logs <- map(log_fn, readLines)

## Job completed
jobs_done <- map_lgl(logs, ~ any(grepl("Job ends", .x)))

table(jobs_done)
# FALSE  TRUE
# 2   178

n_spots <- get_log_info(logs, "ncol")
error_spot <- grepl("Error", n_spots)
n_spots <- ifelse(error_spot, NA, parse_number(unlist(n_spots)))

error <- grepl("Error", logs)
table(grepl("Error: BiocParallel errors", logs))

nnSVG_log_details <- tibble(
    log = names(log_fn),
    Sample = gsub(".*(Br\\d+_.*),.*", "\\1", get_log_info(logs, "Running sample: ")),
    n_spots = n_spots,
    n_genes = parse_number(unlist(get_log_info(logs, "nrow"))),
    error_spot = error_spot,
    error_any = error,
    runtime = unlist(get_log_info(logs, "Reproducibility information:", 3)),
    done = jobs_done
) |>
    separate(runtime, into = c("user", "system", "elapsed"), sep = " +", convert = TRUE, extra = "drop") |>
    separate(log, into = c(NA, "domains", "sample_i"), convert = TRUE, remove = FALSE) |>
    left_join(output_t)

Sys.time()
nnSVG_log_details |> count(domains, done)

# [1] "2022-10-27 16:13:22 EDT"
# domains done      n
# <chr>   <lgl> <int>
# 1 12v13   TRUE     30
# 2 12v16   FALSE     1
# 3 12v16   TRUE     29
# 4 4v16    TRUE     30
# 5 5v9     TRUE     30
# 6 7v12    FALSE     1
# 7 7v12    TRUE     29
# 8 7v13    TRUE     30

log_recap <- nnSVG_log_details |>
    group_by(domains) |>
    summarize(
        logs = n(),
        error_spot = sum(error_spot),
        error_other = sum(error_any) - error_spot,
        output_file = sum(!is.na(file)),
        done = sum(done)
    )

# domains  logs error_spot error_other output_file  done
# <chr>   <int>      <int>       <int>       <int> <int>
# 1 12v13      30          1           2          27    30
# 2 12v16      30          2           1          26    29
# 3 4v16       30          0           0          30    30
# 4 5v9        30          0           0          30    30
# 5 7v12       30          0           1          28    29
# 6 7v13       30          0           0          30    30

#### Identify uncompleted runs ####

## What didn't finish? - ran for 2+ days
nnSVG_log_details |> filter(!done)
# log            domains sample_i Sample     n_spots n_genes error_spot error_any  user system elapsed done  file
# <chr>          <chr>      <int> <chr>        <dbl>   <dbl> <lgl>      <lgl>     <dbl>  <dbl>   <dbl> <lgl> <chr>
# 1 nnSVG_12v16.23 12v16         23 Br6522_mid     307    4039 FALSE      FALSE        NA     NA      NA FALSE NA
# 2 nnSVG_7v12.4   7v12           4 Br3942_ant     645    3980 FALSE      FALSE        NA     NA      NA FALSE NA

## What errored out? not from spot error

nnSVG_log_details |> filter(!error_spot, error_any)
# log            domains sample_i Sample     n_spots n_genes error_spot error_any  user system elapsed done  file
# <chr>          <chr>      <int> <chr>        <dbl>   <dbl> <lgl>      <lgl>     <dbl>  <dbl>   <dbl> <lgl> <chr>
# 1 nnSVG_12v13.10 12v13         10 Br8492_ant      81    1266 FALSE      TRUE         NA     NA      NA TRUE  NA
# 2 nnSVG_12v16.28 12v16         28 Br8667_ant      74    3892 FALSE      TRUE         NA     NA      NA TRUE  NA
# 3 nnSVG_7v12.19  7v12          19 Br6471_ant     376    3884 FALSE      TRUE         NA     NA      NA TRUE  NA

## redo tab
# redo_tab <- nnSVG_log_details |> filter((!error_spot & error_any)|!done)

## Use SGE job tools to re-run
# walk2(redo_tab$domains, redo_tab$sample_i, function(d, i){
#   file.remove((here(log_dir, paste0("01_nnSVG_pairwise_loop_",d,".",i,".txt"))))
#   sgejobs::array_submit(paste0(".01_nnSVG_pairwise_loop_",d,".sh"), i)
# })

## save
save(nnSVG_log_details, file = here(data_dir, "nnSVG_log_details.Rdata"))

#### Plot runtime and other input details ####
runtime_long <- nnSVG_log_details |>
    filter(done) |>
    mutate(prod = n_spots * n_genes) |>
    select(domains, Sample, n_spots, n_genes, elapsed, user, prod) |>
    pivot_longer(!c(domains, Sample, elapsed, user), names_to = "metric", values_to = "n")

runtime_plot <- runtime_long |>
    ggplot(aes(x = n, y = elapsed / 60, color = domains)) +
    geom_point() +
    facet_wrap(~metric, scales = "free_x") +
    labs(y = "elapsed (min)")

ggsave(runtime_plot, filename = here(plot_dir, "nnSVG_runtime_scatter.png"), width = 10)

## cells vs. genes
spots_v_gene_scatter <- nnSVG_log_details |>
    ggplot(aes(x = n_spots, y = n_genes, color = domains, shape = error_any)) +
    geom_point() +
    geom_vline(xintercept = 65, color = "red", linetype = "dashed")

ggsave(spots_v_gene_scatter, filename = here(plot_dir, "nnSVG_spots_v_gene_scatter.png"))

spots_box <- nnSVG_log_details |>
    ggplot(aes(x = domains, y = n_spots, fill = domains)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.25) +
    geom_jitter(aes(color = domains)) +
    geom_hline(yintercept = 65, color = "red", linetype = "dashed")

ggsave(spots_box, filename = here(plot_dir, "nnSVG_n_spots_boxplot.png"))

genes_box <- nnSVG_log_details |>
    ggplot(aes(x = domains, y = n_genes, fill = domains)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.25) +
    geom_jitter(aes(color = domains))

ggsave(genes_box, filename = here(plot_dir, "nnSVG_n_genes_boxplot.png"))

ggsave(
    spots_v_gene_scatter + theme(legend.position = "None") +
        genes_box /
            (spots_box + coord_flip() + theme(legend.position = "None")),
    filename = here(plot_dir, "nnSVG_input_summary.png"), width = 10
)

## low genes all from one sample
nnSVG_log_details |> filter(n_genes < 1000)
# log           domains sample_i Sample     n_spots n_genes error  user system elapsed done  file
# <chr>         <chr>      <int> <chr>        <dbl>   <dbl> <lgl> <dbl>  <dbl>   <dbl> <lgl> <chr>
# 1 nnSVG_4v16.13 4v16          13 Br2720_ant    1390     521 FALSE  574.   5.59    582. TRUE  nnSVG_k16-4v16-Br2720_ant.…
# 2 nnSVG_5v9.13  5v9           13 Br2720_ant    1401     710 FALSE 1127.  12.6    1143. TRUE  nnSVG_k16-5v9-Br2720_ant.R…
# 3 nnSVG_7v12.13 7v12          13 Br2720_ant     873     459 FALSE 1430.  74.4    1507. TRUE  nnSVG_k16-7v12-Br2720_ant.…
# 4 nnSVG_7v13.13 7v13          13 Br2720_ant     926     509 FALSE  581.  11.8     595. TRUE  nnSVG_k16-7v13-Br2720_ant.…

nnSVG_log_details |> filter(Sample == "Br2720_ant")

#### load nnSVG output data ####
nnSVG_output <- lapply(output_fn, function(x) get(load(x)))

nnSVG_all <- do.call("rbind", nnSVG_output) |>
    group_by(domains) |>
    mutate(FDR = p.adjust(pval, method = "BH"))

dim(nnSVG_all)
# [1] 450635     24

save(nnSVG_all, file = here(data_dir, "nnSVG_all.Rdata"))

## Check output again...
nnSVG_all |>
    group_by(domains) |>
    summarize(output_file = length(unique(Sample))) |>
    left_join(log_recap) |>
    mutate(not_done = logs - done) |>
    select(-logs, -done)

# domains output_file error_spot error_other not_done
# <chr>         <int>      <int>       <int>    <int>
# 1 12v13            27          1           2        0
# 2 12v16            26          2           1        1
# 3 4v16             30          0           0        0
# 4 5v9              30          0           0        0
# 5 7v12             28          0           1        1
# 6 7v13             30          0           0        0


#### Br6522_ant quick export for check ####

nnSVG_quick_check <- nnSVG_all |>
    filter(Sample == "Br6522_ant", domains == "5v9") |>
    filter(Sample == "Br6522_ant", rank <= 100) |>
    select(Sample, domains, gene_id, gene_name, rank, pval, FDR) |>
    arrange(rank) |>
    arrange(domains)

nnSVG_quick_check |>
    filter(rank <= 3) |>
    print(n = 25)


write_csv(nnSVG_quick_check, file = here(data_dir, "nnSVG_quick_check_Br6522_ant-5v9.csv"))

# sgejobs::job_single('02_compile_nnSVG_output', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript 02_compile_nnSVG_output.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
