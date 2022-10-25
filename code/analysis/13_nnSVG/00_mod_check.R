
library("SpatialExperiment")
library("tidyverse")
library("here")
library("sessioninfo")

load(here("processed-data", "rdata","spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

k16_domain_pairs <- list(`5v9` = c(5,9),
                         `4v16` = c(4,16),
                         `7v13` = c(7,13),
                         `12v13` = c(12,13),
                         `7v12` = c(7,12),
                         `12v16` = c(12,16))


k16_domains <- unique(unlist(k16_domain_pairs))

pd <- as.data.frame(colData(spe))

mod_tab <- pd |> 
  filter(bayesSpace_harmony_16 %in% k16_domains) |>
  count(sample_id, bayesSpace_harmony_16, bayesSpace_harmony_9)

mod_max <- function(pd = pd, domains, domain_names){
  domain_count <- pd |>
    group_by(sample_id) |>
    filter(bayesSpace_harmony_16 %in% domains)|>
    count(sample_id, bayesSpace_harmony_9) |>
    summarize(total_n = sum(n),
              max_n = max(n),
              max_k9 = bayesSpace_harmony_9[which.max(n)],
              other_n = total_n - max_n) |>
    mutate(domains = domain_names)
  
  return(domain_count)
}

mod_max_df <- map2(k16_domain_pairs, names(k16_domain_pairs), ~mod_max(pd, .x, .y))
mod_max_df <- do.call("rbind", mod_max_df)

mod_max_df |> group_by(domains) |> count(max_k9) |> mutate(all_same_max = n == 30)
# domains max_k9     n all_same_max
# <chr>    <int> <int> <lgl>       
# 1 d12v13       7    11 FALSE       
# 2 d12v13       9    19 FALSE       
# 3 d12v16       4    23 FALSE       
# 4 d12v16       7     7 FALSE       
# 5 d4v16        4    30 TRUE        
# 6 d5v9         8    30 TRUE        
# 7 d7v12        7    30 TRUE        
# 8 d7v13        7    29 FALSE       
# 9 d7v13        9     1 FALSE  

mod_max_df |> group_by(domains) |> summarize(min_other = min(other_n))
# domains min_other
# <chr>       <int>
# 1 d12v13         21
# 2 d12v16         16
# 3 d4v16          13
# 4 d5v9            6
# 5 d7v12           0
# 6 d7v13          54

mod_max_df |> arrange(other_n)
# sample_id   total_n max_n max_k9 other_n domains
# * <chr>         <int> <int>  <int>   <int> <chr>  
# 1 Br6471_ant      432   432      7       0 d7v12  
# 2 Br2720_mid      267   262      7       5 d7v12  
# 3 Br2720_mid      169   163      8       6 d5v9   
# 4 Br6522_mid      350   343      7       7 d7v12  
# 5 Br6471_ant      376   366      8      10 d5v9   
# 6 Br6522_ant      595   585      7      10 d7v12  
# 7 Br2720_post     506   494      8      12 d5v9   
# 8 Br6471_mid      497   484      4      13 d4v16  
# 9 Br8325_mid      600   587      4      13 d4v16  
# 10 Br6522_mid      675   661      4      14 d4v16 


## Compare output from nnSVG with and w/o model
load(here("processed-data", "rdata", "spe", "13_nnSVG", "01_nnSVG_pairwise","nnSVG_k16-5v9-Br2743_ant.RData"), verbose = TRUE)
nnSVG_mod_data <- nnSVG_data
dim(nnSVG_mod_data)
# [1] 2153   23
load(here("processed-data", "rdata", "spe", "13_nnSVG", "01_nnSVG_pairwise_modTest","nnSVG_k16-5v9-Br2743_ant.RData"), verbose = TRUE)
dim(nnSVG_data)
# [1] 2153   23

## all data is different
all_equal(nnSVG_mod_data, nnSVG_data)

colnames(nnSVG_mod_data)
head(nnSVG_mod_data)

nnSVG_mod_data <- nnSVG_mod_data|>
  select(gene_id, 8:21) |>
  rename_with(.cols = !gene_id, function(x){paste0(x,"_mod")}) 

nnSVG_all_data <- left_join(nnSVG_data, nnSVG_mod_data)

#### Compare Results ####
plot_dir <- here("plots", "13_nnSVG", "mod_check")

rank_scatter <- ggplot(nnSVG_all_data, aes(x = rank_mod, y = rank, color = padj !=1)) +
  geom_point() +
  geom_abline(color = "red", linetype = "dashed")

ggsave(rank_scatter, filename = here(plot_dir, "mod_rank_scatter.png"))

padj_scatter <- ggplot(nnSVG_all_data, aes(x = -log10(pval_mod), y = -log10(pval))) +
  geom_point() +
  geom_abline(color = "red", linetype = "dashed")

ggsave(padj_scatter, filename = here(plot_dir, "mod_padj_scatter.png"))

nnSVG_all_data |> select(gene_id, padj, padj_mod) |> filter(padj != 1)



