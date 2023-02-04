library("SingleCellExperiment")
library("spatialLIBD")
library("tidyverse")
library("scMerge")
library("here")
library("sessioninfo")

## Plot setup
plot_dir <- here("plots", "14_spatial_registration_PEC", "04_PEC_check_expression")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "rdata", "spe", "14_spatial_registration_PEC", "04_PEC_check_expression")
if (!dir.exists(data_dir)) dir.create(data_dir)

## Load data
pseudobulk_fn <- list.files(here("processed-data", "rdata", "spe", "14_spatial_registration_PEC"), pattern = "pseudobulk", full.names = TRUE)
names(pseudobulk_fn) <- gsub("pseudobulk_|.rds", "", basename(pseudobulk_fn))

all_pseudobulk <- map2(pseudobulk_fn, names(pseudobulk_fn), function(fn, name) {
    message(Sys.time(), " - loading ", name)
    pb <- readRDS(fn)
    pb$Dataset <- name
    return(pb)
})

map(all_pseudobulk, dim)
map_int(all_pseudobulk, nrow)
# CMC DevBrain-snRNAseq            IsoHuB         SZBDMulti          UCLA-ASD       Urban-DLPFC
# 29460             23405             27461             30912             30440             24273

all_pb_sce <- sce_cbind(all_pseudobulk, colData_names = c("Dataset", "subclass", "ncells", "sampleID"))

dim(all_pb_sce)
# [1] 22138 13297

## blanks out rowData...replace from pb data
rowData(all_pb_sce) <- rowData(all_pseudobulk$CMC)[rownames(all_pb_sce), c("featureid", "gene_name")]
rownames(all_pb_sce) <- rowData(all_pb_sce)$gene_name


## Check out colData
table(all_pb_sce$Dataset, all_pb_sce$subclass)
# Astro Chandelier Endo L2/3 IT L4 IT L5 ET L5 IT L5/6 NP L6 CT L6 IT L6 IT Car3 L6b Lamp5 Lamp5 Lhx6
# CMC                  99         77   52     101    99    35   100      87    95    98         76  82    93         95
# DevBrain-snRNAseq    16         16   16      16    16    13    16      16    16    16         16  14    15         16
# IsoHuB                5          5    5       5     5     1     5       5     4     5          4   3     5          5
# SZBDMulti           487        261   15     567   537   236   556     343   418   535        395 466   469        384
# UCLA-ASD             63         54   50      63    61    39    63      50    55    62         50  52    57         53
# Urban-DLPFC          21         19   21      21    21     7    21      21    21    21         18  21    20         20
#
# Micro/PVM OPC Oligo Pax6 Pvalb Sncg Sst Sst Chodl VLMC Vip
# CMC                      85  97    99   75    99   90  97         5   95  98
# DevBrain-snRNAseq        14  16    16   12    16   14  16         0   16  16
# IsoHuB                    5   5     5    4     5    5   5         0    5   5
# SZBDMulti               285 445   529  135   532  365 505         2   64 527
# UCLA-ASD                 60  63    63   41    62   57  62         2   58  62
# Urban-DLPFC              20  21    21   12    21   18  21         0   21  20

pb_pd <- as.data.frame(colData(all_pb_sce))

## Cell type breakdown of the data
ct_summary <- pb_pd |>
    group_by(Dataset, subclass) |>
    summarize(
        total_cells = sum(ncells),
        n_groups = n()
    ) |>
    group_by(Dataset) |>
    mutate(prop = total_cells / sum(total_cells))

ct_barplot <- ggplot(data = ct_summary, aes(x = Dataset, y = total_cells, fill = subclass)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = ifelse(total_cells > 9000, total_cells, "")),
        size = 2.7,
        position = position_stack(vjust = 0.5)
    ) +
    theme(legend.position = "bottom")

ggsave(ct_barplot, filename = here(plot_dir, "ct_barplot.png"))


ct_prop_barplot <- ggplot(data = ct_summary, aes(x = Dataset, y = prop, fill = subclass)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")),
        size = 2.7,
        position = position_stack(vjust = 0.5)
    ) +
    theme(legend.position = "bottom")

ggsave(ct_prop_barplot, filename = here(plot_dir, "ct_prop_barplot.png"))

ct_group_dotplot <- ggplot(ct_summary, aes(x = Dataset, y = subclass)) +
    geom_point(aes(size = n_groups, color = subclass)) +
    geom_text(aes(label = n_groups)) +
    theme(legend.position = "None")

ggsave(ct_group_dotplot, filename = here(plot_dir, "ct_group_dotplot.png"))


## Genes to explore
SCZ_interactions <- read.csv(here(data_dir, "SCZ_top_interactions.csv"))
genes_of_interest <- unique(unlist(SCZ_interactions))

genes_of_interest[!genes_of_interest %in% rownames(all_pb_sce)]
# [1] "FSHB" missing
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(all_pb_sce)]
length(genes_of_interest)

#### PLOT Expression ####
cat_df <- as.data.frame(colData(all_pb_sce))[, c("sampleID", "subclass"), drop = FALSE]

expression_long <- reshape2::melt(as.matrix(logcounts(all_pb_sce)[genes_of_interest, , drop = FALSE]))
cat <- cat_df[expression_long$Var2, ]
colnames(cat) <- c("sample_id", "cat")
expression_long <- cbind(expression_long, cat)

pdf(here(plot_dir, "PEC_SCZ_interactions.pdf"))
for (i in 1:nrow(SCZ_interactions)) {
    interaction <- SCZ_interactions[i, ]
    interaction_title <- paste(interaction[[1]], "->", interaction[[2]])
    message(i, " ", interaction_title)

    expression_temp <- expression_long |>
        filter(Var1 %in% interaction) |>
        mutate(Var1 = factor(Var1, levels = interaction))

    expression_violin <- ggplot(data = expression_temp, aes(x = cat, y = value, fill = cat)) +
        geom_violin(scale = "width") +
        # geom_jitter(aes(color = cat), size = .3) +
        facet_wrap(~Var1, ncol = 1, scales = "free_y") +
        labs(
            title = interaction_title,
            y = "Expression (logcounts)"
        ) +
        theme_bw() +
        theme(
            legend.position = "None",
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(face = "italic"),
            text = element_text(size = 15)
        ) +
        stat_summary(
            fun = median,
            # fun.min = median,
            # fun.max = median,
            geom = "crossbar",
            width = 0.3
        )
    print(expression_violin)
}
dev.off()
