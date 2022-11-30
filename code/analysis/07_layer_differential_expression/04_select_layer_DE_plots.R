library("SpatialExperiment")
library("scater")
library("spatialLIBD")
library("here")
library("sessioninfo")


plot_dir <-
    here(
        "plots",
        "07_layer_differential_expression",
        "09_select_layer_DE_plots"
    )
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
}

## Read in the custom plotting function
source(
    here(
        "code",
        "analysis",
        "07_layer_differential_expression",
        "custom_plotExpression.R"
    ),
    echo = TRUE,
    max.deparse.length = 500
)

#### K9 violin plots ####
spe_k9 <-
    readRDS(
        file = here(
            "processed-data",
            "rdata",
            "spe",
            "07_layer_differential_expression",
            "sce_pseudo_BayesSpace_k09.rds"
        )
    )
rownames(spe_k9) <- rowData(spe_k9)$gene_name

## establish color pallet
k9_colors <- Polychrome::palette36.colors(9)
names(k9_colors) <- c(1:9)

k9_genes <-
    c("CLDN5", "TAGLN", "MYL9", "ACTA2", "SLC2A1", "HBA1", "EPAS1")

# sce, genes, assay = "logcounts", cat, highlight = "highlight", fill_colors = NULL, title = NULL
k9_plot <-
    custom_plotExpression(
        spe_k9,
        genes = k9_genes,
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k9_colors
    )
ggsave(
    k9_plot,
    filename = here(plot_dir, "k9_expression_meninges.png"),
    height = 10
)


k9_plot_CLDN5 <-
    custom_plotExpression(
        spe_k9,
        genes = c("CLDN5"),
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k9_colors,
        highlight_sample = "Br6522_ant"
    ) +
    labs(x = "Pseudobulk k9 Domains", title = "1 > others p=2.04e-75")
ggsave(
    k9_plot_CLDN5,
    filename = here(plot_dir, "k9_expression_meninges_CLDN5.png"),
    height = 5
)


#### K16 violin plots ####
spe_k16 <-
    readRDS(
        file = here(
            "processed-data",
            "rdata",
            "spe",
            "07_layer_differential_expression",
            "sce_pseudo_BayesSpace_k16.rds"
        )
    )
rownames(spe_k16) <- rowData(spe_k16)$gene_name

## establish color pallet
k16_colors <- Polychrome::palette36.colors(16)
names(k16_colors) <- c(1:16)

## 1a vs. 1b
k16_genes_1ab <- c("SPARC", "MSX1", "RELN", "APOE")
k16_genes %in% rownames(spe_k16)

spe_k16$highlight <- spe_k16$BayesSpace %in% c(1, 2)

k16_1ab_plot <-
    custom_plotExpression(
        spe_k16,
        genes = k16_genes_1ab,
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k16_colors
    )
ggsave(
    k16_1ab_plot,
    filename = here(plot_dir, "k16_expression_1a-1b.png"),
    height = 10
)

spe_k16$highlight <- FALSE

k16_plot_SPARC <-
    custom_plotExpression(
        spe_k16[, spe_k16$BayesSpace %in% c(2, 14)],
        genes = c("SPARC"),
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k16_colors,
        highlight_sample = "Br6522_ant"
    ) +
    labs(x = "Pseudobulk k16 Domains", title = "14 > 2 p=1.52e-13")
ggsave(
    k16_plot_SPARC,
    filename = here(plot_dir, "k16_expression_1a-1b_SPARC.png"),
    width = 6,
    height = 3
)

k16_plot_HTRA1 <-
    custom_plotExpression(
        spe_k16[, spe_k16$BayesSpace %in% c(2, 14)],
        genes = c("HTRA1"),
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k16_colors,
        highlight_sample = "Br6522_ant"
    ) +
    labs(x = "Pseudobulk k16 Domains", title = "2>14 p=3.34e-07")
ggsave(
    k16_plot_HTRA1,
    filename = here(plot_dir, "k16_expression_1a-1b_HTRA1.png"),
    width = 6,
    height = 3
)

## plot both
k16_plot_SPARC_HTRA1 <-
    custom_plotExpression(
        spe_k16[, spe_k16$BayesSpace %in% c(2, 14)],
        genes = c("SPARC", "HTRA1"),
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k16_colors
    ) +
    labs(x = "Pseudobulk k16 Domains")

ggsave(
    k16_plot_SPARC_HTRA1,
    filename = here(plot_dir, "k16_expression_1a-1b_SPARC_HTRA1.png"),
    height = 3,
    width = 11
)


## 6a
k16_6a_plot <- custom_plotExpression(
    spe_k16,
    genes = c("SMIM32", "DACH1", "KIF1A", "GALNT14"),
    assay = "logcounts",
    cat = "BayesSpace",
    fill_colors = k16_colors
)

ggsave(k16_6a_plot, filename = here(plot_dir, "k16_expression_6a.png"))

## 6b
k16_6b_plot <- custom_plotExpression(
    spe_k16,
    genes = c("KRT17", "DIRAS2", "SEMA3E"),
    assay = "logcounts",
    cat = "BayesSpace",
    fill_colors = k16_colors
)

ggsave(k16_6b_plot, filename = here(plot_dir, "k16_expression_6b.png"))


##### Spatial Dot plots ####
# Load full data
load(
    here(
        "processed-data",
        "rdata",
        "spe",
        "01_build_spe",
        "spe_filtered_final_with_clusters.Rdata"
    )
)

## prep objet for plotting
rowData(spe)$gene_search <- rowData(spe)$gene_name
spe$bayesSpace_harmony_9 <- as.factor(spe$bayesSpace_harmony_9)
spe$bayesSpace_harmony_16 <- as.factor(spe$bayesSpace_harmony_16)

## K9 CLDN5
pdf(here(plot_dir, "vis_gene_CLDN5-Br6522_ant.pdf"))
vis_gene(
    spe = spe,
    sampleid = "Br6522_ant",
    geneid = "CLDN5"
)
dev.off()

png(here(plot_dir, "vis_gene_CLDN5-Br6522_ant.png"))
vis_gene(
    spe = spe,
    sampleid = "Br6522_ant",
    geneid = "CLDN5"
)
dev.off()


## K16 SPARC
pdf(here(plot_dir, "vis_gene_SPARC-Br6522_ant.pdf"))
vis_gene(
    spe = spe,
    sampleid = "Br6522_ant",
    geneid = "SPARC"
)
dev.off()

## K16 SPARC
pdf(here(plot_dir, "vis_gene_HTRA1-Br6522_ant.pdf"))
vis_gene(
    spe = spe,
    sampleid = "Br6522_ant",
    geneid = "HTRA1"
)
dev.off()

## filter to just 2 & 14
pdf(here(plot_dir, "vis_gene_SPARC-Br6522_ant_filter.pdf"))
vis_gene(
    spe = spe[, spe$bayesSpace_harmony_16 %in% c(2, 14)],
    sampleid = "Br6522_ant",
    geneid = "SPARC"
)
dev.off()

## K16 SPARC
pdf(here(plot_dir, "vis_gene_HTRA1-Br6522_ant_filter.pdf"))
vis_gene(
    spe = spe[, spe$bayesSpace_harmony_16 %in% c(2, 14)],
    sampleid = "Br6522_ant",
    geneid = "HTRA1"
)
dev.off()

## Visulize clusters 2 & 14
pdf(here(plot_dir, "vis_clus_d16.pdf"))
vis_clus(
    spe = spe,
    clustervar = "bayesSpace_harmony_16",
    colors = k16_colors,
    sampleid = "Br6522_ant"
)
dev.off()

pdf(here(plot_dir, "vis_clust_d16_2-14.pdf"))
vis_clus(
    # spe = spe,
    spe = spe[, spe$bayesSpace_harmony_16 %in% c(2, 14)],
    clustervar = "bayesSpace_harmony_16",
    colors = k16_colors,
    sampleid = "Br6522_ant"
)
dev.off()
