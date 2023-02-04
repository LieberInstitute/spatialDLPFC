library("scater")
library("purrr")
library("SpatialExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")


plot_dir <-
    here(
        "plots",
        "07_layer_differential_expression",
        "04_select_layer_DE_plots"
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

## Load in pseudo-bulked objects
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

## Load the modeling results to extract the p-values
k_list <- c(9, 16)
names(k_list) <- c("k09", "k16")
model_stats <- map(k_list, function(k) {
    get(load(
        file = here(
            "processed-data",
            "rdata",
            "spe",
            "07_layer_differential_expression",
            paste0(
                "modeling_results_BayesSpace_k",
                sprintf("%02d", k),
                ".Rdata"
            )
        ),
        verbose = TRUE
    ))
})

spe_k9$spatialLIBD <- spe_k9$BayesSpace
spe_k16$spatialLIBD <- spe_k16$BayesSpace

# ## for gene matching
# rownames(sce_k9) <- rowData(sce_k9)$gene_id
# rownames(sce_k16) <-  rowData(sce_k16)$gene_id

sig_genes <-
    map2(
        model_stats,
        list(spe_k9, spe_k16),
        # list(spe_k9),
        ~ sig_genes_extract_all(
            n = nrow(.y),
            modeling_results = .x,
            sce_layer = .y
        )
    )


## For custom_plotExpression()
rownames(spe_k9) <- rowData(spe_k9)$gene_name
rownames(spe_k16) <- rowData(spe_k16)$gene_name


#### K9 violin plots ####

## establish color pallet
(k9_colors <-
    spe_k9$BayesSpace_colors[!duplicated(spe_k9$BayesSpace_colors)])
# Sp09D01   Sp09D02   Sp09D03   Sp09D04   Sp09D05   Sp09D06   Sp09D07   Sp09D08   Sp09D09
# "#5A5156" "#E4E1E3" "#F6222E" "#FE00FA" "#16FF32" "#3283FE" "#FEAF16" "#B00068" "#1CFFCE"

k9_genes <-
    c("CLDN5", "TAGLN", "MYL9", "ACTA2", "SLC2A1", "HBA1", "EPAS1")

subset(
    sig_genes$k09,
    model_type == "enrichment" &
        test == "Sp09D01" &
        gene %in% k9_genes
)[, c("top", "gene", "stat", "pval", "fdr", "ensembl")]
# DataFrame with 7 rows and 6 columns
#         top        gene      stat        pval         fdr         ensembl
#   <integer> <character> <numeric>   <numeric>   <numeric>     <character>
# 1         1       CLDN5   26.1024 2.04113e-75 2.49528e-71 ENSG00000184113
# 2         2       TAGLN   25.3402 5.75174e-73 3.51575e-69 ENSG00000149591
# 3         3        MYL9   22.9466 4.47117e-65 1.82200e-61 ENSG00000101335
# 4         4       ACTA2   21.6077 1.53905e-60 4.70372e-57 ENSG00000107796
# 5         6      SLC2A1   20.4544 1.43992e-56 2.93383e-53 ENSG00000117394
# 6         7       EPAS1   20.2733 6.11909e-56 9.36730e-53 ENSG00000116016
# 7        10        HBA1   19.8770 1.46555e-54 1.79164e-51 ENSG00000206172

# sce, genes, assay = "logcounts", cat, highlight = "highlight", fill_colors = NULL, title = NULL
set.seed(20221130) ## to make the jitter reproducible
k9_plot <-
    custom_plotExpression(
        spe_k9,
        genes = k9_genes,
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k9_colors
    )
ggsave(
    k9_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = NULL),
    filename = here(plot_dir, "k9_expression_meninges.png"),
    height = 10,
    width = 10
)

## Locate the pvalue of interest
k9_plot_CLDN5_pval <-
    subset(
        sig_genes$k09,
        model_type == "enrichment" &
            test == "Sp09D01" & gene == "CLDN5"
    )$pval

set.seed(20221130) ## to make the jitter reproducible
k9_plot_CLDN5 <-
    custom_plotExpression(
        spe_k9,
        genes = c("CLDN5"),
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k9_colors
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Pseudobulk k9 Domains", title = paste0("Sp09D01 > others p=", signif(k9_plot_CLDN5_pval, 3)))

ggsave(
    k9_plot_CLDN5,
    filename = here(plot_dir, "k9_expression_meninges_CLDN5.png"),
    height = 4.6
)

ggsave(
    k9_plot_CLDN5,
    filename = here(plot_dir, "k9_expression_meninges_CLDN5.pdf"),
    height = 4.6
)


#### K16 violin plots ####

## establish color pallet
k16_colors <-
    spe_k16$BayesSpace_colors[!duplicated(spe_k16$BayesSpace_colors)]

## 1a vs. 1b
k16_genes_1ab <- c("SPARC", "MSX1", "RELN", "APOE")
stopifnot(all(k16_genes_1ab %in% rownames(spe_k16)))

subset(
    sig_genes$k16,
    model_type == "pairwise" &
        test %in% c("Sp16D02-Sp16D14", "Sp16D14-Sp16D02") &
        gene %in% k16_genes_1ab
)[, c("top", "gene", "stat", "pval", "fdr", "ensembl")]
# DataFrame with 8 rows and 6 columns
#         top        gene      stat        pval         fdr         ensembl
#   <integer> <character> <numeric>   <numeric>   <numeric>     <character>
# 1        89        APOE  3.746661 2.02384e-04 1.79653e-02 ENSG00000130203
# 2      6039        RELN  0.342777 7.31925e-01 8.86246e-01 ENSG00000189056
# 3      9546        MSX1 -3.108805 1.99697e-03 5.92724e-02 ENSG00000163132
# 4      9585       SPARC -7.718414 7.61269e-14 2.43276e-10 ENSG00000113140
# 5         3       SPARC  7.718414 7.61269e-14 2.43276e-10 ENSG00000113140
# 6        42        MSX1  3.108805 1.99697e-03 5.92724e-02 ENSG00000163132
# 7      3549        RELN -0.342777 7.31925e-01 8.86246e-01 ENSG00000189056
# 8      9499        APOE -3.746661 2.02384e-04 1.79653e-02 ENSG00000130203

## Doesn't seem to work with the current version
## It likely did prior to this change in custom_plotExpression()
## https://github.com/LieberInstitute/spatialDLPFC/commit/52118574738d19cc00a3c299a8ed2677d9fb28dd

# spe_k16$highlight <- spe_k16$BayesSpace %in% c("Sp16D01", "Sp16D02")

set.seed(20221130) ## to make the jitter reproducible
k16_1ab_plot <-
    custom_plotExpression(
        spe_k16,
        genes = k16_genes_1ab,
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k16_colors
    )
ggsave(
    k16_1ab_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = NULL),
    filename = here(plot_dir, "k16_expression_1a-1b.png")
)

spe_k16$highlight <- FALSE

k16_plot_SPARC_pval <- subset(
    sig_genes$k16,
    model_type == "pairwise" &
        test %in% c("Sp16D14-Sp16D02") & gene == "SPARC"
)$pval

set.seed(20221130) ## to make the jitter reproducible
k16_plot_SPARC <-
    custom_plotExpression(
        spe_k16[, spe_k16$BayesSpace %in% c("Sp16D02", "Sp16D14")],
        genes = c("SPARC"),
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k16_colors
        # highlight_sample = "Br6522_ant"
    ) +
    # geom_line(aes(group = sample_id))+
    labs(
        x = "Pseudobulk k16 Domains",
        title = "Sp16D02 < Sp16D14",
        subtitle = paste0("p=", signif(k16_plot_SPARC_pval, 3))
    )

ggsave(
    k16_plot_SPARC +
        # theme(axis.text.x = element_text(angle = 90, hjust = 1))
        theme(axis.title.x = element_blank()),
    filename = here(plot_dir, "k16_expression_1a-1b_SPARC.pdf"),
    width = 3,
    height = 4.5
)

k16_plot_HTRA1_pval <- subset(
    sig_genes$k16,
    model_type == "pairwise" &
        test %in% c("Sp16D02-Sp16D14") & gene == "HTRA1"
)$pval

set.seed(20221130) ## to make the jitter reproducible
k16_plot_HTRA1 <-
    custom_plotExpression(
        spe_k16[, spe_k16$BayesSpace %in% c("Sp16D02", "Sp16D14")],
        genes = c("HTRA1"),
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k16_colors
        # highlight_sample = "Br6522_ant"
    ) +
    # geom_line(aes(group = sample_id)) +
    labs(
        x = "Pseudobulk k16 Domains",
        title = "Sp16D02 > Sp16D14",
        subtitle = paste0("p=", signif(k16_plot_HTRA1_pval, 3))
    )
ggsave(
    k16_plot_HTRA1 +
        # + theme(axis.text.x = element_text(angle = 90, hjust = 1)),
        theme(axis.title.x = element_blank(), ),
    filename = here(plot_dir, "k16_expression_1a-1b_HTRA1.pdf"),
    width = 3,
    height = 4.5
)

## Sample close to median change?

## plot both
set.seed(20221130) ## to make the jitter reproducible
k16_plot_SPARC_HTRA1 <-
    custom_plotExpression(
        spe_k16[, spe_k16$BayesSpace %in% c("Sp16D02", "Sp16D14")],
        genes = c("SPARC", "HTRA1"),
        assay = "logcounts",
        cat = "BayesSpace",
        fill_colors = k16_colors
    ) +
    labs(x = "Pseudobulk k16 Domains")

ggsave(
    k16_plot_SPARC_HTRA1 + theme(axis.text.x = element_text(angle = 90, hjust = 1)),
    filename = here(plot_dir, "k16_expression_1a-1b_SPARC_HTRA1.png"),
    height = 3,
    width = 11
)


## 6a
set.seed(20221130) ## to make the jitter reproducible
k16_6a_plot <- custom_plotExpression(
    spe_k16,
    genes = c("SMIM32", "DACH1", "KIF1A", "GALNT14"),
    assay = "logcounts",
    cat = "BayesSpace",
    fill_colors = k16_colors
)

ggsave(
    k16_6a_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = NULL),
    filename = here(plot_dir, "k16_expression_6a.png")
)

## 6b
set.seed(20221130) ## to make the jitter reproducible
k16_6b_plot <- custom_plotExpression(
    spe_k16,
    genes = c("KRT17", "DIRAS2", "SEMA3E"),
    assay = "logcounts",
    cat = "BayesSpace",
    fill_colors = k16_colors
)

ggsave(
    k16_6b_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = NULL),
    filename = here(plot_dir, "k16_expression_6b.png"),
    width = 10,
    height = 5
)

#### PCA plots ####
spe_k16_sub <- spe_k16[, spe_k16$BayesSpace %in% c("Sp16D02", "Sp16D14")]

spe_k16_PCA <- as.data.frame(reducedDims(spe_k16_sub)$runPCA[, c("PC2", "PC3")])
spe_k16_PCA$k16 <- spe_k16_sub$BayesSpace

head(spe_k16_PCA)
dim(spe_k16_PCA)

# pca_name <- paste0(colnames(reducedDims(spe_k16_sub)$runPCA), " (",metadata(spe_k16_sub)$PCA_var_explained,"%)")

pca_k16 <- ggplot(spe_k16_PCA, aes(x = PC3, y = PC2, fill = k16)) +
    geom_point(shape = 21, size = 1.7) +
    theme_bw() +
    scale_fill_manual(values = k16_colors) +
    labs(x = "PC3 (7%)", y = "PC2 (11%)") +
    theme(legend.position = "top", )

ggsave(pca_k16, filename = here(plot_dir, "pca_k16_2v3.png"), height = 3.4, width = 3)
ggsave(pca_k16, filename = here(plot_dir, "pca_k16_2v3.pdf"), height = 3.4, width = 3)

# sgejobs::job_single('04_select_layer_DE_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript 04_select_layer_DE_plots.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
