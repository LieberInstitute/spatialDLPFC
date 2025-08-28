library(SpatialExperiment)
library(sessioninfo)
library(MOFAcellulaR)
library(MOFA2)
library(here)
library(tidyverse)
library(GGally)

spe_path = here('processed-data', 'MFA', 'pseudobulk_rds', 'visium_pb.rds')
sce_path = here('processed-data', 'MFA', 'pseudobulk_rds', 'sn_fine_pb.rds')
out_path = here('processed-data', 'MFA', 'models', 'DLPFC.rds')
plot_path = here('plots', 'MFA', 'DLPFC_heatmap.pdf')
num_factors = 5

dir.create(dirname(plot_path), showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)
set.seed(1)

################################################################################
#   Combine Visium and snRNA-seq objects
################################################################################

#   Load objects as SingleCellExperiments
spe = as(readRDS(spe_path), "SingleCellExperiment")
int_colData(spe)$spatialCoords = NULL
sce = readRDS(sce_path)

spe$cluster = spe$BayesSpace_harmony_09

#   Ensure genes match
shared_genes = intersect(rownames(spe), rownames(sce))
spe = spe[shared_genes, ]
sce = sce[shared_genes, ]
rowRanges(sce) = rowRanges(spe)

#   Simplify SPE object in preparation for combining with SCE
colData(spe) = colData(spe) |>
    as_tibble() |>
    dplyr::rename(donor = subject, cluster = BayesSpace_harmony_09) |>
    select(donor, cluster, age, sex, ncells) |>
    mutate(sample_id = paste(donor, cluster, sep = '_')) |>
    DataFrame()
colnames(spe) = spe$sample_id
reducedDims(spe) = list()
metadata(spe) = list()

#   Simplify SCE object in preparation for combining with SPE
colData(sce) = colData(sce) |>
    as_tibble() |>
    dplyr::rename(donor = BrNum, cluster = cellType_layer) |>
    select(donor, cluster, age, sex, ncells) |>
    mutate(sample_id = paste(donor, cluster, sep = '_')) |>
    DataFrame()
colnames(sce) = sce$sample_id
reducedDims(sce) = list()
metadata(sce) = list()

#   Gather some donor-level info for later
pd = spe[, match(unique(spe$donor), spe$donor)] |>
    colData() |>
    as_tibble() |>
    select(donor, age, sex)

#   Merge Visium and snRNA-seq objects
sce = cbind(spe, sce)

################################################################################
#   Preprocess expression and create a MOFA object
################################################################################

mofa = create_init_exp(
        counts = assays(sce)$counts, coldata = as.data.frame(colData(sce))
    ) |>
    filt_profiles(
        #   Retain all clusters initially
        cts = unique(sce$cluster),
        ncells = 0,
        counts_col = "ncells",
        ct_col = "cluster"
    ) |>
    #   Drop clusters seen in very few samples
    filt_views_bysamples(nsamples = 2) |>
    #   Require a gene to have 5 counts in at least 25% of samples
    filt_gex_byexpr(min.count = 5, min.prop = 0.25) |>
    #   Drop clusters having below 15 genes at this point
    filt_views_bygenes(ngenes = 15) |>
    #   Drop samples below 90% coverage of genes
    filt_samples_bycov(prop_coverage = 0.9) |>
    #   Normalize expression
    tmm_trns(scale_factor = 1000000) |>
    #   Filter to HVGs only
    filt_gex_byhvg(prior_hvg = NULL, var.threshold = 0) |>
    #   Again, drop clusters having below 15 genes at this point
    filt_views_bygenes(ngenes = 15) |>
    pb_dat2MOFA(sample_column = 'donor') |>
    create_mofa()

#   Metadata seems to need to be manually added  
samples_metadata(mofa) = left_join(
    samples_metadata(mofa), pd, by = c("sample" = "donor")
)

################################################################################
#   Fit the MOFA model
################################################################################

#   Use mostly default options except where overriden in
#   https://saezlab.github.io/MOFAcellulaR/articles/get-started.html#fitting-a-mofa-model.
data_opts = get_default_data_options(mofa)
train_opts = get_default_training_options(mofa)
model_opts = get_default_model_options(mofa)
model_opts$spikeslab_weights = FALSE 
model_opts$num_factors = num_factors

mofa = prepare_mofa(
    object = mofa,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
)

model = run_mofa(mofa, out_path, use_basilisk = TRUE)

################################################################################
#    Exploratory plots
################################################################################

#   Test association of various donor-level covariates with each factor
assoc_list = list()
    assoc_list[['sex']] = get_associations(
        model = model,
        metadata = samples_metadata(model),
        sample_id_column = "sample",
        test_variable = "sex",
        test_type = "categorical",
        group = FALSE
    )
assoc_list[['age']] = get_associations(
    model = model,
    metadata = samples_metadata(model),
    sample_id_column = "sample",
    test_variable = "age",
    test_type = "continuous",
    group = FALSE
)

#   Plot a heatmap of summary results, labeling with covariates of interest
p = plot_MOFA_hmap(
    model = model,
    group = FALSE,
    metadata = samples_metadata(model),
    sample_id_column = "sample",
    sample_anns = c("age", "sex"),
    assoc_list = assoc_list
)
pdf(
    plot_path,
    width = as.integer(round(num_factors / 3) + 5),
    height = as.integer(round(length(unique(sce$cluster)) / 4) + 4)
)
print(p)
dev.off()

session_info()
