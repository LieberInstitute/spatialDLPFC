library("spatialLIBD")
library("scran")
library("here")
library("sessioninfo")

## Load SPE data
load(here("processed-data", "rdata", "spe", "spe_final.Rdata"),
    verbose = TRUE)

## Compute the low library size thresholds in a sample-aware way
spe$scran_low_lib_size <-
    factor(
        isOutlier(
            spe$sum_umi,
            type = "lower",
            log = TRUE,
            batch = spe$sample_id
        ),
        levels = c("TRUE", "FALSE")
    )

## Low library size spots per sample
with(colData(spe), addmargins(table(sample_id, scran_low_lib_size)))
#               scran_low_lib_size
# sample_id        TRUE  FALSE    Sum
#   Br2720_ant_2    121   3046   3167
#   Br2720_mid       29   1796   1825
#   Br2720_post      61   4623   4684
#   Br2743_ant       76   4068   4144
#   Br2743_mid      174   4072   4246
#   Br2743_post     309   3553   3862
#   Br3942_ant      199   3842   4041
#   Br3942_mid      557   3371   3928
#   Br3942_post     243   4157   4400
#   Br6423_ant      118   3788   3906
#   Br6423_mid      112   3872   3984
#   Br6423_post      32   3819   3851
#   Br6432_ant_2     79   3905   3984
#   Br6432_mid      142   3369   3511
#   Br6432_post     273   2614   2887
#   Br6471_ant       52   3135   3187
#   Br6471_mid       65   4476   4541
#   Br6471_post      31   4385   4416
#   Br6522_ant       35   4263   4298
#   Br6522_mid       43   3724   3767
#   Br6522_post      49   3861   3910
#   Br8325_ant       87   3442   3529
#   Br8325_mid_2    320   3840   4160
#   Br8325_post     223   4176   4399
#   Br8492_ant       84   4709   4793
#   Br8492_mid      566   3869   4435
#   Br8492_post      51   4567   4618
#   Br8667_ant       37   3621   3658
#   Br8667_mid      326   3939   4265
#   Br8667_post     372   4025   4397
#   Sum            4866 113927 118793

# vis_clus(
#     spe = spe,
#     sampleid = "Br6471_post",
#     clustervar = "scran_low_lib_size",
#     colors = c("FALSE" = "grey90", "TRUE" = "orange"),
#     spatial = FALSE,
#     point_size = 2
# )
#
# vis_clus(
#     spe = spe,
#     sampleid = "Br8492_mid",
#     clustervar = "scran_low_lib_size",
#     colors = c("FALSE" = "grey90", "TRUE" = "orange"),
#     spatial = FALSE,
#     point_size = 2
# )
#
# vis_clus(
#     spe = spe,
#     sampleid = "Br3942_post",
#     clustervar = "scran_low_lib_size",
#     colors = c("FALSE" = "grey90", "TRUE" = "orange"),
#     spatial = FALSE,
#     point_size = 2
# )

## Load current clustering results
spe <- cluster_import(spe,
    here("processed-data", "rdata", "spe", "clustering_results"),
    "")

## One idea: look at BayesSpace_k7_c6 and find the low library size
## spots inside of it, but with a nmads threshold of 1 instead of the
## default 3.
spe$scran_low_lib_size_wm <- rep(NA, ncol(spe))
spe$scran_low_lib_size_wm[spe$spatial.cluster == 6] <-
    isOutlier(
        spe$sum_umi[spe$spatial.cluster == 6],
        type = "lower",
        log = TRUE,
        batch = spe$sample_id[spe$spatial.cluster == 6],
        nmads = 1
    )
spe$scran_low_lib_size_wm <- factor(spe$scran_low_lib_size_wm,
    levels = c("TRUE", "FALSE", "NA"))
spe$scran_low_lib_size_wm[is.na(spe$scran_low_lib_size_wm)] <- "NA"

## Look at how many of spots we would then drop per sample
with(colData(spe), addmargins(table(sample_id, scran_low_lib_size_wm, useNA = 'ifany')))
#               scran_low_lib_size_wm
# sample_id        TRUE  FALSE   <NA>    Sum
#   Br2720_ant_2    288    882   1997   3167
#   Br2720_mid        7     16   1802   1825
#   Br2720_post      15     84   4585   4684
#   Br2743_ant       83    480   3581   4144
#   Br2743_mid       24    106   4116   4246
#   Br2743_post      28    271   3563   3862
#   Br3942_ant       16    127   3898   4041
#   Br3942_mid       72    447   3409   3928
#   Br3942_post      51    134   4215   4400
#   Br6423_ant      140    639   3127   3906
#   Br6423_mid       31    129   3824   3984
#   Br6423_post      28    158   3665   3851
#   Br6432_ant_2    103    293   3588   3984
#   Br6432_mid       21     76   3414   3511
#   Br6432_post      26    132   2729   2887
#   Br6471_ant       21    101   3065   3187
#   Br6471_mid        7     49   4485   4541
#   Br6471_post      37    154   4225   4416
#   Br6522_ant       17     59   4222   4298
#   Br6522_mid       19     88   3660   3767
#   Br6522_post      21     79   3810   3910
#   Br8325_ant       17     75   3437   3529
#   Br8325_mid_2     40    271   3849   4160
#   Br8325_post      28    253   4118   4399
#   Br8492_ant      145    496   4152   4793
#   Br8492_mid       83    361   3991   4435
#   Br8492_post     124    425   4069   4618
#   Br8667_ant       26    105   3527   3658
#   Br8667_mid        9     68   4188   4265
#   Br8667_post      10     47   4340   4397
#   Sum            1537   6605 110651 118793

# vis_clus(
#     spe = spe,
#     sampleid = "Br6471_post",
#     clustervar = "scran_low_lib_size_wm",
#     colors = c(
#         "FALSE" = "grey20",
#         "TRUE" = "orange",
#         "NA" = "grey80"
#     ),
#     spatial = FALSE,
#     point_size = 2
# )
#
# vis_clus(
#     spe = spe,
#     sampleid = "Br8492_mid",
#     clustervar = "scran_low_lib_size_wm",
#     colors = c(
#         "FALSE" = "grey20",
#         "TRUE" = "orange",
#         "NA" = "grey80"
#     ),
#     spatial = FALSE,
#     point_size = 2
# )
#
# vis_clus(
#     spe = spe,
#     sampleid = "Br3942_mid",
#     clustervar = "scran_low_lib_size_wm",
#     colors = c(
#         "FALSE" = "grey20",
#         "TRUE" = "orange",
#         "NA" = "grey80"
#     ),
#     spatial = FALSE,
#     point_size = 2
# )


# vis_clus(
#     spe = spe,
#     sampleid = "Br3942_mid",
#     clustervar = "bayesSpace_harmony_15",
#     colors = setNames(Polychrome::palette36.colors(15), 1:15),
#     spatial = FALSE,
#     point_size = 2
# )
#
# vis_clus(
#     spe = spe,
#     sampleid = "Br3942_post",
#     clustervar = "scran_low_lib_size_wm",
#     colors = c(
#         "FALSE" = "grey20",
#         "TRUE" = "orange",
#         "NA" = "grey80"
#     ),
#     spatial = FALSE,
#     point_size = 2
# )
#
# vis_clus(
#     spe = spe,
#     sampleid = "Br3942_post",
#     clustervar = "bayesSpace_harmony_15",
#     colors = setNames(Polychrome::palette36.colors(15), 1:15),
#     spatial = FALSE,
#     point_size = 2
# )

## Let's look at cluster 12 from BayesSpace k = 15
spe$BayesSpace_k15_c12 <- factor(spe$bayesSpace_harmony_15 == 12,
    levels = c("TRUE", "FALSE"))
vis_grid_clus(
    spe = spe,
    clustervar = "BayesSpace_k15_c12",
    pdf = here("plots", paste0("vis_clus_BayesSpace_k15_c12.pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange"),
    spatial = FALSE,
    point_size = 2,
    height = 24,
    width = 90
)

## Sample-aware low-library size
vis_grid_clus(
    spe = spe,
    clustervar = "scran_low_lib_size",
    pdf = here("plots", paste0("vis_clus_sample_aware_low_lib_size.pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange"),
    spatial = FALSE,
    point_size = 2,
    height = 24,
    width = 90
)

## Low-library size within BayesSpace k = 7 cluster 6
vis_grid_clus(
    spe = spe,
    clustervar = "scran_low_lib_size_wm",
    pdf = here(
        "plots",
        paste0("vis_clus_BayesSpace_k7_c6_low_lib_size.pdf")
    ),
    sort_clust = FALSE,
    colors = c(
        "FALSE" = "grey20",
        "TRUE" = "orange",
        "NA" = "grey80"
    ),
    spatial = FALSE,
    point_size = 2,
    height = 24,
    width = 90
)

## BayesSpace k = 15 clusters in colors
vis_grid_clus(
    spe = spe,
    clustervar = "bayesSpace_harmony_15",
    pdf = here("plots", paste0("vis_clus_BayesSpace_k15.pdf")),
    sort_clust = FALSE,
    colors = setNames(Polychrome::palette36.colors(15), 1:15),
    spatial = FALSE,
    point_size = 2,
    height = 24,
    width = 90
)

## BayesSpace k = 7 clusters in colors
vis_grid_clus(
    spe = spe,
    clustervar = "spatial.cluster",
    pdf = here("plots", paste0("vis_clus_BayesSpace_k7.pdf")),
    sort_clust = FALSE,
    colors = setNames(Polychrome::palette36.colors(15), 1:7),
    spatial = FALSE,
    point_size = 2,
    height = 24,
    width = 90
)

## We noticed that k15 cluster 12 was mostly a subset of k7 cluster 6, which
## is what we wanted to find. So here we were thinking that we could use
## k = 15 cluster 12 as our indicator of "bad spots" that we would need to
## drop.
with(colData(spe), addmargins(
    table(
        BayesSpace_k7_c6 = spatial.cluster == 6,
        BayesSpace_k15_c12 = bayesSpace_harmony_15 == 12
    )
))
#                 BayesSpace_k15_c12
# BayesSpace_k7_c6  FALSE   TRUE    Sum
#            FALSE 110544    107 110651
#            TRUE    4745   3397   8142
#            Sum   115289   3504 118793

## Later we noticed that k = 15 cluster 12 would potentially drop some
## white matter spots in a few samples, plus it was maybe not being too
## strict on dropping bad spots in other samples. Look at the "blue" vs "pink"
## circles on one of the screenshots from today.
## With this, we are now leaning towards dropping the low library size spots
## in a sample-aware way.
with(colData(spe), addmargins(
    table(
        sample_aware_low_lib_size = scran_low_lib_size == "TRUE",
        BayesSpace_k15_c12 = bayesSpace_harmony_15 == 12
    )
))
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE 112809   1118 113927
#                     TRUE    2480   2386   4866
#                     Sum   115289   3504 118793

## Same as the above table, but broken down by sample.
with(colData(spe), addmargins(
    table(
        sample_aware_low_lib_size = scran_low_lib_size == "TRUE",
        BayesSpace_k15_c12 = bayesSpace_harmony_15 == 12,
        sample_id
    )
))
# , , sample_id = Br2720_ant_2
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   2755    291   3046
#                     TRUE       1    120    121
#                     Sum     2756    411   3167
#
# , , sample_id = Br2720_mid
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   1796      0   1796
#                     TRUE      14     15     29
#                     Sum     1810     15   1825
#
# , , sample_id = Br2720_post
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   4618      5   4623
#                     TRUE      40     21     61
#                     Sum     4658     26   4684
#
# , , sample_id = Br2743_ant
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3871    197   4068
#                     TRUE      10     66     76
#                     Sum     3881    263   4144
#
# , , sample_id = Br2743_mid
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   4071      1   4072
#                     TRUE     109     65    174
#                     Sum     4180     66   4246
#
# , , sample_id = Br2743_post
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3553      0   3553
#                     TRUE     108    201    309
#                     Sum     3661    201   3862
#
# , , sample_id = Br3942_ant
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3842      0   3842
#                     TRUE     123     76    199
#                     Sum     3965     76   4041
#
# , , sample_id = Br3942_mid
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3369      2   3371
#                     TRUE     163    394    557
#                     Sum     3532    396   3928
#
# , , sample_id = Br3942_post
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   4157      0   4157
#                     TRUE     162     81    243
#                     Sum     4319     81   4400
#
# , , sample_id = Br6423_ant
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3555    233   3788
#                     TRUE      10    108    118
#                     Sum     3565    341   3906
#
# , , sample_id = Br6423_mid
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3866      6   3872
#                     TRUE      66     46    112
#                     Sum     3932     52   3984
#
# , , sample_id = Br6423_post
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3791     28   3819
#                     TRUE      11     21     32
#                     Sum     3802     49   3851
#
# , , sample_id = Br6432_ant_2
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3837     68   3905
#                     TRUE      32     47     79
#                     Sum     3869    115   3984
#
# , , sample_id = Br6432_mid
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3369      0   3369
#                     TRUE      96     46    142
#                     Sum     3465     46   3511
#
# , , sample_id = Br6432_post
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   2614      0   2614
#                     TRUE     194     79    273
#                     Sum     2808     79   2887
#
# , , sample_id = Br6471_ant
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3100     35   3135
#                     TRUE       5     47     52
#                     Sum     3105     82   3187
#
# , , sample_id = Br6471_mid
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   4472      4   4476
#                     TRUE      30     35     65
#                     Sum     4502     39   4541
#
# , , sample_id = Br6471_post
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   4325     60   4385
#                     TRUE      17     14     31
#                     Sum     4342     74   4416
#
# , , sample_id = Br6522_ant
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   4255      8   4263
#                     TRUE      17     18     35
#                     Sum     4272     26   4298
#
# , , sample_id = Br6522_mid
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3716      8   3724
#                     TRUE      21     22     43
#                     Sum     3737     30   3767
#
# , , sample_id = Br6522_post
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3850     11   3861
#                     TRUE      32     17     49
#                     Sum     3882     28   3910
#
# , , sample_id = Br8325_ant
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3442      0   3442
#                     TRUE      39     48     87
#                     Sum     3481     48   3529
#
# , , sample_id = Br8325_mid_2
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3834      6   3840
#                     TRUE     113    207    320
#                     Sum     3947    213   4160
#
# , , sample_id = Br8325_post
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   4164     12   4176
#                     TRUE      78    145    223
#                     Sum     4242    157   4399
#
# , , sample_id = Br8492_ant
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   4671     38   4709
#                     TRUE      60     24     84
#                     Sum     4731     62   4793
#
# , , sample_id = Br8492_mid
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3868      1   3869
#                     TRUE     287    279    566
#                     Sum     4155    280   4435
#
# , , sample_id = Br8492_post
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   4503     64   4567
#                     TRUE      23     28     51
#                     Sum     4526     92   4618
#
# , , sample_id = Br8667_ant
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3584     37   3621
#                     TRUE       3     34     37
#                     Sum     3587     71   3658
#
# , , sample_id = Br8667_mid
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   3939      0   3939
#                     TRUE     276     50    326
#                     Sum     4215     50   4265
#
# , , sample_id = Br8667_post
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE   4022      3   4025
#                     TRUE     340     32    372
#                     Sum     4362     35   4397
#
# , , sample_id = Sum
#
#                          BayesSpace_k15_c12
# sample_aware_low_lib_size  FALSE   TRUE    Sum
#                     FALSE 112809   1118 113927
#                     TRUE    2480   2386   4866
#                     Sum   115289   3504 118793

## Exploring whether we could justify dropping cluster 12 from k = 15
tapply(spe$sum_umi, spe$bayesSpace_harmony_15, summary)
# $`1`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#      78    2079    2994    3403    4292   46881
#
# $`2`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     124    2120    2971    3411    4220   26262
#
# $`3`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#      26     974    1439    1651    2071   12164
#
# $`4`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#      59    2084    2959    3433    4281   26418
#
# $`5`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    33.0   475.0   657.5   773.8   928.8  3463.0
#
# $`6`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#      99     998    1393    1614    1984   13059
#
# $`7`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    19.0   312.0   515.0   620.9   818.8  3933.0
#
# $`9`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     214    2071    2867    3320    4093   18305
#
# $`10`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   212.0   917.5  1464.0  1704.8  2148.0  9612.0
#
# $`11`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     209    1424    2048    2411    2989   11888
#
# $`12`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     1.0    41.0   123.0   190.4   232.0  2458.0
#
# $`13`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    11.0   342.0   722.0   779.3  1106.0  4244.0
#
# $`14`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     156    1707    2388    2677    3310   18025
#
# $`15`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    70.0   618.0   856.0   990.8  1224.0  6940.0

library("ggplot2")
pdf(here("plots", "BayesSpace_k15_library_size_boxplots.pdf"), width = 14)
ggplot(data.frame(
    cluster = spe$bayesSpace_harmony_15,
    lib_size = log10(spe$sum_umi + 1)
), aes(x = cluster, y = lib_size, group = cluster)) + geom_boxplot() +
    ylab("log10(lib_size = 1)") +
    theme_bw(base_size = 20)
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.2 (2021-11-01)
#  os       macOS Monterey 12.2.1
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2022-02-24
#  rstudio  2021.09.2+382 Ghost Orchid (desktop)
#  pandoc   NA
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version    date (UTC) lib source
#  AnnotationDbi            1.56.2     2021-11-09 [1] Bioconductor
#  AnnotationHub            3.2.1      2022-01-23 [1] Bioconductor
#  assertthat               0.2.1      2019-03-21 [1] CRAN (R 4.1.0)
#  attempt                  0.3.1      2020-05-03 [1] CRAN (R 4.1.0)
#  beachmat                 2.10.0     2021-10-26 [1] Bioconductor
#  beeswarm                 0.4.0      2021-06-01 [1] CRAN (R 4.1.0)
#  benchmarkme              1.0.7      2021-03-21 [1] CRAN (R 4.1.0)
#  benchmarkmeData          1.0.4      2020-04-23 [1] CRAN (R 4.1.0)
#  Biobase                * 2.54.0     2021-10-26 [1] Bioconductor
#  BiocFileCache            2.2.1      2022-01-23 [1] Bioconductor
#  BiocGenerics           * 0.40.0     2021-10-26 [1] Bioconductor
#  BiocIO                   1.4.0      2021-10-26 [1] Bioconductor
#  BiocManager              1.30.16    2021-06-15 [1] CRAN (R 4.1.0)
#  BiocNeighbors            1.12.0     2021-10-26 [1] Bioconductor
#  BiocParallel             1.28.3     2021-12-09 [1] Bioconductor
#  BiocSingular             1.10.0     2021-10-26 [1] Bioconductor
#  BiocVersion              3.14.0     2021-05-19 [1] Bioconductor
#  Biostrings               2.62.0     2021-10-26 [1] Bioconductor
#  bit                      4.0.4      2020-08-04 [1] CRAN (R 4.1.1)
#  bit64                    4.0.5      2020-08-30 [1] CRAN (R 4.1.0)
#  bitops                   1.0-7      2021-04-24 [1] CRAN (R 4.1.0)
#  blob                     1.2.2      2021-07-23 [1] CRAN (R 4.1.0)
#  bluster                  1.4.0      2021-10-26 [1] Bioconductor
#  brio                     1.1.3      2021-11-30 [1] CRAN (R 4.1.1)
#  bslib                    0.3.1      2021-10-06 [1] CRAN (R 4.1.1)
#  cachem                   1.0.6      2021-08-19 [1] CRAN (R 4.1.1)
#  callr                    3.7.0      2021-04-20 [1] CRAN (R 4.1.0)
#  cli                      3.1.1      2022-01-20 [1] CRAN (R 4.1.1)
#  cluster                  2.1.2      2021-04-17 [1] CRAN (R 4.1.2)
#  codetools                0.2-18     2020-11-04 [1] CRAN (R 4.1.2)
#  colorout                 1.2-2      2022-02-07 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace               2.0-2      2021-06-24 [1] CRAN (R 4.1.1)
#  config                   0.3.1      2020-12-17 [1] CRAN (R 4.1.0)
#  cowplot                  1.1.1      2020-12-30 [1] CRAN (R 4.1.1)
#  crayon                   1.4.2      2021-10-29 [1] CRAN (R 4.1.1)
#  curl                     4.3.2      2021-06-23 [1] CRAN (R 4.1.0)
#  data.table               1.14.2     2021-09-27 [1] CRAN (R 4.1.1)
#  DBI                      1.1.2      2021-12-20 [1] CRAN (R 4.1.1)
#  dbplyr                   2.1.1      2021-04-06 [1] CRAN (R 4.1.0)
#  DelayedArray             0.20.0     2021-10-26 [1] Bioconductor
#  DelayedMatrixStats       1.16.0     2021-10-26 [1] Bioconductor
#  desc                     1.4.0      2021-09-28 [1] CRAN (R 4.1.1)
#  devtools               * 2.4.3      2021-11-30 [1] CRAN (R 4.1.1)
#  digest                   0.6.29     2021-12-01 [1] CRAN (R 4.1.1)
#  dockerfiler              0.1.4      2021-09-03 [1] CRAN (R 4.1.1)
#  doParallel               1.0.17     2022-02-07 [1] CRAN (R 4.1.2)
#  dotCall64                1.0-1      2021-02-11 [1] CRAN (R 4.1.0)
#  dplyr                    1.0.8      2022-02-08 [1] CRAN (R 4.1.2)
#  dqrng                    0.3.0      2021-05-01 [1] CRAN (R 4.1.0)
#  DropletUtils             1.14.2     2022-01-09 [1] Bioconductor
#  DT                       0.20       2021-11-15 [1] CRAN (R 4.1.1)
#  edgeR                    3.36.0     2021-10-26 [1] Bioconductor
#  ellipsis                 0.3.2      2021-04-29 [1] CRAN (R 4.1.0)
#  ExperimentHub            2.2.1      2022-01-23 [1] Bioconductor
#  fansi                    1.0.2      2022-01-14 [1] CRAN (R 4.1.1)
#  farver                   2.1.0      2021-02-28 [1] CRAN (R 4.1.0)
#  fastmap                  1.1.0      2021-01-25 [1] CRAN (R 4.1.0)
#  fields                   13.3       2021-10-30 [1] CRAN (R 4.1.1)
#  filelock                 1.0.2      2018-10-05 [1] CRAN (R 4.1.0)
#  foreach                  1.5.2      2022-02-02 [1] CRAN (R 4.1.1)
#  fs                       1.5.2      2021-12-08 [1] CRAN (R 4.1.1)
#  generics                 0.1.2      2022-01-31 [1] CRAN (R 4.1.1)
#  GenomeInfoDb           * 1.30.1     2022-01-30 [1] Bioconductor
#  GenomeInfoDbData         1.2.7      2022-02-07 [1] Bioconductor
#  GenomicAlignments        1.30.0     2021-10-26 [1] Bioconductor
#  GenomicRanges          * 1.46.1     2021-11-18 [1] Bioconductor
#  ggbeeswarm               0.6.0      2017-08-07 [1] CRAN (R 4.1.0)
#  ggplot2                  3.3.5      2021-06-25 [1] CRAN (R 4.1.1)
#  ggrepel                  0.9.1      2021-01-15 [1] CRAN (R 4.1.1)
#  glue                     1.6.1      2022-01-22 [1] CRAN (R 4.1.1)
#  golem                    0.3.1      2021-04-17 [1] CRAN (R 4.1.0)
#  gridExtra                2.3        2017-09-09 [1] CRAN (R 4.1.1)
#  gtable                   0.3.0      2019-03-25 [1] CRAN (R 4.1.1)
#  HDF5Array                1.22.1     2021-11-14 [1] Bioconductor
#  here                   * 1.0.1      2020-12-13 [1] CRAN (R 4.1.0)
#  hms                      1.1.1      2021-09-26 [1] CRAN (R 4.1.1)
#  htmltools                0.5.2      2021-08-25 [1] CRAN (R 4.1.1)
#  htmlwidgets              1.5.4      2021-09-08 [1] CRAN (R 4.1.1)
#  httpuv                   1.6.5      2022-01-05 [1] CRAN (R 4.1.1)
#  httr                     1.4.2      2020-07-20 [1] CRAN (R 4.1.0)
#  igraph                   1.2.11     2022-01-04 [1] CRAN (R 4.1.1)
#  interactiveDisplayBase   1.32.0     2021-10-26 [1] Bioconductor
#  IRanges                * 2.28.0     2021-10-26 [1] Bioconductor
#  irlba                    2.3.5      2021-12-06 [1] CRAN (R 4.1.1)
#  iterators                1.0.14     2022-02-05 [1] CRAN (R 4.1.2)
#  jquerylib                0.1.4      2021-04-26 [1] CRAN (R 4.1.0)
#  jsonlite                 1.7.3      2022-01-17 [1] CRAN (R 4.1.1)
#  KEGGREST                 1.34.0     2021-10-26 [1] Bioconductor
#  knitr                    1.37       2021-12-16 [1] CRAN (R 4.1.1)
#  labeling                 0.4.2      2020-10-20 [1] CRAN (R 4.1.0)
#  later                    1.3.0      2021-08-18 [1] CRAN (R 4.1.1)
#  lattice                  0.20-45    2021-09-22 [1] CRAN (R 4.1.2)
#  lazyeval                 0.2.2      2019-03-15 [1] CRAN (R 4.1.0)
#  lifecycle                1.0.1      2021-09-24 [1] CRAN (R 4.1.1)
#  limma                    3.50.0     2021-10-26 [1] Bioconductor
#  lobstr                   1.1.1      2019-07-02 [1] CRAN (R 4.1.0)
#  locfit                   1.5-9.4    2020-03-25 [1] CRAN (R 4.1.0)
#  lubridate                1.8.0      2021-10-07 [1] CRAN (R 4.1.1)
#  magick                   2.7.3      2021-08-18 [1] CRAN (R 4.1.1)
#  magrittr                 2.0.2      2022-01-26 [1] CRAN (R 4.1.1)
#  maps                     3.4.0      2021-09-25 [1] CRAN (R 4.1.1)
#  Matrix                   1.4-0      2021-12-08 [1] CRAN (R 4.1.1)
#  MatrixGenerics         * 1.6.0      2021-10-26 [1] Bioconductor
#  matrixStats            * 0.61.0     2021-09-17 [1] CRAN (R 4.1.1)
#  memoise                  2.0.1      2021-11-26 [1] CRAN (R 4.1.1)
#  metapod                  1.2.0      2021-10-26 [1] Bioconductor
#  mime                     0.12       2021-09-28 [1] CRAN (R 4.1.1)
#  munsell                  0.5.0      2018-06-12 [1] CRAN (R 4.1.0)
#  pillar                   1.7.0      2022-02-01 [1] CRAN (R 4.1.1)
#  pkgbuild                 1.3.1      2021-12-20 [1] CRAN (R 4.1.1)
#  pkgconfig                2.0.3      2019-09-22 [1] CRAN (R 4.1.0)
#  pkgload                  1.2.4      2021-11-30 [1] CRAN (R 4.1.1)
#  plotly                   4.10.0     2021-10-09 [1] CRAN (R 4.1.1)
#  png                      0.1-7      2013-12-03 [1] CRAN (R 4.1.0)
#  Polychrome               1.3.1      2021-07-16 [1] CRAN (R 4.1.1)
#  prettyunits              1.1.1      2020-01-24 [1] CRAN (R 4.1.0)
#  processx                 3.5.2      2021-04-30 [1] CRAN (R 4.1.0)
#  promises                 1.2.0.1    2021-02-11 [1] CRAN (R 4.1.0)
#  prompt                   1.0.1      2022-02-07 [1] Github (gaborcsardi/prompt@7ef0f2e)
#  ps                       1.6.0      2021-02-28 [1] CRAN (R 4.1.0)
#  purrr                    0.3.4      2020-04-17 [1] CRAN (R 4.1.0)
#  R.methodsS3              1.8.1      2020-08-26 [1] CRAN (R 4.1.0)
#  R.oo                     1.24.0     2020-08-26 [1] CRAN (R 4.1.0)
#  R.utils                  2.11.0     2021-09-26 [1] CRAN (R 4.1.1)
#  R6                       2.5.1      2021-08-19 [1] CRAN (R 4.1.1)
#  rappdirs                 0.3.3      2021-01-31 [1] CRAN (R 4.1.0)
#  RColorBrewer             1.1-2      2014-12-07 [1] CRAN (R 4.1.0)
#  Rcpp                     1.0.8      2022-01-13 [1] CRAN (R 4.1.1)
#  RCurl                    1.98-1.6   2022-02-08 [1] CRAN (R 4.1.2)
#  remotes                  2.4.2      2021-11-30 [1] CRAN (R 4.1.1)
#  restfulr                 0.0.13     2017-08-06 [1] CRAN (R 4.1.0)
#  rhdf5                    2.38.0     2021-10-26 [1] Bioconductor
#  rhdf5filters             1.6.0      2021-10-26 [1] Bioconductor
#  Rhdf5lib                 1.16.0     2021-10-26 [1] Bioconductor
#  rjson                    0.2.21     2022-01-09 [1] CRAN (R 4.1.1)
#  rlang                    1.0.1      2022-02-03 [1] CRAN (R 4.1.1)
#  roxygen2                 7.1.2      2021-09-08 [1] CRAN (R 4.1.1)
#  rprojroot                2.0.2      2020-11-15 [1] CRAN (R 4.1.0)
#  Rsamtools                2.10.0     2021-10-26 [1] Bioconductor
#  rsconnect                0.8.25     2021-11-19 [1] CRAN (R 4.1.1)
#  RSQLite                  2.2.9      2021-12-06 [1] CRAN (R 4.1.1)
#  rsthemes                 0.3.1      2022-02-07 [1] Github (gadenbuie/rsthemes@bbe73ca)
#  rstudioapi               0.13       2020-11-12 [1] CRAN (R 4.1.0)
#  rsvd                     1.0.5      2021-04-16 [1] CRAN (R 4.1.0)
#  rtracklayer              1.54.0     2021-10-26 [1] Bioconductor
#  S4Vectors              * 0.32.3     2021-11-21 [1] Bioconductor
#  sass                     0.4.0.9000 2022-02-07 [1] Github (rstudio/sass@f7a9540)
#  ScaledMatrix             1.2.0      2021-10-26 [1] Bioconductor
#  scales                   1.1.1      2020-05-11 [1] CRAN (R 4.1.0)
#  scater                   1.22.0     2021-10-26 [1] Bioconductor
#  scatterplot3d            0.3-41     2018-03-14 [1] CRAN (R 4.1.0)
#  scran                  * 1.22.1     2021-11-14 [1] Bioconductor
#  scuttle                * 1.4.0      2021-10-26 [1] Bioconductor
#  sessioninfo            * 1.2.2      2021-12-06 [1] CRAN (R 4.1.1)
#  shiny                    1.7.1      2021-10-02 [1] CRAN (R 4.1.1)
#  shinyWidgets             0.6.4      2022-02-06 [1] CRAN (R 4.1.2)
#  SingleCellExperiment   * 1.16.0     2021-10-26 [1] Bioconductor
#  spam                     2.8-0      2022-01-06 [1] CRAN (R 4.1.1)
#  sparseMatrixStats        1.6.0      2021-10-26 [1] Bioconductor
#  SpatialExperiment      * 1.4.0      2021-10-26 [1] Bioconductor
#  spatialLIBD            * 1.6.5      2022-01-12 [1] Bioconductor
#  statmod                  1.4.36     2021-05-10 [1] CRAN (R 4.1.0)
#  stringi                  1.7.6      2021-11-29 [1] CRAN (R 4.1.1)
#  stringr                  1.4.0      2019-02-10 [1] CRAN (R 4.1.1)
#  SummarizedExperiment   * 1.24.0     2021-10-26 [1] Bioconductor
#  suncalc                  0.5.0      2019-04-03 [1] CRAN (R 4.1.0)
#  testthat               * 3.1.2      2022-01-20 [1] CRAN (R 4.1.1)
#  tibble                   3.1.6      2021-11-07 [1] CRAN (R 4.1.1)
#  tidyr                    1.2.0      2022-02-01 [1] CRAN (R 4.1.1)
#  tidyselect               1.1.1      2021-04-30 [1] CRAN (R 4.1.0)
#  usethis                * 2.1.5      2021-12-09 [1] CRAN (R 4.1.1)
#  utf8                     1.2.2      2021-07-24 [1] CRAN (R 4.1.0)
#  vctrs                    0.3.8      2021-04-29 [1] CRAN (R 4.1.0)
#  vipor                    0.4.5      2017-03-22 [1] CRAN (R 4.1.0)
#  viridis                  0.6.2      2021-10-13 [1] CRAN (R 4.1.1)
#  viridisLite              0.4.0      2021-04-13 [1] CRAN (R 4.1.0)
#  withr                    2.4.3      2021-11-30 [1] CRAN (R 4.1.1)
#  xfun                     0.29       2021-12-14 [1] CRAN (R 4.1.1)
#  XML                      3.99-0.8   2021-09-17 [1] CRAN (R 4.1.1)
#  xml2                     1.3.3      2021-11-30 [1] CRAN (R 4.1.1)
#  xtable                   1.8-4      2019-04-21 [1] CRAN (R 4.1.0)
#  XVector                  0.34.0     2021-10-26 [1] Bioconductor
#  yaml                     2.2.2      2022-01-25 [1] CRAN (R 4.1.1)
#  zlibbioc                 1.40.0     2021-10-26 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
