## Automatically style the code in this script:
styler::style_file(
    here::here("analysis", "01a_check_scran_quick_cluster.R"),
    transformers = biocthis::bioc_style()
)

## This script requires R 4.1
# module load conda_R/devel

## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")

## vis
library("spatialLIBD")

## Load SPE data
load(here::here("rdata", "spe", "spe.Rdata"), verbose = TRUE)

## Check cluster 89 from before
load(
    file = here::here("analysis", "clusters_nonzero.rda"),
    verbose = TRUE
)

table(spe$scran_quick_cluster == clusters)
## Clusters did change!
# FALSE  TRUE
# 41244  8755

## old cluster "89" from https://github.com/LTLA/scuttle/issues/7
## is not the same one now, despite the number of spots matching
identical(which(spe$scran_quick_cluster == "89"), which(clusters == "89")) # FALSE
stopifnot(identical(sum(spe$scran_quick_cluster == "89"), sum(clusters == "89")))

## Hm... maybe it's because the samples are in a different order!
## due to https://github.com/LieberInstitute/spatialDLPFC/blob/main/analysis/01_build_SPE.R#L72
stopifnot(sum(table(spe$scran_quick_cluster) - table(clusters)) == 0)
## Ok, that matches. Let's check more closely then.
load(
    file = here::here("analysis", "sce_nonzero.rda"),
    verbose = TRUE
)
m <-
    match(spe$key, with(colData(sce_nonzero), paste0(barcode_id, "_", sample_name)))
## Ok, the quick clusters are the same ones!
stopifnot(identical(spe$scran_quick_cluster, clusters[m]))

## Remove the older R objects since we don't need them
rm(clusters, sce_nonzero)

table(spe$scran_quick_cluster == "89")
# FALSE  TRUE
# 49727   272
spe$quick_cluster_89 <-
    factor(spe$scran_quick_cluster == "89", levels = c("TRUE", "FALSE"))

summary(as.data.frame(colData(spe)[spe$scran_quick_cluster == "89", c(
    "scran_low_lib_size",
    "scran_low_n_features",
    "scran_high_subsets_Mito_percent",
    "scran_discard"
)]))
# scran_low_lib_size scran_low_n_features scran_high_subsets_Mito_percent
# TRUE : 62          TRUE : 86            TRUE :  2
# FALSE:210          FALSE:186            FALSE:270
# scran_discard
# TRUE : 87
# FALSE:185
vis_grid_clus(
    spe = spe,
    clustervar = "quick_cluster_89",
    pdf = here::here("plots", paste0("scuttle_", "quick_cluster_89", ".pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange")
)

## Now check cluster 64 since that's the one giving the error at
## https://github.com/LieberInstitute/spatialDLPFC/blob/main/analysis/01_build_SPE.R#L338
spe$quick_cluster_64 <-
    factor(spe$scran_quick_cluster == "64", levels = c("TRUE", "FALSE"))
summary(as.data.frame(colData(spe)[spe$scran_quick_cluster == "64", c(
    "scran_low_lib_size",
    "scran_low_n_features",
    "scran_high_subsets_Mito_percent",
    "scran_discard"
)]))
# scran_low_lib_size scran_low_n_features scran_high_subsets_Mito_percent
# TRUE :  0          TRUE :  0            TRUE : 56
# FALSE:542          FALSE:542            FALSE:486
# scran_discard
# TRUE : 56
# FALSE:486
vis_grid_clus(
    spe = spe,
    clustervar = "quick_cluster_64",
    pdf = here::here("plots", paste0("scuttle_", "quick_cluster_64", ".pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange")
)

## currently scran::computeSumFactors() is simply scuttle::computePooledFactors()
# function (...)
# {
#     computePooledFactors(...)
# }
# <bytecode: 0x55d9eaf11090>
# <environment: namespace:scran>
##
## computePooledFactors() calls .calculate_pooled_factors()
## https://github.com/LTLA/scuttle/blob/3cb28efbf237ecb9b7a1973ab7e8957371260a32/R/pooledSizeFactors.R#L507
##

x <- counts(spe) ## to match the default
# > class(x)
# [1] "dgCMatrix"
# attr(,"package")
# [1] "Matrix"
sizes=seq(21, 101, 5)
clusters = spe$scran_quick_cluster
ref.clust=NULL
max.cluster.size=3000
positive=TRUE
scaling=NULL
min.mean=NULL
subset.row=NULL
BPPARAM = BiocParallel::MulticoreParam(4)

## Next I adapted the code from
## https://github.com/LTLA/scuttle/blob/3cb28efbf237ecb9b7a1973ab7e8957371260a32/R/pooledSizeFactors.R#L208-L263
## using scuttle::: for the internal functions and BiocParallel:: where necessary

ncells <- ncol(x)
if (is.null(clusters)) {
    clusters <- integer(ncells)
}
clusters <- scuttle:::.limit_cluster_size(clusters, max.cluster.size)
# > table(clusters)
# clusters
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
#  466  634  848  841  513  763  462  544  176  514  760  606  123  576  150  893
#   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32
#  920  161  380  675  877 1309 1329  679 1476  776  284  856  468  462  998  254
#   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  148  697  424  662  579  190  161  235  199  183  111  197  290  229  327  935
#   49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64
#  537  542  366  207  461  104  110  166  272  644  992  393  514  542  269  159
#   65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80
#  650  245  674  142  271  962  359  308  435  198  268  178  586  818  884  684
#   81   82   83   84   85   86   87   88   89   90   91   92   93   94   95
#  271  259  186  432  102 1193  860 1219  579 1141  760  787  896  821  213
### OHHH!!! the clusters change internally!
# > table(clusters) - table(spe$scran_quick_cluster)
# clusters
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
# -297  168  214   -7 -328  250 -394   82 -292 -262 -238  458 -131  292 -158  695
#   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32
#  675 -801 -294  240  518 1041 1058  537  716  170 -260  706  345    0  484 -322
#   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  -28  407  263  479 -118 -472 -418   36    9  -52  -86 -227  179   43 -105  676
#   49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64
# -147 -276   95 -677 -125  -74    8    5 -648 -249  612 -484 -161    5  165 -383
#   65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80
# -285 -216  347  -65   42  596  193  198 -758 -662 -951 -401  -93 -491 -592 -645
#   81   82   83   84   85   86   87   88   89   90   91   92   93   94   95
# -379 -385 -207 -560  -57  651  346  950  307  245  -27  -34 -245  608 -547

## Looking more closely at
## https://github.com/LTLA/scuttle/blob/3cb28efbf237ecb9b7a1973ab7e8957371260a32/R/pooledSizeFactors.R#L451-L481
# > any(table(spe$scran_quick_cluster) > max.cluster.size)
# [1] FALSE
## This line https://github.com/LTLA/scuttle/blob/3cb28efbf237ecb9b7a1973ab7e8957371260a32/R/pooledSizeFactors.R#L461
## likely changes the order of the clusters (and thus meaning) given that the
## first unique cluster is cluster number 2
# > unique(spe$scran_quick_cluster)
#  [1] 2  3  4  5  6  1  30 27 33 31 25 26 29 32 28 58 57 56 59 61 60 78 80 77 79 10 14 7  9  8  11 13 12 37 44 38 39 41
# [39] 35 42 40 36 45 43 34 69 67 65 62 64 70 68 66 63 72 71 89 82 84 83 87 86 88 85 81 17 19 24 23 18 21 15 20 16 22 54
# [77] 53 50 52 49 51 48 46 47 55 73 74 75 76 93 95 91 90 92 94
# 95 Levels: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 ... 95

## Ok, they are the same ones, just re-named!
m_table <- data.frame(
    original = seq_len(length(unique(spe$scran_quick_cluster)))
)
m_table$new <- match(m_table$original, as.integer(unique(spe$scran_quick_cluster)))
identical(factor(m_table$new[as.integer(spe$scran_quick_cluster)]), clusters)
# [1] TRUE
m_table[m_table$new == 64, ]
#    original new
# 85       85  64

spe$quick_cluster_85 <-
    factor(spe$scran_quick_cluster == "85", levels = c("TRUE", "FALSE"))
summary(as.data.frame(colData(spe)[spe$scran_quick_cluster == "85", c(
    "scran_low_lib_size",
    "scran_low_n_features",
    "scran_high_subsets_Mito_percent",
    "scran_discard"
)]))
## OK! Finally this cluster does look bad! I mean, it was weird the previous clusters didn't seem bad before!
# scran_low_lib_size scran_low_n_features scran_high_subsets_Mito_percent scran_discard
#  TRUE :159          TRUE :159            TRUE : 55                       TRUE :159
#  FALSE:  0          FALSE:  0            FALSE:104                       FALSE:  0
vis_grid_clus(
    spe = spe,
    clustervar = "quick_cluster_85",
    pdf = here::here("plots", paste0("scuttle_", "quick_cluster_85", ".pdf")),
    sort_clust = FALSE,
    colors = c("FALSE" = "grey90", "TRUE" = "orange")
)

## I'll tell Aaron about this
## Done at https://github.com/LTLA/scuttle/issues/8
.limit_cluster_size.fix <- function(clusters, max.size)
{
    if (is.null(max.size)) {
        return(clusters)
    }

    new.clusters <- integer(length(clusters))
    counter <- 1L
    for (id in sort(unique(clusters))) {
        current <- id==clusters
        ncells <- sum(current)

        if (ncells <= max.size) {
            new.clusters[current] <- counter
            counter <- counter+1L
            next
        }

        # Size of output clusters is max.size * N / ceil(N), where N = ncells/max.size.
        # This is minimal at the smallest N > 1, where output clusters are at least max.size/2.
        # Thus, we need max.size/2 >= min.size to guarantee that the output clusters are >= min.size.
        mult <- ceiling(ncells/max.size)
        realloc <- rep(seq_len(mult) - 1L + counter, length.out=ncells)
        new.clusters[current] <- realloc
        counter <- counter + mult
    }

    factor(new.clusters)
}
# > factor(rep(c(2, 3, 1), c(5, 4, 3)))
#  [1] 2 2 2 2 2 3 3 3 3 1 1 1
# Levels: 1 2 3
# > scuttle:::.limit_cluster_size(factor(rep(c(2, 3, 1), c(5, 4, 3))), 10)
#  [1] 1 1 1 1 1 2 2 2 2 3 3 3
# Levels: 1 2 3
# > .limit_cluster_size.fix(factor(rep(c(2, 3, 1), c(5, 4, 3))), 10)
#  [1] 2 2 2 2 2 3 3 3 3 1 1 1
# Levels: 1 2 3


### Ok, continuing with the previous code at
## https://github.com/LTLA/scuttle/blob/3cb28efbf237ecb9b7a1973ab7e8957371260a32/R/pooledSizeFactors.R#L214

if (ncells!=length(clusters)) {
    stop("'ncol(x)' is not equal to 'length(clusters)'")
}
indices <- split(seq_along(clusters), clusters)

if (length(indices)==0L || any(lengths(indices)==0L)) {
    stop("zero cells in one of the clusters")
}

# Addigional sanity checks on various parameters.
if (!is.null(scaling) && length(scaling)!=ncol(x)) {
    stop("'length(scaling)' should be equal to 'ncol(x)'")
}

min.mean <- scuttle:::.guessMinMean(x, min.mean=min.mean, BPPARAM=BPPARAM)
# > min.mean
# [1] 0.1

sizes <- sort(as.integer(sizes))
if (anyDuplicated(sizes)) {
    stop("'sizes' are not unique")
}
# > sizes
#  [1]  21  26  31  36  41  46  51  56  61  66  71  76  81  86  91  96 101

# Fragmenting the matrices (and also scaling).
frag.x <- frag.scale <- vector("list", length(indices))
for (i in seq_along(indices)) {
    idx <- indices[[i]]
    if (length(indices) > 1L || !identical(idx, seq_along(idx))) {
        current <- x[,idx,drop=FALSE]
    } else {
        current <- x
    }
    if (!is.null(subset.row)) {
        current <- current[subset.row,,drop=FALSE]
    }
    frag.x[[i]] <- current
    frag.scale[i] <- list(scaling[idx]) # handle NULLs properly.
}

## Note: this is the function that takes like 8-10 min to run
# Computing normalization factors within each cluster.
all.norm <- BiocParallel::bpmapply(FUN=scuttle:::.per_cluster_normalize, x=frag.x, scaling=frag.scale,
    MoreArgs=list(sizes=sizes, min.mean=min.mean, positive=positive),
    BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)

clust.nf <- lapply(all.norm, "[[", i="final.nf")
clust.profile <- lapply(all.norm, "[[", i="ave.cell")

# Adjusting size factors between clusters.
if (is.null(ref.clust)) { ## yes, this is NULL
    non.zeroes <- vapply(clust.profile, FUN=function(x) sum(x>0), FUN.VALUE=0L)
    ref.clust <- which.max(non.zeroes)
    # > ref.clust
    # [1] 79
}
## This is the function that returns the warning (previously the error)
## as noted at https://github.com/LTLA/scuttle/issues/7
rescaling.factors <- scuttle:::.rescale_clusters(clust.profile, ref.col=ref.clust, min.mean=min.mean)
# Warning message:
# In scuttle:::.rescale_clusters(clust.profile, ref.col = ref.clust,  :
#   inter-cluster rescaling factor for cluster 64 is not strictly positive,
# reverting to the ratio of average library sizes

## Let's make the objects required for this function and dig into it
## https://github.com/LTLA/scuttle/blob/3cb28efbf237ecb9b7a1973ab7e8957371260a32/R/pooledSizeFactors.R#L408
mean.prof <- clust.profile
ref.col <- ref.clust
# min.mean already exists

## We continue with
## https://github.com/LTLA/scuttle/blob/3cb28efbf237ecb9b7a1973ab7e8957371260a32/R/pooledSizeFactors.R#L412-L420
## note that `ref.col` is not a character (it's an integer)
# > ref.col
# [1] 79
nclusters <- length(mean.prof)
rescaling <- numeric(nclusters)

## Let's set clust to 64 since that's the one currently prompting the warning
clust <- 64
## Next we can continue with
## https://github.com/LTLA/scuttle/blob/3cb28efbf237ecb9b7a1973ab7e8957371260a32/R/pooledSizeFactors.R#L422-L437
ref.prof <- mean.prof[[ref.col]]
cur.prof <- mean.prof[[clust]]
cur.libsize <- sum(cur.prof)
ref.libsize <- sum(ref.prof)
# > cur.libsize
# [1] 13.62893
# > ref.libsize
# [1] 4919.067
to.use <- (cur.prof/cur.libsize + ref.prof/ref.libsize)/2 * (cur.libsize + ref.libsize)/2 >= min.mean
# > summary(to.use)
#    Mode   FALSE    TRUE
# logical   23848    3165
if (!all(to.use)) {
    cur.prof <- cur.prof[to.use]
    ref.prof <- ref.prof[to.use]
}

# Adjusting for systematic differences between clusters.
rescale.sf <- median(cur.prof/ref.prof, na.rm=TRUE)
# > rescale.sf
# [1] 0

## Since it's 0, that triggers the warning (formerly error)
## https://github.com/LTLA/scuttle/blob/3cb28efbf237ecb9b7a1973ab7e8957371260a32/R/pooledSizeFactors.R#L412-L442
## See https://github.com/LTLA/scuttle/commit/0ed602b33c839b9cad770d5871b7436ebe993caf
## for changes
if (!is.finite(rescale.sf) || rescale.sf <= 0) {
    warning(paste(strwrap(paste0("inter-cluster rescaling factor for cluster ", clust,
        " is not strictly positive, reverting to the ratio of average library sizes")), collapse="\n"))
    rescale.sf <- sum(cur.prof)/sum(ref.prof)
    # > rescale.sf
    # [1] 0.003733694
    # > sum(cur.prof)
    # [1] 13.62893
    # > sum(ref.prof)
    # [1] 3650.254
}

# > summary(cur.prof)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000000 0.000000 0.000000 0.004306 0.004082 0.552340
# > summary(ref.prof)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.4312  0.6260  1.1533  1.1144 57.0261
# > length(cur.prof)
# [1] 3165
# > length(ref.prof)
# [1] 3165
# > addmargins(table(cur.prof == 0))
#
# FALSE  TRUE   Sum
#  1431  1734  3165
#  > addmargins(table(ref.prof == 0))
#
# FALSE  TRUE   Sum
#  3159     6  3165
save(all.norm, clusters, cur.prof, ref.prof, min.mean, file = here::here("rdata", "inspect_scuttle_issue_7.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R Under development (unstable) (2021-02-17 r80017)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-02-17
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version  date       lib source
#  AnnotationDbi            1.53.1   2021-02-04 [2] Bioconductor
#  AnnotationHub            2.23.2   2021-02-05 [2] Bioconductor
#  assertthat               0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  attempt                  0.3.1    2020-05-03 [1] CRAN (R 4.1.0)
#  beachmat                 2.7.6    2021-01-15 [2] Bioconductor
#  beeswarm                 0.2.3    2016-04-25 [1] CRAN (R 4.1.0)
#  benchmarkme              1.0.5    2021-02-09 [1] CRAN (R 4.1.0)
#  benchmarkmeData          1.0.4    2020-04-23 [1] CRAN (R 4.1.0)
#  Biobase                * 2.51.0   2020-10-27 [2] Bioconductor
#  BiocFileCache            1.15.1   2020-11-09 [2] Bioconductor
#  BiocGenerics           * 0.37.1   2021-02-04 [2] Bioconductor
#  BiocManager              1.30.10  2019-11-16 [2] CRAN (R 4.1.0)
#  BiocNeighbors            1.9.4    2020-12-17 [1] Bioconductor
#  BiocParallel             1.25.4   2021-02-04 [2] Bioconductor
#  BiocSingular             1.7.2    2021-01-23 [1] Bioconductor
#  BiocVersion              3.13.1   2020-10-27 [2] Bioconductor
#  Biostrings               2.59.2   2020-12-18 [2] Bioconductor
#  bit                      4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
#  bit64                    4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
#  bitops                   1.0-6    2013-08-17 [2] CRAN (R 4.1.0)
#  blob                     1.2.1    2020-01-20 [2] CRAN (R 4.1.0)
#  bluster                  1.1.5    2021-01-14 [1] Bioconductor
#  cachem                   1.0.4    2021-02-13 [2] CRAN (R 4.1.0)
#  cli                      2.3.0    2021-01-31 [2] CRAN (R 4.1.0)
#  codetools                0.2-18   2020-11-04 [3] CRAN (R 4.1.0)
#  colorout                 1.2-2    2021-02-11 [1] Github (jalvesaq/colorout@726d681)
#  colorspace               2.0-0    2020-11-11 [2] CRAN (R 4.1.0)
#  config                   0.3.1    2020-12-17 [1] CRAN (R 4.1.0)
#  cowplot                  1.1.1    2020-12-30 [1] CRAN (R 4.1.0)
#  crayon                   1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
#  curl                     4.3      2019-12-02 [2] CRAN (R 4.1.0)
#  data.table               1.13.6   2020-12-30 [2] CRAN (R 4.1.0)
#  DBI                      1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
#  dbplyr                   2.1.0    2021-02-03 [2] CRAN (R 4.1.0)
#  DelayedArray             0.17.7   2020-12-26 [2] Bioconductor
#  DelayedMatrixStats       1.13.5   2021-02-04 [2] Bioconductor
#  desc                     1.2.0    2018-05-01 [2] CRAN (R 4.1.0)
#  digest                   0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
#  dockerfiler              0.1.3    2019-03-19 [1] CRAN (R 4.1.0)
#  doParallel               1.0.16   2020-10-16 [2] CRAN (R 4.1.0)
#  dotCall64                1.0-1    2021-02-11 [1] CRAN (R 4.1.0)
#  dplyr                    1.0.4    2021-02-02 [2] CRAN (R 4.1.0)
#  dqrng                    0.2.1    2019-05-17 [1] CRAN (R 4.1.0)
#  DropletUtils             1.11.10  2021-02-04 [1] Bioconductor
#  DT                       0.17     2021-01-06 [2] CRAN (R 4.1.0)
#  edgeR                    3.33.1   2021-01-11 [2] Bioconductor
#  ellipsis                 0.3.1    2020-05-15 [2] CRAN (R 4.1.0)
#  ExperimentHub            1.17.1   2021-02-08 [2] Bioconductor
#  farver                   2.0.3    2020-01-16 [2] CRAN (R 4.1.0)
#  fastmap                  1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
#  fields                   11.6     2020-10-09 [2] CRAN (R 4.1.0)
#  filelock                 1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
#  foreach                  1.5.1    2020-10-15 [2] CRAN (R 4.1.0)
#  fs                       1.5.0    2020-07-31 [2] CRAN (R 4.1.0)
#  generics                 0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
#  GenomeInfoDb           * 1.27.6   2021-02-04 [2] Bioconductor
#  GenomeInfoDbData         1.2.4    2020-11-03 [2] Bioconductor
#  GenomicRanges          * 1.43.3   2021-01-14 [2] Bioconductor
#  ggbeeswarm               0.6.0    2017-08-07 [1] CRAN (R 4.1.0)
#  ggplot2                  3.3.3    2020-12-30 [2] CRAN (R 4.1.0)
#  gh                       1.2.0    2020-11-27 [2] CRAN (R 4.1.0)
#  glue                     1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
#  golem                    0.2.1    2020-03-05 [1] CRAN (R 4.1.0)
#  gridExtra                2.3      2017-09-09 [2] CRAN (R 4.1.0)
#  gtable                   0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  HDF5Array                1.19.4   2021-02-14 [2] Bioconductor
#  here                   * 1.0.1    2020-12-13 [1] CRAN (R 4.1.0)
#  htmltools                0.5.1.1  2021-01-22 [2] CRAN (R 4.1.0)
#  htmlwidgets              1.5.3    2020-12-10 [2] CRAN (R 4.1.0)
#  httpuv                   1.5.5    2021-01-13 [2] CRAN (R 4.1.0)
#  httr                     1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
#  igraph                   1.2.6    2020-10-06 [2] CRAN (R 4.1.0)
#  interactiveDisplayBase   1.29.0   2020-10-27 [2] Bioconductor
#  IRanges                * 2.25.6   2020-12-18 [2] Bioconductor
#  irlba                    2.3.3    2019-02-05 [2] CRAN (R 4.1.0)
#  iterators                1.0.13   2020-10-15 [2] CRAN (R 4.1.0)
#  jsonlite                 1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
#  KEGGREST                 1.31.1   2020-11-23 [2] Bioconductor
#  knitr                    1.31     2021-01-27 [2] CRAN (R 4.1.0)
#  labeling                 0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
#  later                    1.1.0.1  2020-06-05 [2] CRAN (R 4.1.0)
#  lattice                  0.20-41  2020-04-02 [3] CRAN (R 4.1.0)
#  lazyeval                 0.2.2    2019-03-15 [2] CRAN (R 4.1.0)
#  lifecycle                1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
#  limma                    3.47.7   2021-02-15 [2] Bioconductor
#  locfit                   1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
#  magick                   2.6.0    2021-01-13 [2] CRAN (R 4.1.0)
#  magrittr                 2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
#  maps                     3.3.0    2018-04-03 [2] CRAN (R 4.1.0)
#  Matrix                   1.3-2    2021-01-06 [3] CRAN (R 4.1.0)
#  MatrixGenerics         * 1.3.1    2021-02-01 [2] Bioconductor
#  matrixStats            * 0.58.0   2021-01-29 [2] CRAN (R 4.1.0)
#  memoise                  2.0.0    2021-01-26 [2] CRAN (R 4.1.0)
#  metapod                  0.99.5   2020-12-14 [1] Bioconductor
#  mime                     0.10     2021-02-13 [2] CRAN (R 4.1.0)
#  munsell                  0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                   1.4.7    2020-11-20 [2] CRAN (R 4.1.0)
#  pkgconfig                2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  pkgload                  1.1.0    2020-05-29 [2] CRAN (R 4.1.0)
#  plotly                   4.9.3    2021-01-10 [2] CRAN (R 4.1.0)
#  png                      0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
#  Polychrome               1.2.6    2020-11-11 [1] CRAN (R 4.1.0)
#  promises                 1.2.0.1  2021-02-11 [1] CRAN (R 4.1.0)
#  purrr                    0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R.methodsS3              1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
#  R.oo                     1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
#  R.utils                  2.10.1   2020-08-26 [2] CRAN (R 4.1.0)
#  R6                       2.5.0    2020-10-28 [2] CRAN (R 4.1.0)
#  rappdirs                 0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
#  RColorBrewer             1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
#  Rcpp                     1.0.6    2021-01-15 [2] CRAN (R 4.1.0)
#  RCurl                    1.98-1.2 2020-04-18 [2] CRAN (R 4.1.0)
#  remotes                  2.2.0    2020-07-21 [2] CRAN (R 4.1.0)
#  rhdf5                    2.35.0   2020-10-27 [2] Bioconductor
#  rhdf5filters             1.3.3    2020-12-07 [2] Bioconductor
#  Rhdf5lib                 1.13.0   2020-10-27 [2] Bioconductor
#  rjson                    0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
#  rlang                    0.4.10   2020-12-30 [2] CRAN (R 4.1.0)
#  roxygen2                 7.1.1    2020-06-27 [2] CRAN (R 4.1.0)
#  rprojroot                2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  RSQLite                  2.2.3    2021-01-24 [2] CRAN (R 4.1.0)
#  rstudioapi               0.13     2020-11-12 [2] CRAN (R 4.1.0)
#  rsvd                     1.0.3    2020-02-17 [1] CRAN (R 4.1.0)
#  S4Vectors              * 0.29.7   2021-02-04 [2] Bioconductor
#  ScaledMatrix             0.99.2   2021-01-14 [1] Bioconductor
#  scales                   1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  scater                   1.19.9   2021-02-01 [1] Bioconductor
#  scatterplot3d            0.3-41   2018-03-14 [1] CRAN (R 4.1.0)
#  scran                    1.19.13  2021-02-13 [1] Bioconductor
#  scuttle                  1.1.15   2021-02-14 [1] Bioconductor
#  sessioninfo            * 1.1.1    2018-11-05 [2] CRAN (R 4.1.0)
#  shiny                    1.6.0    2021-01-25 [2] CRAN (R 4.1.0)
#  shinyWidgets             0.5.7    2021-02-03 [1] CRAN (R 4.1.0)
#  SingleCellExperiment   * 1.13.10  2021-02-10 [1] Bioconductor
#  spam                     2.6-0    2020-12-14 [2] CRAN (R 4.1.0)
#  sparseMatrixStats        1.3.6    2021-02-04 [2] Bioconductor
#  SpatialExperiment      * 1.1.432  2021-02-17 [1] Github (drighelli/SpatialExperiment@8407ee8)
#  spatialLIBD            * 1.3.5    2021-02-17 [1] Github (LieberInstitute/spatialLIBD@3fbc875)
#  statmod                  1.4.35   2020-10-19 [2] CRAN (R 4.1.0)
#  stringi                  1.5.3    2020-09-09 [2] CRAN (R 4.1.0)
#  stringr                  1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
#  SummarizedExperiment   * 1.21.1   2020-12-12 [2] Bioconductor
#  testthat                 3.0.2    2021-02-14 [2] CRAN (R 4.1.0)
#  tibble                   3.0.6    2021-01-29 [2] CRAN (R 4.1.0)
#  tidyr                    1.1.2    2020-08-27 [2] CRAN (R 4.1.0)
#  tidyselect               1.1.0    2020-05-11 [2] CRAN (R 4.1.0)
#  usethis                  2.0.1    2021-02-10 [2] CRAN (R 4.1.0)
#  vctrs                    0.3.6    2020-12-17 [2] CRAN (R 4.1.0)
#  vipor                    0.4.5    2017-03-22 [1] CRAN (R 4.1.0)
#  viridis                  0.5.1    2018-03-29 [2] CRAN (R 4.1.0)
#  viridisLite              0.3.0    2018-02-01 [2] CRAN (R 4.1.0)
#  withr                    2.4.1    2021-01-26 [2] CRAN (R 4.1.0)
#  xfun                     0.21     2021-02-10 [2] CRAN (R 4.1.0)
#  xml2                     1.3.2    2020-04-23 [2] CRAN (R 4.1.0)
#  xtable                   1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
#  XVector                  0.31.1   2020-12-12 [2] Bioconductor
#  yaml                     2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
#  zlibbioc                 1.37.0   2020-10-27 [2] Bioconductor
#
# [1] /users/lcollado/R/devel
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/library
