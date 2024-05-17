spatialDLPFC
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://zenodo.org/badge/314001778.svg)](https://zenodo.org/badge/latestdoi/314001778)

## Overview

<img src="http://research.libd.org/spatialDLPFC/img/Br6255_ant_Sp09.png" align="left" width="300px" />

Welcome to the `spatialDLPFC` project! This project involves 3 data
types as well as several interactive websites, all of which you are
publicly accessible for you to browse and download.

In this project we studied spatially resolved and single nucleus
transcriptomics data from the dorsolateral prefrontal cortex (DLPFC)
from postmortem human brain samples. From 10 neurotypical controls we
generated spatially-resolved transcriptomics data using using [10x
Genomics
**Visium**](https://www.10xgenomics.com/products/spatial-gene-expression)
across the anterior, middle, and posterior DLPFC (n = 30). We also
generated single nucleus RNA-seq (**snRNA-seq**) data using [10x
Genomics
**Chromium**](https://www.10xgenomics.com/products/single-cell-gene-expression)
from 19 of these tissue blocks. We further generated data from 4
adjacent tissue slices with [10x Genomics **Visium Spatial
Proteogenomics**
(SPG)](https://www.10xgenomics.com/products/spatial-proteogenomics),
that can be used to benchmark spot deconvolution algorithms. This work
is being was performed by the [Keri
Martinowich](https://www.libd.org/team/keri-martinowich-phd/), [Leonardo
Collado-Torres](http://lcolladotor.github.io/), and [Kristen
Maynard](https://www.libd.org/team/kristen-maynard-phd/) teams at the
[Lieber Institute for Brain Development](libd.org) as well as [Stephanie
Hicks](https://www.stephaniehicks.com/)’s group from [JHBSPH’s
Biostatistics
Department](https://publichealth.jhu.edu/departments/biostatistics).

This project involves the GitHub repositories
[LieberInstitute/spatialDLPFC](https://github.com/LieberInstitute/spatialDLPFC)
and
[LieberInstitute/DLPFC_snRNAseq](https://github.com/LieberInstitute/DLPFC_snRNAseq).

If you tweet about this website, the data or the R package please use
the <code>\#spatialDLPFC</code> hashtag. You can find previous tweets
that way as shown
<a href="https://twitter.com/search?q=%23spatialDLPFC&src=typed_query">here</a>.

Thank you for your interest in our work!

## Study Design

<img src="http://research.libd.org/spatialDLPFC/img/study_overview.png" width="1000px" align="left" />

**Study design to generate paired single nucleus RNA-sequencing
(snRNA-seq) and spatially-resolved transcriptomic data across DLPFC**.
(**A**) DLPFC tissue blocks were dissected across the rostral-caudal
axis from 10 adult neurotypical control postmortem human brains,
including anterior (Ant), middle (Mid), and posterior (Post) positions
(n=3 blocks per donor, n=30 blocks total). The same tissue blocks were
used for snRNA-seq (10x Genomics 3’ gene expression assay, n=1-2 blocks
per donor, n=19 samples) and spatial transcriptomics (10x Genomics
Visium spatial gene expression assay, n=3 blocks per donor, n=30
samples). (**B**) Paired snRNA-seq and Visium data were used to identify
data-driven spatial domains (SpDs) and cell types, perform spot
deconvolution, conduct cell-cell communication analyses, and spatially
register companion PsychENCODE snRNA-seq DLPFC data. (**C**)
*t*-distributed stochastic neighbor embedding (t-SNE) summarizing layer
resolution cell types identified by snRNA-seq. (**D**) Tissue block
orientation and morphology was confirmed by hematoxylin and eosin (H&E)
staining and single molecule fluorescent in situ hybridization (smFISH)
with RNAscope (*SLC17A7* marking excitatory neurons in pink, *MBP*
marking white matter (WM) in green, *RELN* marking layer (L)1 in yellow,
and *NR4A2* marking L6 in orange). Scale bar is 2mm. Spotplots depicting
log transformed normalized expression (logcounts) of *SNAP25*, *MBP*,
and *PCP4* in the Visium data confirm the presence of gray matter, WM,
and cortical layers, respectively. (**E**) Schematic of unsupervised SpD
identification and registration using *BayesSpace* SpDs at *k*=7.
Enrichment *t*-statistics computed on *BayesSpace* SpDs were correlated
with manual histological layer annotations from ([Maynard,
Collado-Torres et al., 2021, *Nat
Neuro*](https://doi.org/10.1038/s41593-020-00787-0)) to map SpDs to
known histological layers. The heatmap of correlation values summarizes
the relationship between BayesSpace SpDs and classic histological
layers. Higher confidence annotations (⍴ \> 0.25, merge ratio = 0.1) are
marked with an “X”.

## Interactive Websites

All of these interactive websites are powered by open source software,
namely:

- 🔭 [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w)
- 🔍 [`samui`](http://dx.doi.org/10.1017/S2633903X2300017X)
- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)

We provide the following interactive websites, organized by dataset with
software labeled by emojis:

- Visium (n = 30)
  - 🔭
    [spatialDLPFC_Visium_Sp09](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp09):
    [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) website
    showing the spatially-resolved Visium data (n = 30) with statistical
    results comparing the Sp09 domains.
  - 🔭
    [spatialDLPFC_Visium_Sp16](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp16):
    similar but with the Sp16 domains.
  - 🔭
    [spatialDLPFC_Visium_Sp09_position](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp09_position):
    similar to *spatialDLPFC_Visium_Sp09* but with statistical results
    across the *position* (anterior, middle, posterior) adjusting for
    the Sp09 domains.
  - 🔭
    [spatialDLPFC_Visium_Sp09_position_noWM](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp09_position_noWM):
    similar to *spatialDLPFC_Visium_Sp09_position* but after dropping
    the `SP28D06`, `SP28D16`, `SP28D17`, `SP28D20` and `SP28D28` spots
    which correspond to white matter (hence the `noWM` acronym).
  - 👀
    [spatialDLPFC_Visium_Sp09_pseudobulk](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp09_pseudobulk):
    [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1) website
    showing the pseudo-bulked Sp09 domains spatial data.
  - 👀
    [spatialDLPFC_Visium_Sp16_pseudobulk](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp16_pseudobulk):
    similar to *spatialDLPFC_Visium_Sp09_pseudobulk* but with the Sp16
    domains data.
  - 👀
    [spatialDLPFC_Visium_Sp28_pseudobulk](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp28_pseudobulk):
    similar to *spatialDLPFC_Visium_Sp09_pseudobulk* but with the Sp28
    domains data.
  - 🔍 [spatialDLPFC Visium on
    Samui](https://samuibrowser.com/from?url=data.samuibrowser.com/spatialDLPFC/&s=Br2720_ant&s=Br2720_mid&s=Br2720_post&s=Br2743_ant&s=Br2743_mid&s=Br2743_post&s=Br3942_ant&s=Br3942_mid&s=Br3942_post&s=Br6423_ant&s=Br6423_mid&s=Br6423_post&s=Br6432_ant&s=Br6432_mid&s=Br6432_post&s=Br6471_ant&s=Br6471_mid&s=Br6471_post&s=Br6522_ant&s=Br6522_mid&s=Br6522_post&s=Br8325_ant&s=Br8325_mid&s=Br8325_post&s=Br8492_ant&s=Br8492_mid&s=Br8492_post&s=Br8667_ant&s=Br8667_mid&s=Br8667_post):
    [`samui`](https://github.com/chaichontat/samui) website that allows
    to zoom in at the spot or cell level.
- snRNA-seq (n = 19)
  - 👀
    [spatialDLPFC_snRNA-seq](https://libd.shinyapps.io/spatialDLPFC_snRNA-seq):
    [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1) website
    showing the n = 19 snRNA-seq samples at single nucleus resolution.
- Visium SPG (n = 4)
  - 🔭
    [spatialDLPFC_Visium_SPG](https://libd.shinyapps.io/spatialDLPFC_Visium_SPG):
    [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) website
    showing the spatially-resolved data Visium SPG (n = 4).
  - 🔍 [spatialDLPFC Visium SPG on
    Samui](https://samuibrowser.com/from?url=data.samuibrowser.com/VisiumIF/&s=Br2720_Ant_IF&s=Br6432_Ant_IF&s=Br6522_Ant_IF&s=Br8667_Post_IF):
    [`samui`](https://github.com/chaichontat/samui) website that allows
    to zoom in at the spot or cell level.

### Local `spatialLIBD` apps

If you are interested in running the
[`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) applications
locally, you can do so thanks to the
[`spatialLIBD::run_app()`](http://research.libd.org/spatialLIBD/reference/run_app.html),
which you can also use with your own data as shown in our [vignette for
publicly available datasets provided by 10x
Genomics](http://bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/TenX_data_download.html).

``` r
## Run this web application locally with:
spatialLIBD::run_app()

## You will have more control about the length of the session and memory usage.
## See http://research.libd.org/spatialLIBD/reference/run_app.html#examples
## for the full R code to run https://libd.shinyapps.io/spatialDLPFC_Visium_Sp09
## locally. See also:
## * https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/deploy_app_k09
## * https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/deploy_app_k09_position
## * https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/deploy_app_k09_position_noWM
## * https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/deploy_app_k16
## * https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/analysis_IF/03_spatialLIBD_app

## You could also use spatialLIBD::run_app() to visualize your
## own data given some requirements described
## in detail in the package vignette documentation
## at http://research.libd.org/spatialLIBD/.
```

## Contact

We value public questions, as they allow other users to learn from the
answers. If you have any questions, please ask them at
[LieberInstitute/spatialDLPFC/issues](https://github.com/LieberInstitute/spatialDLPFC/issues)
and refrain from emailing us. Thank you again for your interest in our
work!

## Citing our work

Please cite this [manuscript](https://doi.org/10.1126/science.adh1938)
if you use data from this project.

> A data-driven single-cell and spatial transcriptomic map of the human
> prefrontal cortex Louise A. Huuki-Myers, Abby Spangler, Nicholas J.
> Eagles, Kelsey D. Montgomery, Sang Ho Kwon, Boyi Guo, Melissa
> Grant-Peters, Heena R. Divecha, Madhavi Tippani, Chaichontat
> Sriworarat, Annie B. Nguyen, Prashanthi Ravichandran, Matthew N. Tran,
> Arta Seyedian, PsychENCODE Consortium, Thomas M. Hyde, Joel E.
> Kleinman, Alexis Battle, Stephanie C. Page, Mina Ryten, Stephanie C.
> Hicks, Keri Martinowich, Leonardo Collado-Torres, Kristen R. Maynard
> *Science* 384, eadh1938 (2024).; doi:
> <https://doi.org/10.1126/science.adh1938>

Below is the citation in [`BibTeX`](http://www.bibtex.org/) format.

    @article {Huuki-Myers2024.eadh1938,
        author = {Huuki-Myers, Louise A. and Spangler, Abby and Eagles, Nicholas J. and Montgomery, Kelsey D. and Kwon, Sang Ho and Guo, Boyi and Grant-Peters, Melissa and Divecha, Heena R. and Tippani, Madhavi and Sriworarat, Chaichontat and Nguyen, Annie B. and Ravichandran, Prashanthi and Tran, Matthew N. and Seyedian, Arta and , and Hyde, Thomas M. and Kleinman, Joel E. and Battle, Alexis and Page, Stephanie C. and Ryten, Mina and Hicks, Stephanie C. and Martinowich, Keri and Collado-Torres, Leonardo and Maynard, Kristen R.},
        title = {A data-driven single-cell and spatial transcriptomic map of the human prefrontal cortex},
        year = {2024},
        doi = {10.1126/science.adh1938},
        publisher = {American Association for the Advancement of Science (AAAS)},
        URL = {https://doi.org/10.1126/science.adh1938},
        journal = {Science}
    }

### Cite `spatialLIBD`

Below is the citation output from using `citation('spatialLIBD')` in R.
Please run this yourself to check for any updates on how to cite
**spatialLIBD**.

``` r
print(citation("spatialLIBD")[1], bibtex = TRUE)
#> Pardo B, Spangler A, Weber LM, Hicks SC, Jaffe AE, Martinowich K,
#> Maynard KR, Collado-Torres L (2022). "spatialLIBD: an R/Bioconductor
#> package to visualize spatially-resolved transcriptomics data." _BMC
#> Genomics_. doi:10.1186/s12864-022-08601-w
#> <https://doi.org/10.1186/s12864-022-08601-w>,
#> <https://doi.org/10.1186/s12864-022-08601-w>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {spatialLIBD: an R/Bioconductor package to visualize spatially-resolved transcriptomics data},
#>     author = {Brenda Pardo and Abby Spangler and Lukas M. Weber and Stephanie C. Hicks and Andrew E. Jaffe and Keri Martinowich and Kristen R. Maynard and Leonardo Collado-Torres},
#>     year = {2022},
#>     journal = {BMC Genomics},
#>     doi = {10.1186/s12864-022-08601-w},
#>     url = {https://doi.org/10.1186/s12864-022-08601-w},
#>   }
```

Please note that the `spatialLIBD` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing the package.

### Cite `Samui`

To cite [`samui`](https://github.com/chaichontat/samui) please use:

> Performant web-based interactive visualization tool for
> spatially-resolved transcriptomics experiments Chaichontat Sriworarat,
> Annie Nguyen, Nicholas J. Eagles, Leonardo Collado-Torres, Keri
> Martinowich, Kristen R. Maynard, Stephanie C. Hicks Biological
> Imaging; doi: <https://doi.org/10.1017/S2633903X2300017X>

Below is the citation in [`BibTeX`](http://www.bibtex.org/) format.

    @article{sriworarat_performant_2023,
        title = {Performant web-based interactive visualization tool for spatially-resolved transcriptomics experiments},
        volume = {3},
        issn = {2633-903X},
        url = {https://www.cambridge.org/core/journals/biological-imaging/article/performant-webbased-interactive-visualization-tool-for-spatiallyresolved-transcriptomics-experiments/B66303984D10B9E5A23D3656CB8537C0},
        doi = {10.1017/S2633903X2300017X},
        language = {en},
        urldate = {2024-04-19},
        journal = {Biological Imaging},
        author = {Sriworarat, Chaichontat and Nguyen, Annie and Eagles, Nicholas J. and Collado-Torres, Leonardo and Martinowich, Keri and Maynard, Kristen R. and Hicks, Stephanie C.},
        month = jan,
        year = {2023},
        keywords = {georeferencing, interactive image viewer, multi-dimensional image, single-cell transcriptomics, spatially resolved transcriptomics, web-based browser},
        pages = {e15}
    }

### Cite `VistoSeg`

To cite [`VistoSeg`](http://research.libd.org/VistoSeg/) please use:

> VistoSeg: {Processing utilities for high-resolution images for
> spatially resolved transcriptomics data. Madhavi Tippani, Heena R.
> Divecha, Joseph L. Catallini II, Sang Ho Kwon, Lukas M. Weber, Abby
> Spangler, Andrew E. Jaffe, Thomas M. Hyde, Joel E. Kleinman, Stephanie
> C. Hicks, Keri Martinowich, Leonardo Collado-Torres, Stephanie C.
> Page, Kristen R. Maynard Biological Imaging ; doi:
> <https://doi.org/10.1017/S2633903X23000235>

Below is the citation in [`BibTeX`](http://www.bibtex.org/) format.

    @article{tippani_vistoseg_2023,
        title = {{VistoSeg}: {Processing} utilities for high-resolution images for spatially resolved transcriptomics data},
        volume = {3},
        issn = {2633-903X},
        shorttitle = {{VistoSeg}},
        url = {https://www.cambridge.org/core/journals/biological-imaging/article/vistoseg-processing-utilities-for-highresolution-images-for-spatially-resolved-transcriptomics-data/990CBC4AC069F5EDC62316919398404B},
        doi = {10.1017/S2633903X23000235},
        language = {en},
        urldate = {2024-04-19},
        journal = {Biological Imaging},
        author = {Tippani, Madhavi and Divecha, Heena R. and Catallini, Joseph L. and Kwon, Sang H. and Weber, Lukas M. and Spangler, Abby and Jaffe, Andrew E. and Hyde, Thomas M. and Kleinman, Joel E. and Hicks, Stephanie C. and Martinowich, Keri and Collado-Torres, Leonardo and Page, Stephanie C. and Maynard, Kristen R.},
        month = jan,
        year = {2023},
        keywords = {hematoxylin and eosin, immunofluorescence, MATLAB, segmentation, spatially resolved transcriptomics, Visium, Visium-Spatial Proteogenomics},
        pages = {e23}
    }

## Data Access

We highly value open data sharing and believe that doing so accelerates
science, as was the case between our
[`HumanPilot`](https://doi.org/10.1038/s41593-020-00787-0) and the
external [`BayesSpace`](https://doi.org/10.1038/s41587-021-00935-2)
projects, documented [on this
slide](https://speakerdeck.com/lcolladotor/hca-la-2022?slide=18).

### Processed Data

[`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) also allows
you to access the data from this project as ready to use R objects. That
is, a:

- [`SpatialExperiment`](https://doi.org/10.1093/bioinformatics/btac299)
  object for the Visium samples (n = 30)
- [`SpatialExperiment`](https://doi.org/10.1093/bioinformatics/btac299)
  object for the Visium SPG samples (n = 4)
- [`SingleCellExperiment`](https://www.nature.com/articles/s41592-019-0654-x)
  object for the snRNA-seq samples (n = 19)

You can use the
[`zellkonverter`](https://bioconductor.org/packages/zellkonverter/)
Bioconductor package to convert any of them into Python
[`AnnData`](https://anndata.readthedocs.io/en/latest/) objects. If you
browse our code, you can find examples of such conversions.

If you are unfamiliar with these tools, you might want to check the
[LIBD rstats club](http://research.libd.org/rstatsclub/#.Y4hWlOzMJUM)
(check and search keywords on the
[schedule](https://docs.google.com/spreadsheets/d/1is8dZSd0FZ9Qi1Zvq1uRhm-P1McnJRd_zxdAfCRoMfA/edit?usp=sharing))
videos and resources.

#### Installing spatialLIBD

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `spatialLIBD` from
[Bioconductor](http://bioconductor.org/) with the following code:

``` r
## Install BiocManager in order to install Bioconductor packages properly
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

## Check that you have a valid R/Bioconductor installation
BiocManager::valid()

## Now install spatialLIBD from Bioconductor
## (this version has been tested on macOS, winOS, linux)
BiocManager::install("spatialLIBD")

## If you need the development version from GitHub you can use the following:
# BiocManager::install("LieberInstitute/spatialLIBD")
## Note that this version might include changes that have not been tested
## properly on all operating systems.
```

### R objects

Using `spatialLIBD` you can access the spatialDLPFC transcriptomics data
from the 10x Genomics Visium platform. For example, this is the code you
can use to access the spatially-resolved data. For more details, check
the help file for `fetch_data()`.

``` r
## Check that you have a recent version of spatialLIBD installed
stopifnot(packageVersion("spatialLIBD") >= "1.11.6")

## Download the spot-level data
spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")

## This is a SpatialExperiment object
spe
#> class: SpatialExperiment 
#> dim: 28916 113927 
#> metadata(1): BayesSpace.data
#> assays(2): counts logcounts
#> rownames(28916): ENSG00000243485 ENSG00000238009 ... ENSG00000278817 ENSG00000277196
#> rowData names(7): source type ... gene_type gene_search
#> colnames(113927): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ... TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
#> colData names(155): age array_col ... VistoSeg_proportion wrinkle_type
#> reducedDimNames(8): 10x_pca 10x_tsne ... HARMONY UMAP.HARMONY
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor

## Note the memory size
lobstr::obj_size(spe)
#> 6.97 GB

## Set the cluster colors
colors_BayesSpace <- Polychrome::palette36.colors(28)
names(colors_BayesSpace) <- seq_len(28)

## Remake the logo image with histology information
p09 <- spatialLIBD::vis_clus(
    spe = spe,
    clustervar = "BayesSpace_harmony_09",
    sampleid = "Br6522_ant",
    colors = colors_BayesSpace,
    ... = " spatialDLPFC Human Brain\nSp09 domains -- made with spatialLIBD"
)
p09
```

<a href="https://libd.shinyapps.io/spatialDLPFC_Visium_Sp09"><img src="http://research.libd.org/spatialDLPFC/img/Br6255_ant_Sp09.png" width="800px" align="center" /></a>

``` r
## Repeat but for Sp16
p16 <- spatialLIBD::vis_clus(
    spe = spe,
    clustervar = "BayesSpace_harmony_16",
    sampleid = "Br6522_ant",
    colors = colors_BayesSpace,
    ... = " spatialDLPFC Human Brain\nSp16 domains -- made with spatialLIBD"
)
p16
```

<a href="https://libd.shinyapps.io/spatialDLPFC_Visium_Sp16"><img src="http://research.libd.org/spatialDLPFC/img/Br6255_ant_Sp16.png" width="800px" align="center" /></a>

### Raw data

The source data described in this manuscript are available from the
[National Institute of Mental Health (NIMH) Data
Archive](https://nda.nih.gov/) under [NDA Study
2619](https://nda.nih.gov/study.html?id=2619); doi:
[10.15154/7893-6778](https://doi.org/10.15154/7893-6778).

You can also access all the raw data through
[Globus](http://research.libd.org/globus/) (`jhpce#spatialDLPFC` and
`jhpce#DLPFC_snRNAseq`). This includes all the input FASTQ files as well
as the outputs from tools such as
[`SpaceRanger`](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger)
or
[`CellRanger`](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).
The files are *mostly* organized following the
[LieberInstitute/template_project](https://github.com/LieberInstitute/template_project)
project structure.

## Internal

- JHPCE locations:
  - `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC`
  - `/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq`
- Slack channel:
  [`libd_dlpfc_spatial`](https://jhu-genomics.slack.com/archives/C01EA7VDJNT).

### Files: [`spatialDLPFC`](https://github.com/LieberInstitute/spatialDLPFC)

- `code`: R, python, and shell scripts for running various analyses.
  - `spot_deconvo`: cell-type deconvolution within Visium spots, enabled
    by tools like
    [`tangram`](https://doi.org/10.1038/s41592-021-01264-7),
    [`cell2location`](https://doi.org/10.1038/s41587-021-01139-4),
    [`cellpose`](https://cellpose.readthedocs.io/en/latest/), and
    [`SPOTlight`](https://doi.org/10.1093/nar/gkab043)
  - `spython`: older legacy testing scripts mostly replaced by
    `spot_deconvo`
- `plots`: plots generated by RMarkdown or R analysis scripts in `.pdf`
  or `.png` format
- `processed-data`
  - `images_spatialLIBD`: images used for running `SpaceRanger`
  - `NextSeq`: `SpaceRanger` output files
  - `rdata`: R objects
- `raw-data`
  - `FASTQ`: FASTQ files from `NextSeq` runs.
  - `FASTQ_renamed`: renamed symbolic links to the original FASTQs, with
    consistent nomenclature
  - `Images`: raw images from the scanner in `.tif` format and around 3
    GB per sample.
  - `images_raw_align_json`
  - `psychENCODE`: external data from PsychENCODE (doi:
    [10.7303/syn2787333](https://doi.org/10.7303/syn2787333)).
  - `sample_info`: spreadsheet with information about samples (sample
    ID, sample name, slide serial number, capture area ID)

This GitHub repository is organized along the [*R/Bioconductor-powered
Team Data Science* group
guidelines](https://lcolladotor.github.io/bioc_team_ds/organizing-your-work.html#.Yaf9fPHMIdk).
It aims to follow the
[LieberInstitute/template_project](https://github.com/LieberInstitute/template_project)
structure, though most of the `code/analysis` output is saved at
`processed-data/rdata/spe` directory unlike what’s specified in the
template structure. This is due to historical reasons.

### Files: [`DLPFC_snRNAseq`](https://github.com/LieberInstitute/DLPFC_snRNAseq)

- `code`: R scripts for running various analyses.
- `plots`: plots generated by RMarkdown or R analysis scripts in `.pdf`
  or `.png` format
- `processed-data`
  - `cellranger`: `CellRanger` output files
- `raw-data`
  - `FASTQ`: FASTQ files.
  - `sample_info`: spreadsheet with information about samples (sample
    ID, sample name)

This GitHub repository is organized along the [*R/Bioconductor-powered
Team Data Science* group
guidelines](https://lcolladotor.github.io/bioc_team_ds/organizing-your-work.html#.Yaf9fPHMIdk).
It aims to follow the
[LieberInstitute/template_project](https://github.com/LieberInstitute/template_project)
structure.

### Other related files

- Reference transcriptome from 10x Genomics:
  `/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/`
