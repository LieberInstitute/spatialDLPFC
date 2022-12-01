
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialDLPFC <img src="http://research.libd.org/spatialDLPFC/img/Br6255_ant_Sp09.png" align="right" width="300px" />

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
Chromium](https://www.10xgenomics.com/products/single-cell-gene-expression)
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
the <code>\#spatialLIBD</code> hashtag. You can find previous tweets
that way as shown
<a href="https://twitter.com/search?q=%23spatialLIBD&src=typed_query">here</a>.

Thank you for your interest in our work!

## Study Design

<img src="http://research.libd.org/spatialDLPFC/img/study_overview.png" width="600px" align="right" />

**Study design to generate paired single-nucleus RNA sequencing
(snRNA-seq) and spatially-resolved transcriptomic data across DLPFC**.
**A**. Tissue blocks were dissected across the rostral-caudal axis from
10 adult neurotypical control postmortem human brains of the DLPFC,
including anterior (Ant), middle (Mid, and posterior (Post) positions
(n=3 blocks per donor, n=30 blocks total). **B**. The same tissue blocks
were used for snRNA-seq (10x Genomics 3’ gene expression assay, n=1-2
blocks per donor, n=19 samples) and spatial transcriptomics (10x
Genomics Visium spatial gene expression assay, n=3 blocks per donor,
n=30 samples). **C**. Tissue block orientation and morphology was
confirmed by single molecule fluorescent in situ hybridization (smFISH)
for laminar marker genes with RNAscope (*SLC17A7* marking excitatory
neurons in pink, *MBP* marking white matter in green, *RELN* marking
layer 1 in yellow, and *NR4A2* marking layer 6 in orange) and
hematoxylin and eosin (H&E) staining. Spotplots depicting log
transformed normalized expression (logcounts) of *SNAP25*, *MBP*, and
*PCP4* in the Visium data confirm the presence of gray matter, white
matter, and cortical layers, respectively.

## Citing our work

Please cite this [manuscript](TODO) if you use data from this project.
Below is the citation in [`BibTeX`](http://www.bibtex.org/) format.

    TODO

Below is the citation output from using `citation('spatialLIBD')` in R.
Please run this yourself to check for any updates on how to cite
**spatialLIBD**.

``` r
print(citation("spatialLIBD"), bibtex = TRUE)
#> 
#> To cite package 'spatialLIBD' in publications use:
#> 
#>   Pardo B, Spangler A, Weber LM, Hicks SC, Jaffe AE, Martinowich K,
#>   Maynard KR, Collado-Torres L (2022). "spatialLIBD: an R/Bioconductor
#>   package to visualize spatially-resolved transcriptomics data." _BMC
#>   Genomics_. doi:10.1186/s12864-022-08601-w
#>   <https://doi.org/10.1186/s12864-022-08601-w>,
#>   <https://doi.org/10.1186/s12864-022-08601-w>.
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
#> 
#>   Maynard KR, Collado-Torres L, Weber LM, Uytingco C, Barry BK,
#>   Williams SR, II JLC, Tran MN, Besich Z, Tippani M, Chew J, Yin Y,
#>   Kleinman JE, Hyde TM, Rao N, Hicks SC, Martinowich K, Jaffe AE
#>   (2021). "Transcriptome-scale spatial gene expression in the human
#>   dorsolateral prefrontal cortex." _Nature Neuroscience_.
#>   doi:10.1038/s41593-020-00787-0
#>   <https://doi.org/10.1038/s41593-020-00787-0>,
#>   <https://www.nature.com/articles/s41593-020-00787-0>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex},
#>     author = {Kristen R. Maynard and Leonardo Collado-Torres and Lukas M. Weber and Cedric Uytingco and Brianna K. Barry and Stephen R. Williams and Joseph L. Catallini II and Matthew N. Tran and Zachary Besich and Madhavi Tippani and Jennifer Chew and Yifeng Yin and Joel E. Kleinman and Thomas M. Hyde and Nikhil Rao and Stephanie C. Hicks and Keri Martinowich and Andrew E. Jaffe},
#>     year = {2021},
#>     journal = {Nature Neuroscience},
#>     doi = {10.1038/s41593-020-00787-0},
#>     url = {https://www.nature.com/articles/s41593-020-00787-0},
#>   }
```

Please note that the `spatialLIBD` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing the package.

## Interactive Websites

We provide the following interactive websites:

- Visium
  - [spatialDLPFC_Visium_Sp09](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp09):
    [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) website
    showing the spatially-resolved Visium data (n = 30) with statistical
    results comparing the Sp09 domains.
  - [spatialDLPFC_Visium_Sp16](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp16):
    similar but with the Sp16 domains.
  - [spatialDLPFC_Visium_Sp09_position](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp09_position):
    similar to *spatialDLPFC_Visium_Sp09* but with statistical results
    across the \_position_s (anterior, middle, posterior) adjusting for
    the Sp09 domains.
  - [spatialDLPFC_Visium_Sp09_pseudobulk](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp09_pseudobulk):
    [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1) website
    showing the pseudo-bulked Sp09 domains spatial data.
  - [spatialDLPFC_Visium_Sp16_pseudobulk](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp16_pseudobulk):
    similar to *spatialDLPFC_Visium_Sp09_pseudobulk* but with the Sp16
    domains data.
  - [spatialDLPFC_Visium_Sp28_pseudobulk](https://libd.shinyapps.io/spatialDLPFC_Visium_Sp28_pseudobulk):
    similar to *spatialDLPFC_Visium_Sp09_pseudobulk* but with the Sp16
    domains data.
- snRNA-seq
  - [spatialDLPFC_snRNA-seq](https://libd.shinyapps.io/spatialDLPFC_snRNA-seq):
    [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1) website
    showing the n = 19 snRNA-seq samples at single nucleus resolution.
- Visium SPG
  - [spatialDLPFC_Visium_SPG](https://libd.shinyapps.io/spatialDLPFC_Visium_SPG):
    [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) website
    showing the spatially-resolved data Visium SPG (n = 4) with
    statistical results comparing the Sp09 domains.
  - [spatialDLPFC Visium SPG on Loopy](https://loopybrowser.com/):
    [`loopy`](https://github.com/chaichontat/loopy-browser) website that
    allows to zoom in at the spot or cell level.

If you are interested in running the
[`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) applications
locally, you can do so thanks to the
[`spatialLIBD::run_app()`](http://research.libd.org/spatialLIBD/reference/run_app.html),
which you can also use with your own data as shown in our [vignette for
publicly available datasets provided by 10x
Genomics](http://bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/TenX_data_download.html).

``` r
## Run this web application locally
spatialLIBD::run_app()
## You will have more control about the length of the
## session and memory usage.
## You could also use this function to visualize your
## own data given some requirements described
## in detail in the package vignette documentation
## at http://research.libd.org/spatialLIBD/.
```

All of these websites are powered by open source software, namely:

- [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w)
- [`loopy`](https://github.com/chaichontat/loopy-browser)
- [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)

## Data Access

We highly value open data sharing and believe that doing so accelerates
science, as was the case between our
[`HumanPilot`](https://doi.org/10.1038/s41593-020-00787-0) and
[`BayesSpace`](https://doi.org/10.1038/s41587-021-00935-2) projects,
documented [on this
slide](https://speakerdeck.com/lcolladotor/hca-la-2022?slide=18). We
also value public questions, as they allow other users to learn from the
answers. If you have any questions, please ask them on a public forum
such as
[LieberInstitute/spatialDLPFC/issues](https://github.com/LieberInstitute/spatialDLPFC/issues).

### Processed Data

[`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w) also allows
you to access the data from this project. That is, a:

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
[LIBD rstats
club](https://docs.google.com/spreadsheets/d/1is8dZSd0FZ9Qi1Zvq1uRhm-P1McnJRd_zxdAfCRoMfA/edit?usp=sharing)
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

##### Access the data

Through the `spatialLIBD` package you can access the processed data in
it’s final R format.

### Processed data

Using `spatialLIBD` you can access the Human DLPFC spatial
transcriptomics data from the 10x Genomics Visium platform. For example,
this is the code you can use to access the spatially-resolved data. For
more details, check the help file for `fetch_data()`.

``` r
## Check that you have a recent version of spatialLIBD installed
stopifnot(packageVersion("spatialLIBD") >= "1.11.2")

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
#> colnames(113927): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ... TTGTTTGTATTACACG-1
#>   TTGTTTGTGTAAATTC-1
#> colData names(93): age array_col ... VistoSeg_count VistoSeg_percent
#> reducedDimNames(8): 10x_pca 10x_tsne ... HARMONY UMAP.HARMONY
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor

## Note the memory size
lobstr::object_size(spe)
#> 6.96 GB

## Set the cluster colors
colors_BayesSpace <- Polychrome::palette36.colors(28)
names(colors_BayesSpace) <- seq_len(28)

## Remake the logo image with histology information
p09 <- spatialLIBD::vis_clus(
    spe = spe,
    clustervar = "BayesSpace_harmony_09",
    sampleid = "Br6522_ant",
    colors = colors_BayesSpace,
    ... = " spatialDLPFC Human Brain - Sp09 domains\nMade with github.com/LieberInstitute/spatialDLPFC + spatiaLIBD"
)
p09
```

<img src="http://research.libd.org/spatialDLPFC/img/Br6255_ant_Sp09.png" width="800px" align="center" />

``` r
## Repeat but for Sp16
p16 <- spatialLIBD::vis_clus(
    spe = spe,
    clustervar = "BayesSpace_harmony_16",
    sampleid = "Br6522_ant",
    colors = colors_BayesSpace,
    ... = " spatialDLPFC Human Brain - Sp16 domains\nMade with github.com/LieberInstitute/spatialDLPFC + spatiaLIBD"
)
p16
```

<img src="http://research.libd.org/spatialDLPFC/img/Br6255_ant_Sp16.png" width="800px" align="center" />

### Raw data

You can access all the raw data through
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
- Slack channel: `libd_dlpfc_spatial`.

### Files: [`spatialDLPFC`](https://github.com/LieberInstitute/spatialDLPFC)

- `code`: R, python, and shell scripts for running various analyses.
  - `spot_deconvo`: cell-type deconvolution within Visium spots, enabled
    by tools like `tangram`, `cell2location`, and `cellpose`
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
  - `psychENCODE`: external data from PsychENCODE, originally retrieved
    from
    [here](https://www.synapse.org/#!Synapse:syn30108587.1/datasets/)
  - `sample_info`: spreadsheet with information about samples (sample
    ID, sample name, slide serial number, capture area ID)

This project is organized along the guidelines at
<https://lcolladotor.github.io/bioc_team_ds/organizing-your-work.html#.Yaf9fPHMIdk>.

### Files: [`DLPFC_snRNAseq`](https://github.com/LieberInstitute/DLPFC_snRNAseq)

TODO

### Other related files

- reference transcriptome from 10x Genomics:
  `/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/`

<!-- Google Tag Manager (noscript) -->
<noscript>
<iframe src="https://www.googletagmanager.com/ns.html?id=GTM-T57MPGR" height="0" width="0" style="display:none;visibility:hidden">
</iframe>
</noscript>
<!-- End Google Tag Manager (noscript) -->
<center>
<script type="text/javascript" id="clustrmaps" src="//clustrmaps.com/map_v2.js?d=tE-y1eDw-d-ZFWbu6cdHXSRgAwMRA-y9OobY9j8Krqk&cl=ffffff&w=a"></script>
</center>
