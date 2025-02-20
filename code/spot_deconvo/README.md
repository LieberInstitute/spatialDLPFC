Code for deconvolving cell types within Visium spots (for both IF and non-IF samples). `tangram` and `cell2location` are essentially alternatives used here, where snRNA-seq data is spatially aligned onto Visium spots, and prior information about cell counts for each spot allow the classification of discrete cells in each spot. `cellpose` is currently used to segment the DAPI channel of IF images, ultimately generating cell counts (which can be fed to `tangram`/`cell2location` as input). Fluorescence intensities in other channels of the IF images can be used to deduce cell types for each nucleus (and therefore cell).

# Software Management

At JHPCE, each python tool is managed in its own python virtual environment, which can be accessed by loading a particular [module](https://github.com/LieberInstitute/jhpce_mod_source) prior to opening a python interactive session or running a python script. The modules used are:

- `cellpose/2.0`
- `tangram/1.0.2`
- `cell2location/0.8a0`

# Computational Workflow

Spot deconvolution is performed following the below workflow. The code directories used for each portion are listed in black font on the left. See the README in `02-cellpose` to describe steps used to segment IF images and train a CART to classify cell types as an alternative to spot deconvolution via `tangram`, `cell2location`, and `SPOTlight`. These steps are performed prior to the workflow shown below.

![Computational Workflow](computational_workflow.png)
