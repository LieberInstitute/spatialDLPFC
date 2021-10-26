#!/bin/bash

git submodule add https://github.com/broadinstitute/Tangram tangram_libd/Tangram

#  Since conda environments are user-specific, each user planning to run
#  tangram should run these commands. This produces a "conda environment", a
#  particular version of python and pip along with python packages required as
#  dependencies for tangram.

conda env create -n tangram -f Tangram/environment.yml

conda activate tangram
pip install tangram-sc # install the tangram python package
pip install pyhere     # python's equivalent to 'here' in R
pip install numpy==1.20
conda deactivate

#  To later load the environment:
#      conda activate tangram
#
#  Now when python or python scripts are run, the correct python version and
#  set of packages is automatically used/ available, isolated from any other
#  python installations which may be present.
#
#  To restore the original environment:
#      conda deactivate
