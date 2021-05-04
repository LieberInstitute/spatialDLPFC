#!/bin/bash
#$ -cwd
#$ -N "tangram-cluster-expdata"
#$ -o logs/tangram_expdata.log
#$ -e logs/tangram_expdata.err
#$ -j y
#$ -l gpu,mf=45G,h_vmem=45G

# -l bluejay,mem_free=64G,h_vmem=64G,h_fsize=100G



#  This is also intended to be the directory from which this script is
#  submitted
main_dir=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spython/local_tangram_expdata

ml conda
conda activate tangram

python $main_dir/tangram_exp_arta.py
