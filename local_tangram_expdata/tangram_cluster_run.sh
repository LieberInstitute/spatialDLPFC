#!/bin/bash
#$ -cwd
#$ -N "tangram-cluster-gpu-expdata"
#$ -o logs/tangram_expdata.log
#$ -e logs/tangram_expdata.err
#$ -j y
#$ -l gpu,mf=32G,h_vmem=32G

#  This is also intended to be the directory from which this script is
#  submitted
main_dir=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spython/local_tangram_expdata

conda activate tangram

python $main_dir/tangram_exp_arta.py
