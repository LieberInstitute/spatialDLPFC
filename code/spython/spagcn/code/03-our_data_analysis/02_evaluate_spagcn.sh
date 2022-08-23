#!/bin/bash
#$ -cwd
#$ -N "evaluate_spagcn_02"
#$ -o ../../processed-data/03-our_data_analysis/02_evaluate_spagcn.log
#$ -e ../../processed-data/03-our_data_analysis/02_evaluate_spagcn.log
#$ -l bluejay,mf=10G,h_vmem=10G

module load conda_R/4.1.x
Rscript 02_evaluate_spagcn.R
