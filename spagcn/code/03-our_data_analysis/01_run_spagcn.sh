#!/bin/bash
#$ -cwd
#$ -N "run_spagcn_01"
#$ -o ../../processed-data/03-our_data_analysis/01_run_spagcn_$TASK_ID.log
#$ -e ../../processed-data/03-our_data_analysis/01_run_spagcn_$TASK_ID.log
#$ -l mf=20G,h_vmem=20G
#$ -t 1
#$ -tc 1

module load spagcn/1.2.0
python 01_run_spagcn.py -i $SGE_TASK_ID
