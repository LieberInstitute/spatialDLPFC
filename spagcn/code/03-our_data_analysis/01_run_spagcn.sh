#!/bin/bash
#$ -cwd
#$ -N "run_spagcn_01"
#$ -o ../../processed-data/03-our_data_analysis/01_run_spagcn_$TASK_ID.log
#$ -e ../../processed-data/03-our_data_analysis/01_run_spagcn_$TASK_ID.log
#$ -l mf=60G,h_vmem=60G
#$ -t 3
#$ -tc 1

module load spagcn/1.2.0
python 01_run_spagcn.py -i $SGE_TASK_ID
