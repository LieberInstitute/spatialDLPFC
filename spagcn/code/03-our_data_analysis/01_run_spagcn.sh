#!/bin/bash
#$ -cwd
#$ -N "main_analysis_03"
#$ -o ../../processed-data/02-our_data/03_main_analysis_$TASK_ID.log
#$ -e ../../processed-data/02-our_data/03_main_analysis_$TASK_ID.log
#$ -l mf=20G,h_vmem=20G
#$ -t 1-12
#$ -tc 1

module load spagcn/1.2.0
python 03_main_analysis.py -i $SGE_TASK_ID
