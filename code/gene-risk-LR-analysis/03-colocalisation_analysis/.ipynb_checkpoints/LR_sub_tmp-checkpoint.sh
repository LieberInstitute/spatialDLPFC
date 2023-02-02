#!/bin/bash -l

# Simple SLURM sbatch example
#SBATCH --job-name=coloc-pairs
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=128G
#SBATCH --partition=cpu

set -e -o pipefail
 
ml purge 
ml Anaconda3
ml GEOS/3.9.1-GCC-11.2.0

source activate panpipes
cd /camp/lab/gandhis/home/users/grantpm/LIBD_LR/grantpm/DLPFC_Visium_LIBD/processed-data/MGP_analysis/code/gene-risk-LR-analysis/03-colocalisation_analysis
echo "Start running python"

python3 03_c2l_colocalisation_analysis_top3.py EFNA5 EPHA5

