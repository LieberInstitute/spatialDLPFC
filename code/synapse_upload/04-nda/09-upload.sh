#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=3G
#SBATCH --job-name=09-upload
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../../../processed-data/synapse_upload/04-nda/logs/09-upload.txt
#SBATCH -e ../../../processed-data/synapse_upload/04-nda/logs/09-upload.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load nda-tools/0.2.27
module list

repo_dir=$(git rev-parse --show-toplevel)
meta_dir=$repo_dir/processed-data/synapse_upload/04-nda
data_dir=$repo_dir/processed-data/synapse_upload/04-nda/to_upload

vtcmd \
    -l $data_dir \
    -b \
    -c C5229 \
    -t "LIBD spatial DLPFC data submission" \
    -d "All imaging, genotype, and FASTQ data for the spatialDLPFC project" \
    -u nickeagles77 \
    $meta_dir/rna_seq.csv $meta_dir/snp_array.csv $meta_dir/visium_image.csv $meta_dir/genomics_subject.csv

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
