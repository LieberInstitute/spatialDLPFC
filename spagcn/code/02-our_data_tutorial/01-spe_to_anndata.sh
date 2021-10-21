#$ -cwd
#$ -o ../../processed-data/02-our_data_tutorial/01-spe_to_anndata.log
#$ -e ../../processed-data/02-our_data_tutorial/01-spe_to_anndata.log
#$ -l mf=20G,h_vmem=20G
#$ -N spe_to_anndata_01

module load conda_R/4.1
Rscript 01-spe_to_anndata.R
