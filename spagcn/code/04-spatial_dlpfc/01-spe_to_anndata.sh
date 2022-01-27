#$ -cwd
#$ -o ../../processed-data/04-spatial_dlpfc/01-spe_to_anndata.log
#$ -e ../../processed-data/04-spatial_dlpfc/01-spe_to_anndata.log
#$ -l mf=20G,h_vmem=20G
#$ -N spe_to_anndata_01

module load conda_R/4.1.x
Rscript 01-spe_to_anndata.R
