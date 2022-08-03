#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G
#$ -N "merge_meta"
#$ -o ../../../processed-data/synapse_upload/02-metadata_files/05-merge_meta.log
#$ -e ../../../processed-data/synapse_upload/02-metadata_files/05-merge_meta.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

#   Copy assay and individual metadata files to the final location
dest_dir=../../../processed-data/synapse_upload/02-metadata_files/combined
spatial_dir=../../../processed-data/synapse_upload/02-metadata_files
sc_dir=/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/04_synapse_upload

cp $sc_dir/assay.csv $dest_dir/assay_scrnaSeq.csv
cp $spatial_dir/assay.csv $dest_dir/assay_TODO.csv

cp $spatial_dir/individual.csv $dest_dir/individual.csv

module load conda_R/devel
Rscript 05-merge_meta.R

echo "**** Job ends ****"
date
