#!/bin/bash

## Usage:
# sh 01_pseudobulk_data.sh

## Create the logs directory
mkdir -p logs

for PE_data in CMC-CellHashing_annotated.h5ad DevBrain-snRNAseq_annotated.h5ad IsoHuB-snRNAseq_annotated.h5ad SZBDMulti-Seq_annotated.h5ad UCLA-ASD-snRNAseq_annotated_mismatches_removed.h5ad Urban-DLPFC-snRNAseq_annotated.h5ad; do

    ## Internal script name
    SHORT="01_pseudobulk_data_${PE_data}"

    # Construct shell file
    echo "Creating script 01_pseudobulk_data_${PE_data}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N ${SHORT}
#$ -o logs/${SHORT}.txt
#$ -e logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript -e "options(width = 120); print('${PE_data}'); sessioninfo::session_info()"

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
