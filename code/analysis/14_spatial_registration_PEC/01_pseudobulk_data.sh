#!/bin/bash

## Usage:
# sh 01_pseudobulk_data.sh

## Create the logs directory
mkdir -p logs

for PE_data in CMC/CMC-CellHashing_annotated.h5ad DevBrain-snRNAseq/DevBrain-snRNAseq_annotated.h5ad IsoHuB/IsoHuB-snRNAseq_annotated.h5ad SZBDMulti-Seq/SZBDMulti-Seq_annotated.h5ad UCLA-ASD/UCLA-ASD-snRNAseq_annotated_mismatches_removed.h5ad Urban-DLPFC/Urban-DLPFC-snRNAseq_annotated.h5ad; do
	
	## Use dirname from file
	PE_data2=(${PE_data//// })
	PE_data2=${PE_data2[0]}
	## Internal script name
    SHORT="01_pseudobulk_data_${PE_data2}"
	NAME="pseudobulk_data_${PE_data2}" 

    # Construct shell file
    echo "Creating script 01_pseudobulk_data_${PE_data2}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N ${NAME}
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
module load conda_R/4.2

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_pseudobulk_data.R ${PE_data}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
