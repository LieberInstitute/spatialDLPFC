#!/bin/bash

## Usage:
# sh 05_compute_registration_Dx.sh

## Create the logs directory
mkdir -p logs

for dataset in CMC DevBrain-snRNAseq SZBDMulti-Seq UCLA-ASD MultiomeBrain-DLPFC LIBD PTSDBrainomics; do

    ## Internal script name
    SHORT="05_compute_registration_Dx_${dataset}"
	NAME="compute_registration_Dx_${dataset}"

    # Construct shell file
    echo "Creating script 05_compute_registration_Dx_${dataset}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N ${NAME}
#$ -o logs/${SHORT}.txt
#$ -e logs/${SHORT}.txt
#$ -hold_jid pseudobulk_data_${dataset}
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
module load conda_R/4.2.x

## List current modules for reproducibility
module list

## Sumbit Rscript
Rscript 05_compute_registration_Dx.R ${dataset}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
