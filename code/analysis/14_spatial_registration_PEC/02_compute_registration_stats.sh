#!/bin/bash

## Usage:
# sh 02_compute_registration_stats.sh

## Create the logs directory
mkdir -p logs

for dataset in CMC DevBrain-snRNAseq IsoHuB SZBDMulti-Seq UCLA-ASD LIBD PTSDBrainomics; do

    ## Internal script name
    SHORT="02_compute_registration_stats_${dataset}"
	NAME="compute_registration_stats_${dataset}"

    # Construct shell file
    echo "Creating script 02_compute_registration_stats_${dataset}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
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
Rscript 02_compute_registration_stats.R ${dataset}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
