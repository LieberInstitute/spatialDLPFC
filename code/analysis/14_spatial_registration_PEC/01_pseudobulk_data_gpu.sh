#!/bin/bash

## Usage:
# sh 01_pseudobulk_data_gpu.sh

## Create the logs directory
mkdir -p logs

for PE_data in CMC DevBrain-snRNAseq IsoHuB MultiomeBrain-DLPFC SZBDMulti-Seq UCLA-ASD LIBD PTSDBrainomics; do

    ## Internal script name
    SHORT="01_pseudobulk_data_${PE_data}"
	NAME="pseudobulk_data_${PE_data}"

    # Construct shell file
    echo "Creating script 01_pseudobulk_data_${PE_data}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l caracol,mem_free=50G,h_vmem=50G,h_fsize=100G
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
module load conda_R/4.2.x

## List current modules for reproducibility
module list

## GPU check
USAGE_CUTOFF=10
NUM_GPUS=1

avail_gpus=\$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader | cut -d " " -f 1 | awk -v usage="\$USAGE_CUTOFF" '\$1 < usage {print NR - 1}')

#  Simply exit with an error if there are no GPUs left
if [[ -z \$avail_gpus ]]; then
    echo "No GPUs are available."
    exit 1
fi

export CUDA_VISIBLE_DEVICES=\$(echo "\$avail_gpus" | head -n \$NUM_GPUS | paste -sd ",")

## Edit with your job command
Rscript 01_pseudobulk_data.R version5 ${PE_data}_annotated.h5ad

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

##    call="qsub .${SHORT}.sh"
##    echo $call
##    $call
done
