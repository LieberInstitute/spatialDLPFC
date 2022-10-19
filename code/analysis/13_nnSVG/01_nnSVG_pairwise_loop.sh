#!/bin/bash

## Usage:
# sh 01_nnSVG_pairwise_loop.sh

## Create the logs directory
mkdir -p logs

for domains in 5v9 4v16 7v13 12v13 7v12 12v16; do

    ## Internal script name
    SHORT="01_nnSVG_pairwise_loop_${domains}"
    NAME="nnSVG_pairwise_loop_${domains}"

    # Construct shell file
    echo "Creating script 01_nnSVG_pairwise_loop_${domains}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -pe local 5
#$ -N ${NAME}
#$ -o logs/${SHORT}.\$TASK_ID.txt
#$ -e logs/${SHORT}.\$TASK_ID.txt
#$ -m e
#$ -t 1-30
#$ -tc 10

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
Rscript 01_nnSVG_pairwise.R ${domains} $SGE_TASK_ID

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
