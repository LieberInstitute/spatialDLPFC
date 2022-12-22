#!/bin/bash
#$ -cwd
#$ -N "nonIF"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=60G,h_vmem=60G,h_fsize=50G

cell_group="broad"
n_cells=0

#   Determine the log name for this run
log_dir="../../../processed-data/spot_deconvo/04-spotlight"
if [[ $n_cells -eq 0 ]]; then
    log_file="${log_dir}/02-nonIF_${cell_group}_full.log"
else
    log_file="${log_dir}/02-nonIF_${cell_group}_N${n_cells}.log"
fi

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

#   Use "0" to actually indicate the full dataset
if [[ $n_cells -eq 0 ]]; then
    echo "Running at $cell_group resolution without subsetting."
else
    echo "Running at $cell_group resolution, subsetting to $n_cells cells."
fi

module load conda_R/4.2.x
Rscript 02-nonIF.R -c $cell_group -n $n_cells

echo "**** Job ends ****"
date
} > $log_file 2>&1

