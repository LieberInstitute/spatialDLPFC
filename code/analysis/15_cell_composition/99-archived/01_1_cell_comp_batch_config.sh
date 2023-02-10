#!/bin/bash
## TODO: edit this section
#$ -l mem_free=5G,h_vmem=5G,h_fsize=1G
##$ -cwd
#$ -m e
#$ -M bguo6@jhu.edu

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"


## Load the R module
module load conda_R/4.2.x

## List current modules for reproducibility
module list


# Create a directory to store your log files
# Import to debuging
logPath=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/analysis/15_cell_composition/logs/
mkdir -p $logPath

# Run simulation with simulation parameters using R batch --args flag
###      The $ sign helps to fetch the simulation parameters passed from start_sim.R\
###      If the parameter is a string (a vector of string), you need to quote the string with single quote `
R CMD BATCH --no-save "--args deconv_method='$deconv_method'" /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/deconv_method/01_2_cell_comp_EDA.R  $logPath${deconv_method}.Rout


echo "**** Job ends ****"
date
