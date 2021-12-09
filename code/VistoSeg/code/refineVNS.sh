#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -o /users/hdivecha/VistoSeg/code/Logs/$TASK_ID.txt
#$ -e /users/hdivecha/VistoSeg/code/Logs/$TASK_ID.txt
#$ -m e
#$ -M heenadivecha@gmail.com
#$ -t 2
#$ -tc 4


 
echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(awk 'BEGIN {FS="\t"} {print $1}' /users/hdivecha/VistoSeg/code/listOffiles.txt | awk "NR==${SGE_TASK_ID}")"
echo "****"


##module load matlab/R2019a

toolbox='/users/hdivecha/VistoSeg/code'


fname=$(awk 'BEGIN {FS="\t"} {print $1}' /users/hdivecha/VistoSeg/code/listOffiles.txt | awk "NR==${SGE_TASK_ID}")
echo $fname
M=$(awk 'BEGIN {FS="\t"} {print $3}' /users/hdivecha/VistoSeg/code/listOffiles.txt | awk "NR==${SGE_TASK_ID}")
echo $M

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), refineVNS('$fname','$M')"


echo "**** Job ends ****"
date

