---
mode: 'agent'
description: 'Generate a corresponding shell script for a given R script'
---

In this project, R scripts are executed via a corresponding shell script of the
same name. The shell script contains SLURM-recognized syntax and submits the
R script (via the `Rscript` command) as a job to a SLURM-based computing cluster
to execute.

All corresponding shell scripts should have the following structure. Note the
brackets "<base name of R script here>", denoting the name of the corresponding
R script to be filled in. Note also the parent directory name of the R script
that must be filled in:

## Content template for shell scripts

```
#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=<base name of R script here>
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/<parent directory to R script>/logs/<base name of R script here>.txt
#SBATCH -e ../../processed-data/<parent directory to R script>/logs/<base name of R script here>.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.5

## List current modules for reproducibility
module list

Rscript <base name of R script here>.R

echo "**** Job ends ****"
date
```

## Concrete example

Suppose an R script `01_first.R` exists under `code/02_QC`, and you are asked to
generate the corresponding shell script. You will create
`code/02_QC/01_first.sh` whose log path will be
`../../processed-data/02_QC/logs/01_first.txt` (and of course fill in the other
placeholders from the template).

## Interpreting the prompt

If you are given just the name of an R script, generate the corresponding shell
script for that R script. If you are given an empty prompt but the R script is
attached as context, generate a shell script for that attached R script.
