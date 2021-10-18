#$ -cwd
#$ -o ../../processed-data/02-our_data/02-run_tutorial.log
#$ -e ../../processed-data/02-our_data/02-run_tutorial.log
#$ -l mf=20G,h_vmem=20G
#$ -N run_tutorial_02

module load spagcn/1.2.0
python 02-run_tutorial.py
