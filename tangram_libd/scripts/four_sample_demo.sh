#$ -cwd
#$ -o four_sample_demo_out/$TASK_ID.log
#$ -e four_sample_demo_out/$TASK_ID.log
#$ -l gpu,mf=64G,h_vmem=64G
#$ -t 1-4
#$ -tc 1

conda activate tangram
python four_sample_demo.py -i $SGE_TASK_ID
