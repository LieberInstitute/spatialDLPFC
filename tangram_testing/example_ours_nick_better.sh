#$ -cwd
#$ -o example_ours_nick_better_out/example_ours_nick_better_out.log
#$ -e example_ours_nick_better_out/example_ours_nick_better_out.log
#$ -l gpu,mf=64G,h_vmem=64G

#  This is also intended to be the directory from which this script is
#  submitted
main_dir=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spython/tangram_testing

conda activate tangram

python $main_dir/example_ours_nick_better.py
