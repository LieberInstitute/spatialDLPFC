#$ -cwd
#$ -o example_orig_nick_out/example_orig_nick_out.log
#$ -e example_orig_nick_out/example_orig_nick_out.log
#$ -l gpu,mf=32G,h_vmem=32G

#  This is also intended to be the directory from which this script is
#  submitted
main_dir=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spython/tangram_testing

cd $main_dir/Tangram/example
conda activate tangram

python $main_dir/example_orig_nick.py
