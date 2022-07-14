#   Symbolically link high-resolution histology images for the 12 NN DLPFC
#   Visium samples into raw-data, for use in cell-type spot-level deconvolution
#   with Tangram

#  Get absolute path to 'spython' repo
base_dir=$(git rev-parse --show-toplevel)

sample_path=$base_dir/tangram_libd/processed-data/03_nn_run/brain_samples.txt

src_dir=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X
dest_dir=$base_dir/tangram_libd/raw-data/03_nn_run/hires_histology
mkdir -p $dest_dir

for id in $(cat $sample_path); do
    ln -s $src_dir/$id/tissue_hires_image.png $dest_dir/${id}.png
done
