module load loopy/1.0.0-next.8

repo_dir=$(git rev-parse --show-toplevel)
spot_diameter_micro=55 # 55-micrometer diameter for Visium spot

img_path=$repo_dir/raw-data/Images/VisiumIF/VistoSeg/V10B01-087_A1.tif
json_path=$repo_dir/processed-data/01_spaceranger_IF/Br2720_Ant_IF/outs/spatial/scalefactors_json.json

#   Grab the spot diameter in pixels from the spaceranger JSON
spot_diameter_px=$(cat $json_path | tr , '\n' | tr -d "{  }" | grep "spot_diameter_fullres" | cut -d ":" -f 2)

#   Calculate meters per pixel for a spot, assuming the result is smaller than
#   e-6 (a "0" is prepended)
m_per_px="0$(echo "$spot_diameter_micro / $spot_diameter_px" | bc -l)e-6"

#   Produce the directory that can be dragged onto loopybrowser.com to display
#   image
mkdir -p $repo_dir/processed-data/spot_deconvo/07-loopy

loopy image $img_path \
    --scale $m_per_px \
    --channels Lipofuscin,DAPI,GFAP,NeuN,OLIG2,TMEM119 \
    -o $repo_dir/processed-data/spot_deconvo/07-loopy

#   Download directory 'V10B01-087_A1' and drag into local browser
