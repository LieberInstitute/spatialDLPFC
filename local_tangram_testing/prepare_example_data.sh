#  Download and unzip data required for the newer tutorial:
#  https://github.com/broadinstitute/Tangram/blob/master/example/1_tutorial_tangram.ipynb
mkdir data/
wget https://storage.googleapis.com/tommaso-brain-data/tangram_demo/mop_sn_tutorial.h5ad.gz -O data/mop_sn_tutorial.h5ad.gz
wget https://storage.googleapis.com/tommaso-brain-data/tangram_demo/slideseq_MOp_1217.h5ad.gz -O data/slideseq_MOp_1217.h5ad.gz
wget https://storage.googleapis.com/tommaso-brain-data/tangram_demo/MOp_markers.csv -O data/MOp_markers.csv
gunzip -f data/mop_sn_tutorial.h5ad.gz
gunzip -f data/slideseq_MOp_1217.h5ad.gz
