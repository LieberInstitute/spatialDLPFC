#   Performed from the transfer node, interactively, for maximal speed:
#       jhpce-transfer01.jhsph.edu

cd /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
log_path=processed-data/synapse_upload/03-download/05-download_v3.log

module load synapse/2.6.0
python code/synapse_upload/03-download/05-download_v3.py > $log_path
