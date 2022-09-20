#   Performed from the transfer node, interactively, for maximal speed:
#       jhpce-transfer01.jhsph.edu

cd /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
log_path=processed-data/synapse_upload/03-download/03-update_ucla_v2.log

module load synapse/2.6.0
python code/synapse_upload/03-download/03-update_ucla_v2.py > $log_path
