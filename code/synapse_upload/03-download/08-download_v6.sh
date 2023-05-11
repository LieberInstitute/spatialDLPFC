#   Performed from the transfer node, interactively, for maximal speed:
#       jhpce-transfer01.jhsph.edu

cd /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
log_path=processed-data/synapse_upload/03-download/08-download_v6.log
out_dir=raw-data/psychENCODE/version6
credential_path=/users/neagles/synapse_credentials.yaml

mkdir -p $out_dir

{
module load synapse/2.6.0

synapse \
    -p $(grep "^token:" $credential_path | cut -d " " -f 2 | tr -d '"') \
    get \
        -r syn51123274 \
        --downloadLocation $out_dir

echo "Done at $(date)."
} > $log_path 2>&1
