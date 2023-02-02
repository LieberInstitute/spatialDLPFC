#!/bin/bash

for PAIR in "EFNA5 FYN" "EFNA5 EPHA5" #"MBP PLP1" "EFNA5 FYN" "EFNA5 EPHA5"  "AQP4 GFAP" "SEMA6D PLXN4" "FSHB GPR20" "BMP3 BMPR1B" "LRFN5 PTPRD" "MDK LRP1" "DKK1 LRP1" "FYN PRKD1" "FYN ITGB7" "FYN NOX4" "FYN DCC" "FYN MST1R" "FYN LPP" "LRFN5 PTPRF" "FYN PTPRF" "FYN GRIN2A" "FYN GRIA1"

do
    cat << EOF > LR_sub_tmp.sh 
#!/bin/bash -l

# Simple SLURM sbatch example
#SBATCH --job-name=coloc-pairs
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=128G
#SBATCH --partition=cpu

set -e -o pipefail
 
ml purge 
ml Anaconda3
ml GEOS/3.9.1-GCC-11.2.0

source activate panpipes
echo $CONDA_DEFAULT_ENV
cd /camp/lab/gandhis/home/users/grantpm/LIBD_LR/grantpm/DLPFC_Visium_LIBD/processed-data/MGP_analysis/code/gene-risk-LR-analysis/03-colocalisation_analysis
echo "Start running python"

python3 03_c2l_colocalisation_analysis_top6.py $PAIR

EOF

    sbatch LR_sub_tmp.sh
    echo "Submitted job for $PAIR"
    sleep 0.5

done