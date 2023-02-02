#!/bin/bash -l

#####################################################################
#####   Gene risk LR analysis written by Melissa Grant-Peters   #####
#####            Lieber Institute spatial DLPFC project         #####
#
#
# In order to run this full analysis
# 1) load a python environment with requirements.txt
# 2) make sure that scripts accessing data (analyses in 03, 04 and 05) are consistent with the user's hierarchy
# 3) run script with: bash run_all.sh
# 4) Generated files will be stored in processed-data or plots directories

######################################################################
echo "Running gene risk ligand-receptor analysis for schizophrenia targets"


echo "Formatting OpenTargets gene list..."
cd 00-OpenTargets_SCZ_risk_genes
python3 00-OpenTargets_format_gene_list.py

echo "Calculating occurence of ligands and receptors in risk genes..."
cd ../01-LR_occurrence_bootstrapped
python3 01-LR_occurence_bootstrapped.py

echo "Determining risk ligand-receptor interactions for Schizophrenia..."
cd ../02-SCZ_LRs
python3 02-SCZ_LRs.py

echo "Calculating cell type neighborhood of interactions..."
cd ../03-colocalisation_analysis
python3 03_prepare_data.py
python3 03_c2l_colocalisation_analysis_top3.py EFNA5 EPHA5
python3 03_c2l_colocalisation_analysis_top6.py EFNA5 EPHA5
python3 03_c2l_colocalisation_analysis_top3.py EFNA5 FYN
python3 03_c2l_colocalisation_analysis_top6.py EFNA5 FYN
python3 03_tangram_colocalisation_analysis.py

echo "Determining cell type specificity of targets..."
cd ../04-specificity_LR_targets
python3 04-specificity_LR_targets.py
python3 04-visualise_specificity_EFNA5-EPHA5-FYN.py

echo "Calculating rates of intracellular coexpression"
cd ../05-intracellular_coexpression
python3 05-intracellular_coexpression.py

echo "Done."
