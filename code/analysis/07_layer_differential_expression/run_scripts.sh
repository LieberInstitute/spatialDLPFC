#!/bin/bash

echo "**** Job starts ****"
date


cd /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/analysis/07_layer_differential_expression

rm logs/01*.txt
qsub 01_create_pseudobulk_data.sh

rm logs/02*.txt
qsub 02_explore_expr_variability.sh

rm logs/03*.txt
qsub 03_model_BayesSpace.sh

echo "**** Job ends ****"
date
