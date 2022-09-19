## Documenting the analysis_IF folder for spatialDLPFC study

internal filepath: /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/analysis_IF
github: https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/analysis_IF

Visium_IF was performed on 4 samples as a part of the spatialDLPFC study. This folder contains the analysis performed for those samples. There are two subfolders here:

## 1. 01_build_spe_IF: 
internal: /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/analysis_IF/01_build_spe_IF
Contains the script used to generate the basic spe object (01_build_basic_spe.R) and the bash script which was used to run the .R script(01_build_basic_spe.sh). output of the bash script is saved in the logs folder in a .txt file (01_build_basic_spe.txt)   

## 2. 02_shinyapp
internal: /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/analysis_IF/02_shinyapp
Contains the script used to generate the script for generating and deploying the spatialLIBD shiny app. app.R script was used to generate the app and deploy.R script was used to deploy the app. Some huge files were ignored using the .gitignore file.
