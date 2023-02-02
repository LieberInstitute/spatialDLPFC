###-------------------------------------------- OVERVIEW
This folder contains all the the analyses relating to ligand-receptor role of genes with positive genetic association in SCZ for LIBD DLPFC project. 
Analysis approach conceived by Melissa Grant-Peters under the supervision of Mina Ryten
Code written by Melissa Grant-Peters

###-------------------------------------------- DATA
Publicly available data: 
- OpenTargets platform 
- GTEx
- Omnipath DB
Newly generated data:
- snRNAseq and 10x Visium from DLPFC

###-------------------------------------------- ANALYSES
00-OpenTargets_SCZ_risk_genes
    Formats gene lists exported from OpenTargets platform for this analysis
    
01-LR_occurence_bootstrapped
    Compares the incidence of ligands and receptors in a disease risk gene list to occurance in random CNS-expressed gene list (obtained from GTEx). Random CNS-expressed gene list is bootstrapped. 
    
02-SCZ_LRs
    Characterises the incidence of ligands and receptors in SCZ risk gene list (% of genes, molecule type). 
    Outputs a list of ligand-receptor interactions where both LRs are associated with SCZ risk. 

03-colocalisation_analysis
    Runs in parallellised fashion an analysis for:
        - Which 10x Visium spots have co-expression (counts>0) for a given pair of genes
        - What is the cell type neighbourhood for these spots (uses cell2location deconvoluted outputs) 
                co-localisation of two cell types is recorded in an adjacency matrix and presented as a heatmap
        - Takes the ratio of these results in relation to the spots with counts == 0 for both genes inquired. 
                Final result is a heatmap

04-snRNAseq_specificity
    Assesses cell-type specificity for the genes outputted in 02-SCZ_LRs using a tau statistic and SPM statistic applied to snRNAseq data

requirements
    Contains a pip freeze of the contents of the environment used for this analysis. This is not limited to the libraries actively used and may include other libraries as well. 