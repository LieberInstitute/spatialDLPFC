library(edgeR)
library(scater)
library(scran)

#load pseudobulked object 

#http://bioconductor.org/books/3.14/OSCA.multisample/multi-sample-comparisons.html#creating-pseudo-bulk-samples
# label <- 1
# current <- spe_pseudo[,label == spe_pseudo$bayesSpace_harmony_3]
# y <- DGEList(counts(current), samples=colData(current))


#loop through with pseudoBulkDGE
dim(spe_pseudo)
spe_pseudo <- spe_pseudo[,spe_pseudo$ncells >= 10]
dim(spe_pseudo)

head(model.matrix(~factor(region) + factor(subject),as.data.frame(colData(spe_pseudo)))) #;ppl at rclub session on model matrix
#possibly have to loop through other coefficients other than region such as subject, age, 
de.results <- pseudoBulkDGE(spe_pseudo, 
                            label=spe_pseudo$BayesSpace,
                            design=~factor(region) + factor(subject),
                            coef = "factor(region)middle" #comes from topTable from limma, specifies the coefficient you want to do the t-test on
                            
                            
)

#this is pairwise comparison. can see how they did it before with limma for the pilot study. 


head(de.results)
head(de.results[[1]])
# DataFrame with 6 rows and 5 columns
# logFC    logCPM         F    PValue       FDR
# <numeric> <numeric> <numeric> <numeric> <numeric>
#   ENSG00000243485        NA        NA        NA        NA        NA
# ENSG00000238009        NA        NA        NA        NA        NA
# ENSG00000239945        NA        NA        NA        NA        NA
# ENSG00000241860        NA        NA        NA        NA        NA
# ENSG00000229905        NA        NA        NA        NA        NA
# ENSG00000237491  0.267683   7.50472  0.566731  0.452311  0.899136

dim(de.results[[1]])
table(de.results[[1]]$FDR < 0.05) #shows we 4 differentially expressed genes between middle and anterior in cluster one 
