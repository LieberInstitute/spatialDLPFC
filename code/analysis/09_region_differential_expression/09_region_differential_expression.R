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

de.results <- pseudoBulkDGE(spe_pseudo, 
                            label=spe_pseudo$bayesSpace_harmony_3,
                            design=~factor(region) + factor(subject),
                            coef = 
                            
)