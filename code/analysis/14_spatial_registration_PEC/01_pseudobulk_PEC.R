### read in more datasets, format, pseudobulk and save objects. 

library('zellkonverter')
library(SingleCellExperiment)
library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(rafalib)
library(scuttle)
library(limma)
library(RColorBrewer)
library(lattice)
library(edgeR)

############################
#load CMC annotated snRNAseq data
###########################
sce = readH5AD('/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version2/CMC/CMC-snRNAseq_annotated.h5ad')

# identify annotation/cluster labels
rowData(sce)$gene_name <- rownames(rowData(sce)) #save gene name as column of rowData
rownames(sce) <- rowData(sce)$featureid # have to make row names of object the ensembl id instead of gene names
colData(sce)$anno <- as.factor(colData(sce)$anno)
colData(sce)$demux_type <- as.factor(colData(sce)$demux_type)
colnames(colData(sce))[9] <- "sample_id"
names(assays(sce)) <- "counts"

#save sce here because it takes so long to convert the h5AD file
saveRDS(sce,file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/CMC_sce.rds")

message("Make psuedobulk object")
spe_pseudo <- aggregateAcrossCells(
  sce,
  DataFrame(
    BayesSpace = sce$anno,
    sample_id = sce$sample_id
  )
)
colnames(spe_pseudo) <- paste0(spe_pseudo$sample_id, "_", spe_pseudo$anno)

message("Filter lowly expressed genes")
rowData(spe_pseudo)$filter_expr <- filterByExpr(spe_pseudo)
summary(rowData(spe_pseudo)$filter_expr)

spe_pseudo <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr), ]

message("Normalize expression")
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
## We don't need this 'x' object anymore
rm(x)

#save pseudobulked object
save(spe_pseudo, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/CMC_spe_pseudo.rda")

#######################
###load DevBrain dataset 
########################
sce.dev.brain <- readH5AD(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/psychENCODE/version2/DevBrain/DevBrain-snRNAseq_annotated.h5ad")

# identify annotation/cluster labels
rowData(sce.dev.brain)$gene_name <- rownames(rowData(sce.dev.brain)) #save gene name as column of rowData
rownames(sce.dev.brain) <- rowData(sce.dev.brain)$featureid # have to make row names of object the ensembl id instead of gene names
colData(sce.dev.brain)$anno <- as.factor(colData(sce.dev.brain)$anno)
colnames(colData(sce.dev.brain))[7] <- "sample_id"
names(assays(sce.dev.brain)) <- "counts"

#save sce here because it takes so long to convert the h5AD file
saveRDS(sce.dev.brain,file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/DevBrain_sce.rds")

message("Make psuedobulk object")
spe_pseudo <- aggregateAcrossCells(
  sce.dev.brain,
  DataFrame(
    Anno = sce.dev.brain$anno,
    sample_id = sce.dev.brain$sample_id
  )
)
colnames(spe_pseudo) <- paste0(spe_pseudo$sample_id, "_", spe_pseudo$anno)

message("Filter lowly expressed genes")
rowData(spe_pseudo)$filter_expr <- filterByExpr(spe_pseudo)
summary(rowData(spe_pseudo)$filter_expr)

spe_pseudo <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr), ]

message("Normalize expression")
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))
## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
## We don't need this 'x' object anymore
rm(x)

#save pseudobulked object
saveRDS(spe_pseudo, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/DevBrain_spe_pseudo.rds")

