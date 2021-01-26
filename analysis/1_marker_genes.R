library(here)
library(tidyverse)
library(ggplot2)
library(SingleCellExperiment)
library(Matrix)
library(rjson)
library(Seurat)
library(rtracklayer)
library(scran)
library(scater)
library(readbitmap)
library(grid)


sample_names <- c("DLPFC_Br2743_ant_manual_alignment", "DLPFC_Br2743_mid_manual_alignment","DLPFC_Br2743_post_manual_alignment","DLPFC_Br3942_ant_manual_alignment","DLPFC_Br3942_mid_manual_alignment","DLPFC_Br3942_post_manual_alignment","DLPFC_Br6423_ant_manual_alignment","DLPFC_Br6423_mid_manual_alignment","DLPFC_Br6423_post_manual_alignment","DLPFC_Br8492_ant_manual_alignment","DLPFC_Br8492_mid_manual_alignment","DLPFC_Br8492_post_manual_alignment")
sce_list <- vector("list", length = length(sample_names))
names(sce_list) <- sample_names

# get GTF, this seems like what they used
gtf = rtracklayer::import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
gtf = gtf[gtf$type	== "gene"]
names(gtf) = gtf$gene_id


for (i in seq_along(sample_names)) {
  
  # ---------
  # load data
  # ---------
  
  # select sample
  sample_name <- sample_names[i]
  #sample_name <- sample_names[3]
  
  # path to Space Ranger output files
  if (Sys.info()["sysname"] == "Linux") {
    # files on JHPCE cluster
    dir_outputs <- "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/outputs/NextSeq"
  } else if (Sys.info()["sysname"] == "Darwin") {
    # copy of files on Mac laptop
    dir_outputs <- "~/DLPFC"
  }
  
  # note: using "filtered" barcodes list containing only spots over tissue
  dir_matrix <- file.path(dir_outputs, sample_name, "outs")
  
  # barcodes
  file_barcodes <- file.path(dir_matrix, "filtered_feature_bc_matrix/barcodes.tsv.gz")
  df_barcodes <- read.csv(file_barcodes, sep = "\t", header = FALSE, 
                          col.names = c("barcode_id"))
  # features
  file_features <- file.path(dir_matrix, "filtered_feature_bc_matrix/features.tsv.gz")
  df_features <- read.csv(file_features, sep = "\t", header = FALSE, 
                          col.names = c("gene_id", "gene_name", "feature_type"))
  
  # > dim(df_features)
  # [1] 36601     3
  # # add more gene information from gtf file 
  
  x <- intersect(gtf$gene_id,df_features$gene_id)
  gtf = gtf[x]
  seqlevels(gtf)[seq_len(25)] = paste0("chr", seqlevels(gtf)[seq_len(25)])

  
  # counts
  file_counts <- file.path(dir_matrix, "filtered_feature_bc_matrix.h5")
  #counts <- readMM(file = file_counts)
  counts <- Read10X_h5(file_counts)
  
  
  # spatial scale factors
  dir_spatial <- file.path(dir_outputs, sample_name, "outs", "spatial")
  file_scale <- file.path(dir_spatial, "scalefactors_json.json")
  scalefactors <- fromJSON(file = file_scale)
  
  # spatial coordinates
  file_tisspos <- file.path(dir_spatial, "tissue_positions_list.csv")
  df_tisspos <- read.csv(file_tisspos, header = FALSE, 
                         col.names=c("barcode","tissue","row","col","imagerow","imagecol"))
  df_tisspos$imagerow <-df_tisspos$imagerow * scalefactors$tissue_lowres_scalef    # scale tissue coordinates for lowres image
  df_tisspos$imagecol <- df_tisspos$imagecol * scalefactors$tissue_lowres_scalef
  

  
   # check dimensions
  dim(df_barcodes)
  dim(df_features)
  dim(counts)
  # note df_tisspos contains all spots (not filtered) - need to match later
  dim(df_tisspos)
  
  # image paths
  # library(readbitmap)
  # image_path <- file.path(dir_spatial,"tissue_lowres_image.png")
  # image_cl<- read.bitmap(image_path)
  # dims = dim(image_cl)
  # dims = as.data.frame(t(dims))
  # colnames(dims) = c("height", "width", "channel")
  # 
  # library(grid)
  # grobs<-rasterGrob(image_cl, width=unit(1,"npc"), height=unit(1,"npc"))
  # images_tibble <- tibble(sample=sample_names[i], grob=grobs)
  # images_tibble$height = dims$height
  # images_tibble$width = dims$width
  # images_tibble
  
  
  # ---------------------------
  # create SingleCellExperiment
  # ---------------------------
  
  # note: check and/or re-order rows to make sure barcode IDs match in df_barcodes and df_tisspos
  ord <- match(df_barcodes$barcode_id, df_tisspos$barcode)
  df_tisspos_ord <- df_tisspos[ord, ]
  dim(df_tisspos_ord)
  stopifnot(nrow(df_barcodes) == nrow(df_tisspos_ord))
  stopifnot(all(df_barcodes$barcode_id == df_tisspos_ord$barcode))

  col_data <- cbind(df_barcodes, df_tisspos_ord[, -1])
  head(col_data)
  
  sce <- SingleCellExperiment(
    rowData = df_features, 
    colData = col_data, 
    assays = c(counts = counts), 
    metadata = list(scalefactors = scalefactors)
  )
  
  sce
  
  #sum UMIs
  sce$sum_umi <- colSums(counts)
  sce$sum_gene = colSums(counts > 0)
  
  #adding colData
  sce$sample_name <- sample_names[i]
  #sce$sample_name <- sample_names[3]
  
  # store object
  sce_list[[i]] <- sce
  #sce_list[[3]] <- sce
  
  rm(sce)
}

sce_list
save(sce_list, file = "~/DLPFC/sce_list.rda")
load(file = "~/DLPFC/sce_list.rda")

image_paths <- paste0(dir_outputs,"/",sample_names, "/outs/spatial/tissue_lowres_image.png")
images_cl <- lapply(image_paths, read.bitmap) 
dims = t(sapply(images_cl, dim))
colnames(dims) = c("height", "width", "channel")
dims = as.data.frame(dims)


## ------------------------------------------------------------------------
grobs <- lapply(images_cl, rasterGrob, width=unit(1,"npc"), height=unit(1,"npc"))
images_tibble <- tibble(sample=sample_names, grob=grobs)
images_tibble$height = dims$height
images_tibble$width = dims$width
images_tibble

# sample                             grob       height width
# <chr>                              <list>      <int> <int>
#   1 DLPFC_Br2743_ant_manual_alignment  <rastrgrb>    600   504
# 2 DLPFC_Br2743_mid_manual_alignment  <rastrgrb>    600   517
# 3 DLPFC_Br2743_post_manual_alignment <rastrgrb>    600   452
# 4 DLPFC_Br3942_ant_manual_alignment  <rastrgrb>    600   504
# 5 DLPFC_Br3942_mid_manual_alignment  <rastrgrb>    600   517
# 6 DLPFC_Br3942_post_manual_alignment <rastrgrb>    600   452
# 7 DLPFC_Br6423_ant_manual_alignment  <rastrgrb>    600   504
# 8 DLPFC_Br6423_mid_manual_alignment  <rastrgrb>    600   517
# 9 DLPFC_Br6423_post_manual_alignment <rastrgrb>    600   405
# 10 DLPFC_Br8492_ant_manual_alignment  <rastrgrb>    600   504
# 11 DLPFC_Br8492_mid_manual_alignment  <rastrgrb>    600   517
# 12 DLPFC_Br8492_post_manual_alignment <rastrgrb>    600   561

##combine to one SCE object 
sce <- do.call(cbind, sce_list)

##add image data to meta data
metadata(sce)$image <- images_tibble

#save
save(sce, file = "~/DLPFC/sce_combined.rda")
load(file = "~/DLPFC/sce_combined.rda")


#quality control (scran) start here 1/25/21
qcstats <- perCellQCMetrics(sce)
qcfilter <- quickPerCellQC(qcstats)
colSums(as.matrix(qcfilter))

sce$scran_discard <- factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
sce$scran_low_lib_size <- factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
sce$low_n_features <- factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))



for (i in seq_along(sample_names)) {
  
  # select sample
  sce <- sce_list[[i]]
  
  # identify mitochondrial genes
  is_mito <- grepl("(^MT-)|(^mt-)", rowData(sce)$gene_name)
  table(is_mito)
  rowData(sce)$gene_name[is_mito]
  
  # calculate QC metrics using scater package
  sce <- addPerCellQC(sce, subsets = list(mito = is_mito))
  
  colData(sce)
  
  # store
  sce_list[[i]] <- sce
}

human_markers <- c("SNAP25", "MBP","PCP4", "RELN","AQP4","CUX2","CCK","HPCAL1", "NR4A2","RORB")

colors <- c("navy", "dodgerblue2")

pdf('DLPFC/marker_genes.pdf', useDingbats = FALSE)

for (i in seq_along(sample_names)) {
  
  # select sample
  sce <- sce_list[[i]]
  
  for (j in seq_along(human_markers)) {
    
    # identify marker gene
    ix_marker <- which(toupper(rowData(sce)$gene_name) == toupper(human_markers[j]))
    stopifnot(length(ix_marker) == 1)
    colData(sce)$counts_marker <- counts(sce)[ix_marker, ]
    
    
    # plot UMI counts for marker gene
    
    p <- ggplot(as.data.frame(colData(sce)), 
                aes(x = pxl_row_in_fullres, y = pxl_col_in_fullres, color = counts_marker)) + 
      geom_point(size = 1.0) + 
      coord_fixed() + 
      scale_y_reverse() + 
      scale_color_gradient(low = "gray95", high = colors[1]) + 
      ggtitle(paste0("UMI counts: ", human_markers[j], ": ", sample_names[i])) + 
      labs(color = "counts") + 
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    print(p)
  }
}
dev.off()


# human_markers <- c("SNAP25", "MBP","PCP4", "RELN","AQP4","CUX2","CCK","HPCAL1")
# 
# colors <- c("navy", "dodgerblue2")
pdf('DLPFC/marker_genes_by_gene.pdf', useDingbats = FALSE)
for (i in seq_along(human_markers)) {
  
  # select sample
  
  
  for (j in seq_along(sample_names)) {
    sce <- sce_list[[j]]
    
    # identify marker gene
    ix_marker <- which(toupper(rowData(sce)$gene_name) == toupper(human_markers[i]))
    stopifnot(length(ix_marker) == 1)
    colData(sce)$counts_marker <- counts(sce)[ix_marker, ]
    
    
    # plot UMI counts for marker gene
    
    p <- ggplot(as.data.frame(colData(sce)), 
                aes(x = pxl_row_in_fullres, y = pxl_col_in_fullres, color = counts_marker)) + 
      geom_point(size = 1.0) + 
      coord_fixed() + 
      scale_y_reverse() + 
      scale_color_gradient(low = "gray95", high = colors[1]) + 
      ggtitle(paste0("UMI counts: ", human_markers[i], ": ", sample_names[j])) + 
      labs(color = "counts") + 
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    print(p)
  }
}
dev.off()
