## Automatically style the code in this script:
styler::style_file(
    here::here("analysis", "02_marker_genes.R"),
    transformers = biocthis::bioc_style()
)

## This script requires R 4.1
# module load conda_R/devel

## utils
library("here")
library("sessioninfo")

## reading the data
library("SpatialExperiment")

##analysis
library("scran")

## vis
library("spatialLIBD")
library("RColorBrewer")
nb.cols <- 25
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

## Load SPE data
load(here::here("processed-data", "rdata", "spe", "spe_090821.Rdata"), verbose = TRUE)

## Find marker genes
human_markers <-
    c(
        "SNAP25",
        "MBP",
        "PCP4",
        "RELN",
        "AQP4",
        "CUX2",
        "CCK",
        "HPCAL1",
        "NR4A2",
        "RORB"
    )

## Locate the marker genes
human_markers_search <- rowData(spe)$gene_search[match(human_markers, rowData(spe)$gene_name)]

## Create plots directory
dir.create(here::here("plots", "human_markers"), showWarnings = FALSE)

for (i in human_markers_search) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "human_markers", paste0(gsub("; ", "_", i), ".pdf")),
        assayname = "counts"
    )
}

sample_ids <- unique(colData(spe)$sample_id)




#To see which clusters are present in each sample
for (i in seq_along(sample_ids)){
    print(sample_ids[i])
    print(table(colData(spe)$SNN_k50_k4[which(colData(spe)$sample_id==sample_ids[i])]))
}



cluster_colNames <- c(paste0("SNN_k50_k", seq_len(28)))
cluster_colNames <-cluster_colNames[4:28]

#plot all samples (12) for all clustervars (25)
pdf(file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/plots/vis_grid_clus_spe.pdf")
for(i in seq_along(cluster_colNames)){
    for(j in seq_along(sample_ids)){
            if (enough_ram()) {
                ## Obtain the necessary data
                if (!exists("spe")) spe <- fetch_data("spe")
                
                ## Subset to two samples of interest
                #spe_sub <- spe[, spe$sample_id %in% "Br2743_ant"]
                
                
                ## Obtain the plot list
                p_list <-
                    vis_grid_clus(
                        spe[,spe$sample_id == sample_ids[j]],
                        cluster_colNames[i],
                        spatial = TRUE,
                        return_plots = TRUE,
                        sort_clust = TRUE,
                        colors =  c("#b2df8a", "#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "gold", "#a65628",
                                    "#999999", "black", "grey", "white", "purple"),
                        height = 24,
                        width = 60
                        
                    )
                
                ## Clean up
                rm(spe)
                
                ## Visualize the spatial adjacent replicates for position = 0 micro meters
                ## for subject 3
                print(cowplot::plot_grid(plotlist = p_list, ncol = 1)) #print around here??
            }
    }
}

dev.off()

pdf(file = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/plots/vis_clus_spe.pdf")
for (i in seq_along(sample_ids)){
    for(j in seq_along(cluster_colNames)){
        my_plot <- vis_clus(
            spe = spe,
            clustervar = cluster_colNames[j],
            sampleid = sample_ids[i],
            colors =  mycolors,
            ... = paste0(" ",cluster_colNames[j])
        )
        print(my_plot)
    }
    
}

dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

