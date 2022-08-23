import os,sys,csv,re

#  Some of our images exceed the default maximum number of pixels
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(pow(2,40))

import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import cv2
from scanpy import read_10x_h5
from anndata import AnnData
import pyhere
import getopt

analysis_name = '04-spatial_dlpfc'

out_dir_processed = pyhere.here('spagcn', 'processed-data', analysis_name)
out_dir_plots = pyhere.here('spagcn', 'plots', analysis_name)
input_adata_path = pyhere.here(out_dir_processed, 'spe_anndata.h5ad')
NUM_META_COLUMNS = 5
N_CLUSTERS = 7

#  Column names expected to be present in adata.obs for X and Y array and pixel
#  coordinates, respectively
array_names = ("array_row", "array_col")
pixel_names = ("pxl_col_in_fullres", "pxl_row_in_fullres")

###############################################################################
#  Helper functions
###############################################################################

#  Given an AnnData and Ensembl ID, return the gene symbol
def get_symbol(anndata, ens_id):
    return anndata.var.gene_name[anndata.var.gene_id == ens_id][0]

#  Convert meta-gene name string from 'spg.find_meta_gene' to string
#  containing gene symbol(s) instead of Ensembl IDs. Also correct name,
#  accounting for a bug: a gene name may appear up to 3 times
#  (e.g. 'gene1+gene2-gene1-gene1' should really be 'gene2-gene1')
def parse_metaname(adata, meta_name, convert=True):
    #  Get lists of genes whose expressions are high and low with respect to
    #  outside the spatial domain
    additive = [x.split('-')[0] for x in meta_name.split('+')]
    subtractive = [x.split('+')[0] for x in meta_name.split('-')][1:]
    
    #  Due to (arguably) a bug, it's possible for a gene to be added then
    #  subtracted. Remove the first instance of any "duplicate" gene
    overlaps = [x for x in additive if x in subtractive]
    for x in list(np.unique(overlaps)):
        additive.remove(x)
        subtractive.remove(x)
    
    if convert:
        #  Convert from Ensembl IDs to gene symbols
        additive = [get_symbol(adata, x) for x in additive]
        subtractive = [get_symbol(adata, x) for x in subtractive]
    
    #  Form a new "meta_name" string of the same format as the input
    if len(subtractive) > 0:
        new_meta_name = " + ".join(additive) + " - " + " - ".join(subtractive)
    else:
        new_meta_name = " + ".join(additive)
    
    return new_meta_name

#  Given an AnnData, a particular number of clusters, and other info, return an
#  output AnnData containing additional 'obs' columns with info about cluster
#  predictions
def find_domains(adata_in, x_array, y_array, adj, l, n_clusters, seed):
    adata = adata_in.copy()
    
    #Search for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=200, r_seed=seed, t_seed=seed, n_seed=seed)
    
    #  Train the model; a while loop is used because the training process is
    #  nondeterministic (even with the seeds set) and not guaranteed on a
    #  given run to find the desired number of clusters
    found_clusters = 0
    while found_clusters != n_clusters:
        #Run
        clf=spg.SpaGCN()
        clf.set_l(l)
        clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
        y_pred, prob=clf.predict()
        
        found_clusters = len(np.unique(y_pred))
    
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    
    #Do cluster refinement(optional)
    #shape="hexagon" for Visium data, "square" for ST data.
    adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    
    return adata

#  Given an AnnData, a column name in 'adata.obs', a color list, plot title, 
#  and output file path, create a PNG file containing a plot
def plot_adata(adata, colname, title, plot_color, out_file):
    ax = sc.pl.scatter(
        adata, alpha=1, x="y_pixel", y="x_pixel", color=colname, title=title,
        color_map=plot_color, show=False, size=100000/adata.shape[0]
    )
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    plt.savefig(out_file, dpi=600)
    plt.close()

#  Given an AnnData, 'meta_name' as output as from the first index of
#  'spg.find_meta_gene', data frame containing cluster indices by spot,
#  domain index 'target', and number of meta-gene columns to produce, this
#  function modifies 'cluster_list' to replace NAs present in the meta-gene
#  columns with gene names for spots belonging to the given domain.
def fill_meta_column(adata, meta_name, cluster_list, target, NUM_META_COLUMNS):
    #  Form a list of meta genes
    meta_list = parse_metaname(adata, meta_name, convert=False)
      
    if (meta_list[0] == ' '): # really asking if the first meta gene is subtracted
        meta_list = meta_list[1:]
        assert meta_list[0] == '-'
    else:
        meta_list = '+' + meta_list
    
    meta_list = meta_list.replace(' + ', ' +').replace(' - ', ' -').split(' ')
    assert all([x[0] in ['+', '-'] for x in meta_list])
    
    if len(meta_list) > NUM_META_COLUMNS:
        print('Warning: dropping one or more meta genes for this domain when forming CSV, since we are only taking the first ' + str(NUM_META_COLUMNS) + '.')
    
    #  Populate each meta gene column for rows associated with this
    #  cluster
    for i in range(min(NUM_META_COLUMNS, len(meta_list))):
        cluster_list['meta_gene_' + str(i + 1)][cluster_list.raw_cluster == target] = meta_list[i]
    

#  Given a raw Anndata, processed AnnData, X and Y array coordinates, and a
#  target spatial domain, return a DataFrame of info about SVGs. Under special
#  conditions, None may be returned instead of a DataFrame, if no valid SVGs
#  can be found.
def get_svgs(raw, adata, x_array, y_array, target):
    #Set filtering criterials
    min_in_group_fraction=0.8
    min_in_out_group_ratio=1
    min_fold_change=1.5
    
    #Search radius such that each spot in the target domain has approximately 10 neighbors on average
    adj_2d=spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    start, end= np.quantile(adj_2d[adj_2d!=0],q=0.001), np.quantile(adj_2d[adj_2d!=0],q=0.1)
    r=spg.search_radius(target_cluster=target, cell_id=adata.obs.index.tolist(), x=x_array, y=y_array, pred=adata.obs["pred"].tolist(), start=start, end=end, num_min=10, num_max=14,  max_run=100)
    
    #Detect neighboring domains
    nbr_domians=spg.find_neighbor_clusters(target_cluster=target,
                                       cell_id=raw.obs.index.tolist(), 
                                       x=raw.obs["x_array"].tolist(), 
                                       y=raw.obs["y_array"].tolist(), 
                                       pred=raw.obs["pred"].tolist(),
                                       radius=r,
                                       ratio=1/2)
    
    #  I've observed that domains where no neighbors are found are spatially
    #  scattered and are likely not biologically meaningful or interpretable.
    #  Stop here in that case, and don't find SVGs                               
    if nbr_domians is None:
        print('Found no neighbors for domain', target, '. No SVGs will be found for this domain.')
        return None
    
    nbr_domians=nbr_domians[0:3]
    de_genes_info=spg.rank_genes_groups(input_adata=raw,
                                    target_cluster=target,
                                    nbr_list=nbr_domians, 
                                    label_col="pred", 
                                    adj_nbr=True, 
                                    log=True)
    
    #  Take significant genes only                                
    de_genes_info=de_genes_info[(de_genes_info["pvals_adj"]<0.05)]
    if de_genes_info.shape[0] == 0:
        print("No significant SVGs found for domain " + str(target) + "!")
        return None
    
    filtered_info=de_genes_info
    
    #  Take genes that exceed the expression seen in neighboring domains
    if np.count_nonzero(filtered_info["in_out_group_ratio"] > min_in_out_group_ratio) == 0:
        print("Can't find any genes such that more spots in domain " + str(target) + " express the gene than its neighbors!")
        return None
        
    filtered_info = filtered_info[filtered_info["in_out_group_ratio"] > min_in_out_group_ratio]
    
    #  Relax filtering criteria if required to find at least one SVG
    while np.count_nonzero(
        (filtered_info["in_group_fraction"] > min_in_group_fraction) &
        (filtered_info["fold_change"] > min_fold_change)
        ) == 0:
        #  Relax both parameters simultaneously
        min_in_group_fraction -= 0.1
        min_fold_change = (min_fold_change + 1) / 2
        
        if min_in_group_fraction <= 0:
            print('No "min_in_group_fraction" exists that finds any SVGs for domain ' + str(target) + "!")
            return None
        
        if min_fold_change <= 1.01:
            print('No "fold_change" > 1.01 exists that finds any SVGs for domain ' + str(target) + "!")
            return None
        
        print('Warning: lowering "min_in_group_fraction" to ', min_in_group_fraction, 'and "min_fold_change" to ', min_fold_change, 'for domain', target, 'to find at least 1 SVG.')
    
    #  Filter significant genes
    filtered_info=filtered_info[(filtered_info["in_group_fraction"] > min_in_group_fraction) &
                                (filtered_info["fold_change"] > min_fold_change)]
    
    filtered_info=filtered_info.sort_values(by="in_group_fraction", ascending=False)
    filtered_info["target_dmain"]=target
    filtered_info["neighbors"]=str(nbr_domians)
    
    return filtered_info
    
###############################################################################
#  Main Analysis
###############################################################################

#  Set seed
seed = 100
random.seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

#  Receive command-line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:], "i:s:", ["image_path=", "sample_name="])
except getopt.GetoptError:
    print('python 02-run_spagcn.py -i <path to image> -s <sample ID>')
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-i', '--image_path='):
        image_path = arg
    elif opt in ('-s', '--sample_name='):
        sample_name = arg
    else:
        print('Unknown argument.')
        print('python 02-run_spagcn.py -i <path to image> -s <sample ID>')
        sys.exit(2)

adata = sc.read_h5ad(input_adata_path)

#  Subset to just this sample
adata = adata[adata.obs.sample_id == sample_name, :]

print('Continuing with sample ' + sample_name + '...')
out_dir_processed = pyhere.here(out_dir_processed, sample_name)
out_dir_plots = pyhere.here(out_dir_plots, sample_name)
if not os.path.exists(out_dir_processed):
    os.mkdir(out_dir_processed)
if not os.path.exists(out_dir_plots):
    os.mkdir(out_dir_plots)

#  Rename and set some variables for compatibility with original code from
#  tutorial
adata.obs["x_array"]=adata.obs[array_names[0]]
adata.obs["y_array"]=adata.obs[array_names[1]]

adata.obs["x_pixel"]= adata.obs[pixel_names[0]]
adata.obs["y_pixel"]= adata.obs[pixel_names[1]]

x_array=adata.obs["x_array"].astype(np.int32).tolist()
y_array=adata.obs["y_array"].astype(np.int32).tolist()
x_pixel=adata.obs["x_pixel"].astype(np.int32).tolist()
y_pixel=adata.obs["y_pixel"].astype(np.int32).tolist()

#  Read in histology image
img=cv2.imread(image_path)

#Test coordinates on the image
img_new=img.copy()
for i in range(len(x_pixel)):
    x=x_pixel[i]
    y=y_pixel[i]
    img_new[int(x-20):int(x+20), int(y-20):int(y+20),:]=0

cv2.imwrite(
    str(pyhere.here(out_dir_plots, "coord_test.jpg")),
    img_new
)

#Calculate adjacent matrix

#  Are the values for 's' and 'b' reasonable?
s=1
b=49

adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
np.savetxt(
    pyhere.here(out_dir_processed, "adj.csv"),
    adj,
    delimiter=','
)

adata.var_names_make_unique()
spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)

#Normalize and take log for UMI
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

#  Is the value for 'p' reasonable?

p=0.5 
#Find the l value given p
l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    
#  Compute/ determine spatial domains
adata = find_domains(adata, x_array, y_array, adj, 1, N_CLUSTERS, seed)

#  With Louvain clustering, the exact 'n_clusters' specified might not be
#  found.
actual_n_clusters = len(adata.obs.pred.unique())
refined_n_clusters = len(adata.obs.refined_pred.unique())
if actual_n_clusters != refined_n_clusters:
    print(
        'Warning:', actual_n_clusters,
        'raw clusters were found, but this changed to' ,refined_n_clusters,
        'after refining!'
    )

this_out_dir_plots = pyhere.here(out_dir_plots, str(actual_n_clusters) + "_clusters")
this_out_dir_processed = pyhere.here(out_dir_processed, str(actual_n_clusters) + "_clusters")
if not os.path.exists(this_out_dir_plots):
    os.mkdir(this_out_dir_plots)
if not os.path.exists(this_out_dir_processed):
    os.mkdir(this_out_dir_processed)

#  Save result AnnData
out_file = pyhere.here(this_out_dir_processed, 'results.h5ad')
adata.write_h5ad(out_file)

#  Form and write a CSV file containing rows of the form:
#  (barcode + sample ID, raw cluster num, refined cluster num)
#  This is designed for compatibility with the R function
#  spatialLIBD::import_cluster
cluster_list = adata.obs[['pred', 'refined_pred']]
cluster_list.index += '_' + sample_name
cluster_list.index.name = 'key'
cluster_list.rename(
    columns={'pred': 'raw_cluster', 'refined_pred': 'refined_cluster'},
    inplace=True
)

#  Plot spatial domains, raw and refined
plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]

adata.uns["pred_colors"]=list(plot_color[:actual_n_clusters])
out_file = pyhere.here(this_out_dir_plots, "domains.png")
plot_adata(adata, "pred", "Raw spatial domains", plot_color, out_file)

adata.uns["refined_pred_colors"]=list(plot_color[:refined_n_clusters])
out_file = pyhere.here(this_out_dir_plots, "refined_domains.png")
plot_adata(adata, "refined_pred", "Refined spatial domains", plot_color, out_file)

#  Read in raw data and subset to current sample
raw = sc.read_h5ad(input_adata_path)
raw.var_names_make_unique()
raw = raw[raw.obs.sample_id == sample_name, :]

raw.obs["pred"] = adata.obs["pred"].astype('category') # The 'adata' instead of 'raw' here is intentional
raw.obs["x_array"]=raw.obs[array_names[0]]
raw.obs["y_array"]=raw.obs[array_names[1]]
raw.obs["x_pixel"]= raw.obs[pixel_names[0]]
raw.obs["y_pixel"]= raw.obs[pixel_names[1]]

#Convert sparse matrix to non-sparse
raw.X=(raw.X.A if issparse(raw.X) else raw.X)
raw.raw=raw
sc.pp.log1p(raw)

assert all(cluster_list.raw_cluster.values == raw.obs.pred.values)

#  Initialize empty meta-gene columns
for i in range(NUM_META_COLUMNS):
    cluster_list['meta_gene_' + str(i + 1)] = ''
    
#  Find meta genes for each domain found (due to a bug, domain names are not
#  necessarily sequential, hence why 'range' is not used here)
for target in adata.obs.pred.cat.categories:
    #  Determine SVGs
    filtered_info = get_svgs(raw, adata, x_array, y_array, target)
    if filtered_info is None:
        continue
    
    print("SVGs for domain ", str(target),":", filtered_info["genes"].tolist())
    
    #  Order SVGs by "in_group_fraction" and within that, "pvals_adj"
    filtered_info=filtered_info.sort_values(by="pvals_adj", ascending=True)
    filtered_info=filtered_info.sort_values(by="in_group_fraction", ascending=False)
    start_gene = filtered_info.genes.values[0]
    
    meta_name, meta_exp=spg.find_meta_gene(input_adata=raw,
                        pred=raw.obs["pred"].tolist(),
                        target_domain=target,
                        start_gene=start_gene,
                        mean_diff=0,
                        early_stop=True,
                        max_iter=3,
                        use_raw=False)
    
    raw.obs["meta"]=meta_exp
    
    #  Fill in meta-gene columns in 'cluster_list' that are associated with
    #  domain 'target'
    fill_meta_column(raw, meta_name, cluster_list, target, NUM_META_COLUMNS)
    
    color_self = clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
    
    #  Plot "start gene" in meta gene set
    title = get_symbol(raw, start_gene) + ' (marker for domain ' + str(target) + ')'
    out_file = pyhere.here(this_out_dir_plots, "start_meta_gene_domain_" + str(target) + ".png")
    raw.obs["exp"]=raw.X[:,raw.var.index==start_gene]
    plot_adata(raw, "exp", title, color_self, out_file)
    
    #  Plot entire meta gene set
    raw.obs["exp"]=raw.obs["meta"]
    out_file = pyhere.here(this_out_dir_plots, "all_meta_genes_domain_" + str(target) + ".png")
    plot_adata(raw, "exp", parse_metaname(raw, meta_name), color_self, out_file)


out_file = pyhere.here(this_out_dir_processed, 'clusters.csv')
cluster_list.to_csv(out_file)
