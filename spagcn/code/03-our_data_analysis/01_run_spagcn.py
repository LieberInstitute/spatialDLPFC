import os,sys,csv,re
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

out_dir_processed = pyhere.here("spagcn/processed-data/03-our_data_analysis")
out_dir_plots = pyhere.here("spagcn/plots/03-our_data_analysis")
input_adata_path = pyhere.here("spagcn/processed-data/02-our_data_tutorial/spe_anndata.h5ad")
start_gene = "ENSG00000131095" # Finding meta-genes requires a "start gene"; this is GFAP here
NUM_META_COLUMNS = 5

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
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=seed, t_seed=seed, n_seed=seed)
    
    clf=spg.SpaGCN()
    clf.set_l(l)
    
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
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

###############################################################################
#  Main Analysis
###############################################################################

#  Recieve the '-i' argument, an integer in [1, 12] corresponding to the index
#  in the 12 samples present in the AnnData object
try:
    opts, args = getopt.getopt(sys.argv[1:], "i:", ["index="])
except getopt.GetoptError:
    print('python 03-main_analysis.py -i <sample index in 1-12>')
    sys.exit(2)

for opt, arg in opts:
    assert opt in ('-i', '--index='), opt
    sample_index = int(arg)

adata = sc.read_h5ad(input_adata_path)

#  Subset to just this sample
sample_name = adata.obs.sample_id.unique()[sample_index - 1]
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
adata.obs["x_array"]=adata.obs["array_row"]
adata.obs["y_array"]=adata.obs["array_col"]

adata.obs["x_pixel"]= adata.obs["pxl_col_in_fullres"]
adata.obs["y_pixel"]= adata.obs["pxl_row_in_fullres"]

x_array=adata.obs["x_array"].astype(np.int32).tolist()
y_array=adata.obs["y_array"].astype(np.int32).tolist()
x_pixel=adata.obs["x_pixel"].astype(np.int32).tolist()
y_pixel=adata.obs["y_pixel"].astype(np.int32).tolist()

#  Read in histology image
img=cv2.imread(
    str(pyhere.here("spagcn/raw-data/02-our_data_tutorial/" + sample_name + ".tif"))
)

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

#  Set seed
seed = 100
random.seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

adata_orig = adata

#  Try different numbers of clusters: 5, 7, 9 (might not be exact)
for n_clusters in range(5, 11, 2):
    adata = adata_orig.copy()
    
    #  Compute/ determine spatial domains
    adata = find_domains(adata, x_array, y_array, adj, 1, n_clusters, seed)
    
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
    cluster_list.index.name = 'Barcode'
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
    raw.obs["x_array"]=raw.obs["array_row"]
    raw.obs["y_array"]=raw.obs["array_col"]
    raw.obs["x_pixel"]= raw.obs["pxl_col_in_fullres"]
    raw.obs["y_pixel"]= raw.obs["pxl_row_in_fullres"]
    
    #Convert sparse matrix to non-sparse
    raw.X=(raw.X.A if issparse(raw.X) else raw.X)
    raw.raw=raw
    sc.pp.log1p(raw)
    
    cl_copy = cluster_list.copy()
    assert all(raw.obs.pred == cl_copy.raw_cluster)
    
    for i in range(NUM_META_COLUMNS):
        cl_copy['meta_gene_' + str(i + 1)] = ''
        
    #  Find meta genes for each domain found
    for target in range(actual_n_clusters):
        meta_name, meta_exp=spg.find_meta_gene(input_adata=raw,
                            pred=raw.obs["pred"].tolist(),
                            target_domain=target,
                            start_gene=start_gene,
                            mean_diff=0,
                            early_stop=True,
                            max_iter=3,
                            use_raw=False)
        
        raw.obs["meta"]=meta_exp
        
        #  Form a list of meta genes
        meta_list = parse_metaname(raw, meta_name, convert=False)
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
            cl_copy['meta_gene_' + str(i + 1)][cl_copy.raw_cluster == target] = meta_list[i]
        
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

    #  At this point, the meta gene columns should be fully populated for as
    #  many meta genes as exist for this domain
    for i in range(min(NUM_META_COLUMNS, len(meta_list))):
        assert '' not in cl_copy['meta_gene_' + str(i + 1)]
    
    out_file = pyhere.here(this_out_dir_processed, 'clusters.csv')
    cl_copy.to_csv(out_file)
