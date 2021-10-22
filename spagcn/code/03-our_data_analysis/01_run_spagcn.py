import os,csv,re
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
    #  subtracted twice. Remove the first instance of any "duplicate" gene
    overlaps = [x for x in additive if x in subtractive]
    for x in overlaps:
        additive.remove(x)
        subtractive.remove(x)
    
    if convert:
        #  Convert from Ensembl IDs to gene symbols
        additive = [get_symbol(adata, x) for x in additive]
        subtractive = [get_symbol(adata, x) for x in subtractive]
    
    #  Form a new "meta_name" string of the same format as the input
    if len(subtractive) > 0:
        new_meta_name = "+".join(additive) + "-" + "-".join(subtractive)
    else:
        new_meta_name = "+".join(additive)
    
    return new_meta_name

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
    str(pyhere.here(out_dir_plots, "coord_test_" + sample_name + ".jpg")),
    img_new
)

#Calculate adjacent matrix

#  Are the values for 's' and 'b' reasonable?
s=1
b=49

adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
np.savetxt(
    pyhere.here(out_dir_processed, "adj_" + sample_name + ".csv"),
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

#  Is this appropriate for our data?
n_clusters=7

#Set seed
r_seed=t_seed=n_seed=100

#Search for suitable resolution
res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

clf=spg.SpaGCN()
clf.set_l(l)

#Set seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)

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

#Save results
adata.write_h5ad(
    pyhere.here(out_dir_processed, "results_" + sample_name + ".h5ad")
)

#Set colors used
plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]

#Plot spatial domains
domains="pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig(
    pyhere.here(out_dir_plots, "domains_" + sample_name + ".png"), dpi=600
)
plt.close()

domains="refined_pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig(
    pyhere.here(out_dir_plots, "refined_domains_" + sample_name + ".png"),
    dpi=600
)
plt.close()

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

#Use domain 2 as an example
target=2
meta_name, meta_exp=spg.find_meta_gene(input_adata=raw,
                    pred=raw.obs["pred"].tolist(),
                    target_domain=target,
                    start_gene=start_gene,
                    mean_diff=0,
                    early_stop=True,
                    max_iter=3,
                    use_raw=False)

raw.obs["meta"]=meta_exp

#  Plot "start gene" in meta gene set
color_self = clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
g=start_gene
title = get_symbol(raw, g) + ' (marker for domain ' + str(target) + ')'
raw.obs["exp"]=raw.X[:,raw.var.index==g]
ax=sc.pl.scatter(raw,alpha=1,x="y_pixel",y="x_pixel",color="exp",title=title,color_map=color_self,show=False,size=100000/raw.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig(
    pyhere.here(out_dir_plots, "sample_results", "start_meta_gene_" + sample_name + "_domain_" + str(target) + ".png"),
    dpi=600
)
plt.close()

#  Plot entire meta gene set
raw.obs["exp"]=raw.obs["meta"]
ax=sc.pl.scatter(raw,alpha=1,x="y_pixel",y="x_pixel",color="exp",title=parse_metaname(raw, meta_name),color_map=color_self,show=False,size=100000/raw.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig(
    pyhere.here(out_dir_plots, "sample_results", "all_meta_genes_" + sample_name + "_domain_" + str(target) + ".png"),
    dpi=600
)
plt.close()
