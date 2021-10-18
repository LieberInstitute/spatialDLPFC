#  Based on:
#  https://github.com/jianhuupenn/SpaGCN/blob/8fb25b2ec918702b1afb0a80662ae84939f87d24/tutorial/tutorial.md

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
import pyhere # this was added; analagous to 'here' in R

first_sample_name = '151507'


os.mkdir(pyhere.here("spagcn/processed-data/02-our_data/sample_results"))

adata = sc.read_h5ad(
    pyhere.here("spagcn/processed-data/02-our_data/spe_anndata.h5ad")
)

#  Take just the first sample for testing
adata = adata[adata.obs.sample_id == first_sample_name, :]

#  Rename and set some variables for consistency with the remaining tutorial
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
    str(pyhere.here("spagcn/raw-data/02-our_data/" + first_sample_name + ".tif"))
)

#Test coordinates on the image
img_new=img.copy()
for i in range(len(x_pixel)):
    x=x_pixel[i]
    y=y_pixel[i]
    img_new[int(x-20):int(x+20), int(y-20):int(y+20),:]=0

cv2.imwrite(
    str(pyhere.here("spagcn/processed-data/02-our_data/coord_test.jpg")),
    img_new
)

#Calculate adjacent matrix

#  Are the values for 's' and 'b' reasonable?
s=1
b=49

adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
np.savetxt(
    pyhere.here("spagcn/processed-data/02-our_data/adj.csv"),
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
    pyhere.here("spagcn/processed-data/02-our_data/results.h5ad")
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
    pyhere.here("spagcn/processed-data/02-our_data/pred.png"), dpi=600
)
plt.close()

domains="refined_pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig(
    pyhere.here("spagcn/processed-data/02-our_data/refined_pred.png"),
    dpi=600
)
plt.close()

#Read in raw data
raw = sc.read_h5ad(
    pyhere.here("spagcn/processed-data/02-our_data/spe_anndata.h5ad")
)
raw.var_names_make_unique()

#  Take just the first sample for testing
raw = raw[raw.obs.sample_id == first_sample_name, :]

raw.obs["pred"] = adata.obs["pred"].astype('category') # The 'adata' instead of 'raw' here is intentional
raw.obs["x_array"]=raw.obs["array_row"]
raw.obs["y_array"]=raw.obs["array_col"]
raw.obs["x_pixel"]= raw.obs["pxl_col_in_fullres"]
raw.obs["y_pixel"]= raw.obs["pxl_row_in_fullres"]

#Convert sparse matrix to non-sparse
raw.X=(raw.X.A if issparse(raw.X) else raw.X)
raw.raw=raw
sc.pp.log1p(raw)

#Use domain 0 as an example
target=0

#Set filtering criterials (are these appropriate for our data?)
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

nbr_domians=nbr_domians[0:3]
de_genes_info=spg.rank_genes_groups(input_adata=raw,
                                target_cluster=target,
                                nbr_list=nbr_domians, 
                                label_col="pred", 
                                adj_nbr=True, 
                                log=True)
#Filter genes
de_genes_info=de_genes_info[(de_genes_info["pvals_adj"]<0.05)]
filtered_info=de_genes_info
filtered_info=filtered_info[(filtered_info["pvals_adj"]<0.05) &
                            (filtered_info["in_out_group_ratio"]>min_in_out_group_ratio) &
                            (filtered_info["in_group_fraction"]>min_in_group_fraction) &
                            (filtered_info["fold_change"]>min_fold_change)]
filtered_info=filtered_info.sort_values(by="in_group_fraction", ascending=False)
filtered_info["target_dmain"]=target
filtered_info["neighbors"]=str(nbr_domians)
print("SVGs for domain ", str(target),":", filtered_info["genes"].tolist())

color_self = clr.LinearSegmentedColormap.from_list('pink_green', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
for g in filtered_info["genes"].tolist():
    raw.obs["exp"]=raw.X[:,raw.var.index==g]
    ax=sc.pl.scatter(raw,alpha=1,x="y_pixel",y="x_pixel",color="exp",title=g,color_map=color_self,show=False,size=100000/raw.shape[0])
    ax.set_aspect('equal', 'box')
    ax.axes.invert_yaxis()
    plt.savefig(
        pyhere.here(
            "spagcn/processed-data/02-our_data/sample_results/" + g + ".png"
        ), dpi=600
    )
    plt.close()
