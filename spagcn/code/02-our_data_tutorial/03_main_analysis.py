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

adata = sc.read_h5ad(
    pyhere.here("spagcn/processed-data/02-our_data/spe_anndata.h5ad")
)

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
    str(pyhere.here("spagcn/raw-data/02-our_data/" + sample_name + ".tif"))
)

#Test coordinates on the image
img_new=img.copy()
for i in range(len(x_pixel)):
    x=x_pixel[i]
    y=y_pixel[i]
    img_new[int(x-20):int(x+20), int(y-20):int(y+20),:]=0

cv2.imwrite(
    str(pyhere.here("spagcn/plots/02-our_data/coord_test.jpg")),
    img_new
)
