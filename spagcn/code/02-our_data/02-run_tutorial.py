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

adata = sc.read_h5ad(
    pyhere.here("spagcn/processed-data/02-our_data/spe_anndata.h5ad")
)

#  Rename and set some variables for consistency with the remaining tutorial
adata.obs["x_array"]=adata.obs["array_row"]
adata.obs["y_array"]=adata.obs["array_col"]

adata.obs["x_pixel"]= adata.obs["pxl_row_in_fullres"]
adata.obs["y_pixel"]= adata.obs["pxl_col_in_fullres"]

x_array=adata.obs["x_array"].tolist()
y_array=adata.obs["y_array"].tolist()
x_pixel=adata.obs["x_pixel"].tolist()
y_pixel=adata.obs["y_pixel"].tolist()

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
