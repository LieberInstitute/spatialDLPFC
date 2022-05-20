import numpy as np
from cellpose import models
from cellpose.io import imread
import pyhere
import pandas as pd

model_type='nuclei'
channel = 1 # DAPI
cell_diameter = None

#   Path to excel sheet containing sample info; directory containing .tif
#   Visium-IF images
sample_info_path = pyhere.here(
    'raw-data', 'sample_info', 'Visium_IF_DLPFC_MasterExcel_01262022.xlsx'
)
img_dir = pyhere.here('raw-data', 'Images', 'VisiumIF', 'VistoSeg')

#   Determine paths to IF images; read them in
sample_info = pd.read_excel(sample_info_path, header = 1)[:4]
sample_ids = sample_info['Slide SN #'] + '_' + sample_info['Array #']

img_paths = [str(x) + '.tif' for x in pyhere.here(img_dir, sample_ids)]
imgs = [imread(f)[channel, :, :] for f in img_paths]

#   Initialize the model and process images using it
model = models.Cellpose(gpu = True, model_type = model_type)
masks, flows, styles, diams = model.eval(imgs, diameter=cell_diameter)

#   TODO: save masks
