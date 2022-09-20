sample_ids = ["Br2720_Ant_IF", "Br6432_Ant_IF", "Br6522_Ant_IF", "Br8667_Post_IF"];
img_ids = ["V10B01-087_A1", "V10B01-087_B1", "V10B01-087_C1", "V10B01-087_D1"];

% img_ids(getenv('SGE_TASK_ID')
img_path = strcat("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/VisiumIF/VistoSeg/", img_ids(1), ".tif");

N = 5;

cd '/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/spot_deconvo/06-vistoseg/VistoSeg/code'

VNS(img_path, N)
