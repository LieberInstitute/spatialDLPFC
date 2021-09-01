cd /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/spotspotcheck-master/

img = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_2.tif';
im = imread(img);
mask = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_2_nuclei.mat';
load(mask);
BW = mask_dark_blue;
jsonname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/outputs/NextSeq/DLPFC_Br3942_post_manual_alignment/outs/spatial/scalefactors_json.json';
posname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/outputs/NextSeq/DLPFC_Br3942_post_manual_alignment/outs/spatial/tissue_positions_list.csv';
[posPath,~] = fileparts(posname);

w = jsondecode(fileread(jsonname));
R = ceil(w.spot_diameter_fullres/2);
tbl = readtable(posname) ;
    
count = [];
    if size(tbl, 2) == 7
        count = table2array(tbl(:, 7));
    else
        tic
        disp('counting nuclei')
        count = countSpots(BW, R, tbl, posPath);
        toc
    end
