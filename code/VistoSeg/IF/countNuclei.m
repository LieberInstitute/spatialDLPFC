%mask = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_2_nuclei.mat';
%jsonname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/outputs/NextSeq/DLPFC_Br3942_post_manual_alignment/outs/spatial/scalefactors_json.json';
%posname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/outputs/NextSeq/DLPFC_Br3942_post_manual_alignment/outs/spatial/tissue_positions_list.csv';

function [count,prop] = countNuclei(mask,jsonname,posname) 

disp('loading data')

 BW = load(mask);
 O = fieldnames(BW);
   
[posPath,~] = fileparts(posname);

w = jsondecode(fileread(jsonname));
R = ceil(w.spot_diameter_fullres/2);
tbl = readtable(posname) ;
count = [];
prop = [];

    if size(tbl, 2) > 6
        b = 7;
        for C = 1:numel(O)
        count.(O{C}) = table2array(tbl(:, b));
        prop.(O{C}) = table2array(tbl(:, b+1));
        b = b+3;
        end
    else
        tbl.Properties.VariableNames = {'barcode','tissue','row','col','imagerow','imagecol'};
        [count,prop] = countSpots(BW, R, tbl, posPath);
        
    end
end

