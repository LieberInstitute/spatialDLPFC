function refineSegmentation(fname)
load([fname(1:end-4),'_nuclei_WS.mat'])
he = imread(fname);

BW = bwareafilt(mask_dark_blue,[30 4000]); 
stats =  struct2table(regionprops(BW, rgb2gray(he), 'Area', 'BoundingBox', 'Centroid', 'Circularity', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Perimeter', 'MeanIntensity' , 'WeightedCentroid'));

% idx = find([stats.Eccentricity]<0.9);
% BW2 = ismember(bwlabel(BW),idx);
% stats =  struct2table(regionprops(BW2, rgb2gray(he), 'Area', 'BoundingBox', 'Centroid', 'Circularity', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Perimeter', 'MeanIntensity' , 'WeightedCentroid'));

mask_dark_blue = BW;
save([fname(1:end-4),'_nuclei_WS_final.mat'],'mask_dark_blue','-v7.3')
blue_nuclei = he .* uint8(mask_dark_blue);
imwrite(blue_nuclei,[fname(1:end-4),'_nuclei_WS_final.tif'])
writetable(stats,[fname(1:end-4),'_refineVNS_metric.csv'])

end

