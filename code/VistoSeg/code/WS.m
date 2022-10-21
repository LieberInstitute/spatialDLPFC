function WS(fname,M)

load([fname(1:end-4),'_mask.mat'])
mask_dark_blue = mask{M};
he = imread(fname);

D = -bwdist(~mask_dark_blue);
mask = imextendedmin(D,1);
% figure
% imshowpair(mask_dark_blue,mask,'blend')

D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = mask_dark_blue;
bw3(Ld2 == 0) = 0;
%imshow(bw3)

stats =  struct2table(regionprops(bw3, rgb2gray(he), 'Area', 'BoundingBox', 'Centroid', 'Circularity', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Perimeter', 'MeanIntensity' , 'WeightedCentroid'));
writetable(stats,[fname(1:end-4),'_refineVNS_metric.csv'])
mask_dark_blue = bw3;
save([fname(1:end-4),'_nuclei_WS.mat'],'mask_dark_blue','-v7.3')
blue_nuclei = he .* uint8(mask_dark_blue);
imwrite(blue_nuclei,[fname(1:end-4),'_nuclei.tif'])
end

