function refineVNS(fname,M)
disp('loading data')
tic
%Img1 = imread(fname);
mask_name = [fname(1:end-4),'_mask.mat'];
load(mask_name)
toc

disp('refining segmentations')
tic
he = imread(fname);
lab_he = rgb2lab(he);
L_blue = lab_he(:,:,1) .* double(mask{M});
L_blue = rescale(L_blue);
idx_light_blue = imbinarize(nonzeros(L_blue));

blue_idx = find(mask{M});
mask_dark_blue = mask{M};
mask_dark_blue(blue_idx(idx_light_blue)) = 0;

blue_nuclei = he .* uint8(mask_dark_blue);
toc 

disp('saving final segmentations')
tic
imwrite(blue_nuclei,[fname(1:end-4),'_nuclei.tif'])
save([fname(1:end-4),'_nuclei.mat'],'mask_dark_blue','-v7.3')
toc
end


