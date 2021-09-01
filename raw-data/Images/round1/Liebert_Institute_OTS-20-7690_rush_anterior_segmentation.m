%%%%%%%%%%%% IMG1 %%%%%%%%%%%%%%%
Img1 = imread('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_1.tif');
Img1_smooth = imgaussfilt(Img1,4);
Img1_smooth_adj = imadjust(Img1_smooth, [.2 .3 0; .6 .7 1],[]);

he = Img1;
lab_he = rgb2lab(he);
ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 5;
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3);

parfor i = 1:nColors
mask{i} = pixel_labels==i;
cluster{i} = he .* uint8(mask{i});
end

%I = 1;imshow(mask(I))%to check the matching mask

	for i = 1:5
nuclei_mask = mask{i};
imwrite(cluster{i},['/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_1_cluster',num2str(i),'.tif']) 
save(['/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_1_segmentation',num2str(i),'.mat'],'nuclei_mask')
    end
%% IMG2 %%
img2 = imread('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_2.tif');
img2_smooth = imgaussfilt(img2,4);
img2_smooth_adj = imadjust(img2_smooth, [.2 .3 0; .6 .7 1],[]);


he = img2_smooth_adj;
lab_he = rgb2lab(he);
ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 5;
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3);

parfor i = 1:nColors
mask{i} = pixel_labels==i;
cluster{i} = he .* uint8(mask{i});
end

nuclei_mask = mask{4};
imwrite(cluster{4},'/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_2_cluster4Nuclei.tif') 
save('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_2_nucleisegmentation.mat','nuclei_mask')

%% IMG3 %%

Img3 = imread('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_anterior_3.tif');
Img3_smooth = imgaussfilt(Img3,4);
Img3_smooth_adj = imadjust(Img3_smooth, [.2 .3 0; .6 .7 1],[]);

he = Img3_smooth_adj;
lab_he = rgb2lab(he);
ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 5;
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3);

parfor i = 1:nColors
mask{i} = pixel_labels==i;
cluster{i} = he .* uint8(mask{i});
end

nuclei_mask = mask{4};
imwrite(cluster{4},'/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_anterior_3_cluster4Nuclei.tif') 
save('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_anterior_3_nucleisegmentation.mat','nuclei_mask')

%% IMG4 %%

img4 = imread('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Liebert_Institute_OTS-20-7748_rush_anterior_3.tif');
img4_smooth = imgaussfilt(img4,4);
img4_smooth_adj = imadjust(img4_smooth, [.2 .3 0; .6 .7 1],[]);

he = img4_smooth_adj;
lab_he = rgb2lab(he);
ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 5;
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3);

parfor i = 1:nColors
mask{i} = pixel_labels==i;
cluster{i} = he .* uint8(mask{i});
end

nuclei_mask = mask{5};
imwrite(cluster{5},'/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_3_cluster5Nuclei.tif') 
save('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_3_nucleisegmentation.mat','nuclei_mask')
