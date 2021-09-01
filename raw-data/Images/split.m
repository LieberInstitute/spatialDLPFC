fname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior.tif';
numimgs = size(imfinfo(fname),1); 

parfor i = 1:numimgs
    tic
    I{i}.image = imread(fname,i);
    toc
disp(num2str(i))
end


save('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior.mat','I', '-v7.3');

% for i = 1:numimgs
% filename = ['/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush-001_mid_',num2str(i)];
% img = imresize(I{i}.image,0.5);
% imwrite(img,[filename,'.png']);
% imwrite(img,[filename,'.tif']);
% 
% x(i) = size(I{i}.image,1);
% y(i) = size(I{i}.image,2);
% end

Img = I{1}.image;

[x,y,z] = size(Img);

Img1 = Img(:,1:round(y/4),:);
Img2 = Img(:,round(y/4):round(y/4)*2,:);
Img3 = Img(:,75000:111000,:);
Img4 = Img(:,111000:end,:);

IMG1 = imresize(Img1,0.7);
IMG2 = imresize(Img2,0.7);
IMG3 = imresize(Img3,0.7);
IMG4 = imresize(Img4,0.7);

for i = 3:4
save(['/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_',num2str(i),'.mat'],['Img',num2str(i)], '-v7.3');
eval(['imwrite(IMG',num2str(i),',''/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush_posterior_',num2str(i),'.tif'')']);
end


%% image 2 %%%
fname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush-001_mid.tif';
numimgs = size(imfinfo(fname),1); 

parfor i = 1:numimgs
    tic
    I{i}.image = imread(fname,i);
    toc
disp(num2str(i))
end

save('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush-001_mid.mat','I', '-v7.3');

Img = I{1}.image;

[x,y,z] = size(Img);

Img1 = Img(:,1:round(y/4),:);
IMG1 = imresize(Img1,0.7);
Img2 = Img(:,round(y/4):round(y/4)*2,:);
IMG2 = imresize(Img2,0.7);
Img3 = Img(:,round(y/4)*2:round(y/4)*3,:);
IMG3 = imresize(Img3,0.7);
Img4 = Img(:,round(y/4)*3:end,:);
IMG4 = imresize(Img4,0.7);


for i = 1:4
save(['/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush-001_mid_',num2str(i),'.mat'],['Img',num2str(i)], '-v7.3');
eval(['imwrite(IMG',num2str(i),',''/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7748_rush-001_mid_',num2str(i),'.tif'')']);
end


%% image 3 %%%
fname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Lieber_Institute_OTS-20-7690_rush_anterior.tif';
numimgs = size(imfinfo(fname),1); 

parfor i = 1:numimgs
    tic
    I{i}.image = imread(fname,i);
    toc
disp(num2str(i))
end

save('/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Lieber_Institute_OTS-20-7690_rush_anterior.mat','I', '-v7.3');
Img = I{1}.image;

[x,y,z] = size(Img);

Img1 = Img(:,1:round(y/4),:);
Img2 = Img(:,round(y/4):round(y/4)*2,:);
Img3 = Img(:,round(y/4)*2:round(y/4)*3,:);
Img4 = Img(:,round(y/4)*3:end,:);

IMG1 = imresize(Img1,0.7);
IMG2 = imresize(Img2,0.7);
IMG3 = imresize(Img3,0.7);
IMG4 = imresize(Img4,0.7);

for i = 1:4
save(['/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_',num2str(i),'.mat'],['Img',num2str(i)], '-v7.3');
eval(['imwrite(IMG',num2str(i),',''/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Liebert_Institute_OTS-20-7690_rush_anterior_',num2str(i),'.tif'')']);
end

