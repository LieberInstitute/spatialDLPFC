function splitSlide(fname)
%fname = '/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/Lieber_Institute_OTS-20-7748_rush_posterior.tif';
numimgs = size(imfinfo(fname),1); %numimgs is the number of images in the multiplane tif file
N = 4; %number of capture areas
disp(['The multiplane tif has ',num2str(numimgs),' images'])

parfor i = 1:numimgs
    tic
    I{i}.image = imread(fname,i); 
    disp(['Imported image ',num2str(i), ' of the multiplane tif'])
    toc
end

%save tif image in mat format
tic
disp('Saving the multiplane tif to mat file')
save([fname(1:end-4),'.mat'],'I', '-v7.3');
toc

Img = I{1}.image; %whole slide image

clear I

[~,x,~] = size(Img);
[path1,name1,ext1] = fileparts(fname);

tic
disp('Splitting whole slide into individual capture areas')

Img1 = Img(:,1:round(x/N),:);
%imshow(Img1)
save([fullfile(path1,name1),'_A1.mat'],'Img1','-v7.3');
IMG1 = imresize(Img1,0.7);
imwrite(IMG1,[fullfile(path1,name1),'_A1',ext1])
clear Img1 IMG1

Img2 = Img(:,round(x/N):round(x/N)*2,:);
%imshow(Img2)
save([fullfile(path1,name1),'_B1.mat'],'Img2','-v7.3');
IMG2 = imresize(Img2,0.7);
imwrite(IMG2,[fullfile(path1,name1),'_B1',ext1])
clear Img2 IMG2

Img3 = Img(:,round(x/N)*2:round(x/N)*3,:);
%imshow(Img3)
save([fullfile(path1,name1),'_C1.mat'],'Img3','-v7.3');
IMG3 = imresize(Img3,0.7);
imwrite(IMG3,[fullfile(path1,name1),'_C1',ext1])
clear Img3 IMG3

Img4 = Img(:,round(x/N)*3:end,:);
%imshow(Img4)
save([fullfile(path1,name1),'_D1.mat'],'Img4','-v7.3');
IMG4 = imresize(Img4,0.7);
imwrite(IMG4,[fullfile(path1,name1),'_D1',ext1])
clear Img4 IMG4

toc 
