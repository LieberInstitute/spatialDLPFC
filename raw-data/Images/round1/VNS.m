function VNS(fname, N)
%N = 5;
Img1 = imread(fname);
disp('adjusting contrast in image')
Img1_smooth = imgaussfilt(Img1,4);
Img1_smooth_adj = imadjust(Img1_smooth, [.2 .3 0; .6 .7 1],[]);

he = Img1_smooth_adj;
disp('converting adjusted r*g*b to L*a*b space')
lab_he = rgb2lab(he);
ab = lab_he(:,:,2:3);
ab = im2single(ab);
disp('running Kmeans')
tic
pixel_labels = imsegkmeans(ab,N,'NumAttempts',3);
toc

disp('printing clusters')
parfor i = 1:N
mask{i} = pixel_labels==i;
cluster{i} = he .* uint8(mask{i});
imwrite(cluster{i},[fname(1:end-4),'_cluster',num2str(i),'.tif']) 
end

disp('saving data')
save([fname(1:end-4),'_mask.mat'],'mask','-v7.3')
save([fname(1:end-4),'_cluster.mat'],'cluster','-v7.3')
save([fname(1:end-4),'.mat'],'Img1','-v7.3')

%I = 1;imshow(mask(I))%to check the matching mask

   
