function randomFields(fname)
he = imread(fname);

%load([fname(1:end-4),'_nuclei_WS_final.mat'])
load([fname(1:end-4),'_nuclei.mat'])
%load([fname(1:end-4),'_mask.mat'])
%mask_dark_blue = mask{M};

[x,y,~] = size(he);
xn = randi([round(x/3) round(x/3)*3],3,1);
yn = randi([round(y/3) round(y/3)*3],3,1);

A = he(xn(1):xn(1)+499, yn(1):yn(1)+499, :);
B = he(xn(2):xn(2)+499, yn(2):yn(2)+499, :);
C = he(xn(3):xn(3)+499, yn(3):yn(3)+499, :);
temp = [A;zeros(10,500,3);B;zeros(10,500,3);C];
imshow(temp)

A1 = mask_dark_blue(xn(1):xn(1)+499, yn(1):yn(1)+499);
B1 = mask_dark_blue(xn(2):xn(2)+499, yn(2):yn(2)+499);
C1 = mask_dark_blue(xn(3):xn(3)+499, yn(3):yn(3)+499);
temp1 = [A1;ones(10,500);B1;ones(10,500);C1];
figure, imshow(temp1)

close all
imwrite(temp,[fname(1:end-4),'_nuclei_raw.png'])
imwrite(temp1,[fname(1:end-4),'_nuclei_seg.png'])
end

