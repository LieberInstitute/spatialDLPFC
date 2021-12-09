%filename = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/segmentations/V10B01-086_B1_nuclei_Run1.mat';
filename = '/Users/heenadivecha/Desktop/V10B01-086_B1_nuclei_Run1.mat';
load(filename)

BW = mask_dark_blue;

%check image
%imshow(BW(14000:14500, 8000:9000))
%temp = regionprops(BW(14000:14500, 8000:9000));
stats = regionprops(BW);

mask_dark_blue = bwareafilt(BW, [50 max([stats.Area])]);

save(mask_dark_blue, '/Users/heenadivecha/Desktop/V10B01-086_B1_nuclei_Run1.mat')

