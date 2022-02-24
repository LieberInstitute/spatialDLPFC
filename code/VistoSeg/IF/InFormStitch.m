function InFormStitch(filename,O,P,fname)

% O{1} = 'DAPI';
% O{2} = 'Abeta';
% O{3} = 'pTau';
% O{4} = 'GFAP';
% O{5} = 'MAP2';
% O{6} = 'Lipofuscin';

%filename = '/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/InForm/VIFAD1_V10A27-004/20210204_Visium-IF_Scan1*component_data.tif';

myfiles = dir(filename);

disp('Extracting coordinates')
tic
loc = cellfun(@(x) strsplit(x,'_'), {myfiles.name}',  'UniformOutput', false);
temp = cellfun(@(x) strsplit(x{P},','), loc,  'UniformOutput', false);
X = cellfun(@(x) str2double(x{1}(2:end)), temp);
%Y = cellfun(@(x) str2double(x{2}(1:end-1)), temp);
X = sort(unique(X));
%Y = sort(unique(Y));

toc

for C = 1:numel(O)
disp(['Stitching ', O{C}])
tic
Iy = [];
for x1 = 1:numel(X)
myfilesx = dir(fullfile(myfiles(1).folder,['*',num2str(X(x1)),',*component_data.tif']));
Ix = [];
for K = 1:numel(myfilesx)
fname1 = fullfile(myfilesx(K).folder,myfilesx(K).name);
temp = imread(fname1,C); 
Ix = [Ix;temp];
end
Iy = [Iy,Ix];
end
img.(O{C}) = Iy;
toc
%disp([O{C},' stitched'])
end

if ~exist('fname','var')
[~,fname,~] = fileparts(filename);
fname = strsplit(fname,'*');
fname = fname{1};
end

disp('Saving mat file')
save(fullfile(myfiles(1).folder,[fname,'.mat']) , '-struct','img','-v7.3')
end
