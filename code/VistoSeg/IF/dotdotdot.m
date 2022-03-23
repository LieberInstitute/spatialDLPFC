 function dotdotdot(fname)

%fname = '/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/Images/VisiumIF/VistoSeg/V10B01-087_A1.mat';
img = load(fname);
[Y, X] = size(img.DAPI);

disp(fieldnames(img))
OUTimg = [];

O = fieldnames(img);
for i = 1:numel(O)
	  channe_i = ['rescale(img.',O{i},')'];
       
  if contains(channe_i, 'AF')
      channel = eval(channe_i);
  	  Lip = i;
  elseif contains(channe_i,'DAPI')
	  channel = medfilt3(eval(channe_i),[3 3 3]);
	  DAPI = i;
  else
      channel = imhmin(eval(channe_i),2*std2(eval(channe_i)));% suppress background noise in RNA scope channels.
  end
  
   thresh = graythresh(channel); %for rosehip
    
   %  if thresh<0.1
   %      thresh = 0.1;
   %  end
   %  if contains(channe_i,'Cy5')
   %     thresh = 0.2;
   %  end
     
     
    BWc = imbinarize(channel,thresh);
 
%if contains(channe_i,'DAPI')
    
    BWc = imfill(BWc,'holes');
    x = imcomplement(channel);
    x = imhmin(x,std(channel(:)));
	L = watershed(x);
	BWc(L==0) = 0;
     
	bw3=max(BWc,[],3);
    D = -bwdist(~bw3);
    mask = imextendedmin(D,2);
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    bw3(Ld2 == 0) = 0;
	BWc(bw3==0) = 0;	
	
else	
	 x = imcomplement(channel);
	 x = imhmin(x,2*std(channel(:)));
	 L = watershed(x);
	 BWc(L==0) = 0;  
end	 

v = ['segmentation.',O{i}];
eval([v '= BWc;']);v = ['segmentation.',O{i}];
eval([v '= BWc;']);

%OUTimg = [OUTimg ,[max(rescale(img.(O{i})[10001:10500,1001:10500]),[],3),ones(Y,5); ones(5,X), ones(5,5); max(BWc,[],3), ones(Y,5)]];
end
