function [count,prop] = countSpots(BW, R, tbl, posPath)
  
count = [];
prop = [];

O = fieldnames(BW);

    nSpots = size(tbl, 1);
    disp([num2str(nSpots),' Visium spots detected'])
    
    disp('Building spot grid')
    crow = round(table2array(tbl(:, 5)));
    ccol = round(table2array(tbl(:, 6)));
    mask = zeros(size(BW.(O{1})));
    mask(sub2ind(size(mask),crow,ccol)) = 1;
    mask = bwdist(mask) <= R;
    mask = bwlabel(mask);
    
for C = 1:numel(O)
    [BW.(O{C}),numROI] = bwlabel(BW.(O{C}));
	disp([num2str(numROI),' ', O{C},' ROIS detected'])
    count.(O{C}) = zeros(nSpots, 1);
    prop.(O{C}) = zeros(nSpots, 1);
end

tic
disp('counting dots per Visium spot')

for i = 1:nSpots 
    idx = mask(crow(i), ccol(i));
    spot = regionprops(mask==idx);
      for C = 1:numel(O)
        signal = struct2table(regionprops(mask==idx & BW.(O{C})>0));
        points = signal.Centroid;
        if isempty(points)
            isincircle = 0;
        else 
        isincircle = sum((points - [ccol(i) crow(i)]).^2,2)<= R^2;
        end
        %check
%         [tempx,tempy] = find(mask == idx);
%         temp = BW.(O{C})(min(tempx):max(tempx),min(tempy):max(tempy));
%         imshow(temp)
%         viscircles(size(temp)/2, repelem(R, 1), 'Color', 'r');

        count.(O{C})(i) = length(find(isincircle));
        prop.(O{C})(i) = sum([signal.Area])/spot.Area;
      end
      
      if mod(i,100) == 0
        disp([num2str(i),' spots finished in time ', num2str(toc),'s'])
      end

end

for C = 1:numel(O)
 temp = [count.(O{C}) prop.(O{C})];  
 tbl = [tbl array2table(temp, 'VariableNames', {['N',O{C}],['P',O{C}]})];        
end

if ~exist(posPath, 'dir')
   mkdir(posPath);
end
    
    disp('writing table')
    writetable(tbl, fullfile(posPath, 'tissue_spot_counts.csv'), 'Delimiter', ',');

end
