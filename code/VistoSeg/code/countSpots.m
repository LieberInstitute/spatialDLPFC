function count = countSpots(BW, R, tbl, posPath)
    
count = [];


    nSpots = size(tbl, 1);
	[~, numNuclei] = bwlabel(BW);
	disp([num2str(numNuclei), ' cells detected'])
	disp([num2str(nSpots),' spots detected'])
		disp('counting nuclei per spot')
    crow = round(table2array(tbl(:, 5)));
    ccol = round(table2array(tbl(:, 6)));
    mask = zeros(size(BW));
    count = zeros(nSpots, 1);
    for i = 1:nSpots
        mask(crow(i), ccol(i)) = 1;
    end
    mask = bwdist(mask) <= R;
    mask = bwlabel(mask);
    BW = bwlabel(BW);
    for i = 1:nSpots
        idx = mask(crow(i), ccol(i));
        %tmpBW = BW;
        %tmpBW(mask~=idx) = 0;
        %[~, c] = bwlabel(tmpBW);
        spot = BW(mask==idx & BW>0);
        c = length(unique(spot));
        count(i) = c;
        if mod(i,100) == 0
        disp([num2str(i),' spots finished in time ', num2str(toc),'s'])
        end
        
    end
    tbl = [tbl array2table(count)];
    tbl.Properties.VariableNames = {'barcode','tissue','row','col','imagerow','imagecol','count'};
    if ~exist(posPath, 'dir')
        mkdir(posPath);
    end
    
    disp('writing table')
    writetable(tbl, fullfile(posPath, 'tissue_spot_counts.csv'), 'Delimiter', ',');

end
