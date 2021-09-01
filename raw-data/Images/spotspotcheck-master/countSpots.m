function count = countSpots(BW, R, tbl, path)

    nSpots = size(tbl, 1);
    disp([num2str(nSpots),' spots detected'])
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
    if ~exist(path, 'dir')
        mkdir(path);
    end
    
    disp('writing table')
    writetable(tbl, fullfile(path, 'tissue_spot_counts.csv'), 'Delimiter', ',');

end
