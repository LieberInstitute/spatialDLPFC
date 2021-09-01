function clickPlot(im, mask, R, tbl, count)

    imsize = size(im);
    spotRow = table2array(tbl(:, 5));
    spotCol = table2array(tbl(:, 6));
    nSpots = size(spotCol, 1);
    figure('units', 'normalized', 'position', [0 0 1 1]);
    ax = imagesc(im);
    set(ax, 'UserData', true);
    set(ax, 'ButtonDownFcn', {@click, im, mask})
    hold on;
    viscircles([spotCol spotRow], repelem(R, nSpots), 'Color', 'r');
    if ~isempty(count)
        text(spotCol, spotRow, string(count), 'Color', 'red', 'FontSize', 20);
    end
    uipanel('units', 'normalized', 'position', [0.095 0.04 0.15 0.05], 'Title', 'Barcode Lookup');
    lookupDialog = uicontrol('Style', 'edit', 'units', 'normalized', 'position', [0.1 0.05 0.075 0.02], 'String', '');
    lookupButton = uicontrol('Style', 'push', 'units', 'normalized', 'position', [0.185 0.05 0.05 0.02], 'String', 'Lookup');
    lookupButton.Callback = {@lookup, lookupDialog, tbl, R};
    zoomoutButton = uicontrol('Style', 'push', 'units', 'normalized', 'position', [0.25 0.05 0.05 0.02], 'String', 'Zoom Out');
    zoomoutButton.Callback = {@zoomout, imsize};
end

function click(hObj, ~, im, mask)

    if get(hObj, 'UserData')
        set(hObj, 'CData', mask);
        set(hObj, 'UserData', false);
    else
        set(hObj, 'CData', im);
        set(hObj, 'UserData', true);
    end

end

function lookup(~, ~, lookupDialog, tbl, R)

    barcodes = table2array(tbl(:,1));
    findThis = lookupDialog.String;
    idx = contains(barcodes, findThis);
    try
        row = table2array(tbl(idx, 5));
        row = row(1);
        col = table2array(tbl(idx, 6));
        col = col(1);
        xlim([col-R-5 col+R+5]);
        ylim([row-R-5 row+R+5]);
    catch
        errordlg("Could not find this barcode");
    end
    

end

function zoomout(~, ~, imsize)

    xlim([1 imsize(2)]);
    ylim([1 imsize(1)]);

end