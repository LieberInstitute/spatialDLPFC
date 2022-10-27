function clickPlot(im, mask, R, tbl, count, O)

    imsize = size(im.(O{1}));
    spotRow = table2array(tbl(:, 5));
    spotCol = table2array(tbl(:, 6));
    nSpots = size(spotCol, 1);
    figure('Name','spotspotcheck','NumberTitle','off','units', 'normalized', 'position', [0 0 1 1]);
    ax = imagesc(im.(O{1}));
    set(ax, 'UserData', true);
    cmap = colormap(gca);
    hold on;
    
    daspect([1 1 1])
    viscircles([spotCol spotRow], repelem(R, nSpots), 'Color', 'r');
    
    uipanel('units', 'normalized', 'position', [0.095 0.04 0.15 0.05], 'Title', 'Barcode Lookup');
    lookupDialog = uicontrol('Style', 'edit', 'units', 'normalized', 'position', [0.1 0.05 0.075 0.02], 'String', '');
    lookupButton = uicontrol('Style', 'push', 'units', 'normalized', 'position', [0.185 0.05 0.05 0.02], 'String', 'Lookup');
    lookupButton.Callback = {@lookup, lookupDialog, tbl, R};
    zoomoutButton = uicontrol('Style', 'push', 'units', 'normalized', 'position', [0.3 0.05 0.1 0.02], 'String', 'Zoom Out');
    zoomoutButton.Callback = {@zoomout, imsize};
    channelDD = uicontrol('Style','popupmenu','units', 'normalized','position', [0.5 0.05 0.1 0.02], 'String', O);
    channelDD.Callback = {@selection,im,O,ax,cmap,count,spotCol,spotRow} ;
    set(ax, 'ButtonDownFcn', {@click, im, mask, O, channelDD, cmap})

     if ~isempty(count)
         C = get(channelDD,'Value');
         text(spotCol, spotRow, string(count.(O{C})), 'Color', 'red', 'FontSize', 20, 'Clipping','on');
     end
   
end

function click(hObj, ~, im, mask, O, channelDD, cmap)

    C = get(channelDD,'Value');
    if get(hObj, 'UserData')
        set(hObj, 'CData', mask.(O{C}));
        set(hObj, 'UserData', false);
        colormap(gca, gray(2))
        %if exist('count','var')
        %text(spotCol, spotRow, string(count.(O{C})), 'Color', 'red', 'FontSize', 20);
        %end
    else
        set(hObj, 'CData', im.(O{C}));
        set(hObj, 'UserData', true);
        colormap(gca, cmap)
        %if exist('count','var')
        %text(spotCol, spotRow, string(round(prop.(O{C}),1)), 'Color', 'red', 'FontSize', 20);
        %end
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

function selection(hObj,~,im,O,ax,cmap,count,spotCol,spotRow)
C = get(hObj,'Value');
set(ax, 'CData', im.(O{C}));
colormap(gca, cmap)
set(ax, 'UserData', true);
if ~isempty(count)
    text(spotCol, spotRow, string(count.(O{C})), 'Color', 'red', 'FontSize', 20);
end

end