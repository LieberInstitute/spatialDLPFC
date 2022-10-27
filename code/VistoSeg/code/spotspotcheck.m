function spotspotcheck()

    fig = uifigure('Name', 'spotspotcheck', 'position', [360 198 260 100]);
    launchButton = uibutton(fig, 'push', 'position', [30 20 200 20], 'Text', 'Start');
    countCheck = uicheckbox(fig, 'position', [30 70 200 20], 'Text', 'Get Cell Counts', 'Value', 0);
    launchButton.ButtonPushedFcn = {@start, countCheck};
    
end

function start(~, ~, countCheck)

    com.mathworks.mwswing.MJFileChooserPerPlatform.setUseSwingDialog(1)
    getImg = '*.tif';
    if getenv('HOMEPATH')
        homepath = getenv('HOMEPATH');
    elseif getenv('HOME')
        homepath = getenv('HOME');
    end
    if ~exist(fullfile(homepath, 'SpotSpot_paths'), 'dir')
        mkdir(fullfile(homepath, 'SpotSpot_paths'));
    end
    if exist(fullfile(homepath, 'SpotSpot_paths'), 'file')
        try
            pathinfo = fileread(fullfile(homepath, 'SpotSpot_paths', 'SpotSpot_paths.txt'));
            pathinfo = regexp(pathinfo, '\n', 'split');
            getImg = fullfile(pathinfo{1}, '*.tif');
        catch
        end
    end
    [imgFile, imgPath] = uigetfile(getImg, 'Select Visium Image');
    getMask = fullfile(imgPath, '*.mat');
    getJSON = fullfile(imgPath, '*.json');
    getPositions = fullfile(imgPath, '*.txt;*.csv'); 
    [maskFile, maskPath] = uigetfile(getMask, 'Select Segmented Mat File');
    [jsonFile, jsonPath] = uigetfile(getJSON, 'Select Scale Factors JSON File');
    [posFile, posPath] = uigetfile(getPositions, 'Select Tissue Positions/Spot Counts File');
    disp('loading.....')
    if ischar(imgPath)
        ip = strrep(imgPath, '\', '/');
        fid = fopen(fullfile(homepath, 'SpotSpot_paths', 'SpotSpot_paths.txt'), 'w');
        fprintf(fid, ip);
        fprintf(fid, '\n');
        fclose(fid);
        clear ip
    end

    BW = load(fullfile(maskPath, maskFile));
    O = fieldnames(BW);
    %O1 =[2,1,4,6,5,3]; %only for AD
    O1 = 1:numel(O);
    for C = 1:numel(O)
    im.(O{C}) = imread(fullfile(imgPath,imgFile),O1(C));
    end
    
    jsonname = fullfile(jsonPath, jsonFile);
    w = jsondecode(fileread(jsonname));
    R = ceil(w.spot_diameter_fullres/2);
    tbl = readtable(fullfile(posPath, posFile));
    count = [];
    %prop = [];
    
    if size(tbl, 2) > 6
        b = 7;
        for C = 1:numel(O)
        count.(O{C}) = table2array(tbl(:, b));
        %prop.(O{C}) = table2array(tbl(:, b+1));
        b = b+3;
        end
    end
    if countCheck.Value
        [count,~] = countSpots(BW, R, tbl, posPath);
    end
    clickPlot(im, BW, R, tbl, count, O);    

end
