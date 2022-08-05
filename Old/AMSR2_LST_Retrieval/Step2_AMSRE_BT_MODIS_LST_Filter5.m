%% AMSRE BT and MODIS LST Filter.
%{
2019/01/10
  The modisPixelPercentage threshold was set 0.6 in the current experiment, which lead to the very
few matched AMSRE and upscaled MODIS pixels in some places like the northeastern China, and the
regression model, therefore, can not be constructed in these places. The strategy borrowing the
regression models in similar land cover zones is verified to be improper becasue they result in
huge errors.
  A new strategy is planned to set the threshold based on the heterogeneity of the landscape. In the
gentle slope or flat area with homogeneous land cover type, such as the northeastern China, the
threshold could be set as small as possible, while in the steep terrain or highly heterogeneous
area, the threshold could be set larger.
  The AMSRE BT images has no corresponding QC layer to control their qualites, so there are numbers
of abnormal pixels in the retrieved AMSRE LST images due to the bad AMSRE BT pixels. Try to develope
a stratedy to handle this problem. Using the STD of AMSRE LST block pixels or the time series of one
pixel.
%}

%% Flags and Preseted parameters.
% Flag assign the orbit. 1 means ascending (A, daytime), 2 means descending (D, night-time).
flg1 = 2;
% Flag control whether export QC controlled AMSRE BT and MODIS LST images. 1 means Yes, 0 means No.
flg2 = 1;
% Flag control whether export scatter plot of MODIS LST and AMSRE BT image in each channel. 1 means
%   Yes, 0 means No.
flg3 = 0;
% Flag control whether save the key variables into matlab mat files.
flg4 = 0;

orbit = {'A', 'D'};
dayNight = {'Day', 'Night'};

% Year of the operational data. Avaiable years: 2003 to 2011.
operationalYear = 2011;
yearString = num2str(operationalYear);

% QC that indicate the quality of MODIS LST value, include:
%   [0, 2, 3, 8, 17, 25, 65, 73, 81, 89, 129, 137, 145, 153].
% The meanings of QC can be found in 'F:\PaperFusionLST\Doc\MYD11A1_QC.csv'
qcControl = [0, 17];
% !!! The QC codes in HDF were modified automatically and inexplicably when converted to TIF. !!!
% qcControl = [65535, 34];

% Threshold that picks out the avaliable AMSRE pixels, which is the percentage of QC controlled
%   MODIS pixels in the extent of one AMSRE pixel.
modisPixelPercentage = 0.6;

%% Paths.
addpath('F:\PaperFusionLST\Matlab\Functions\');
dataRootPath = 'F:\PaperFusionLST\Data\';
figureRootPath = 'F:\PaperFusionLST\Figure\';

% ------------------------------------------ Input paths -------------------------------------------
% MODIS LST.
modisLstFolder = 'MYD11A1\';
modisLstMosaicFolder = 'MYD11A1_2_LST_MosaicCN_TIF\';
modisLstYearFolder = ['MYD11A1_', yearString, 'XXXX\'];
modisLstMosaicPath = [dataRootPath, modisLstFolder, modisLstMosaicFolder, modisLstYearFolder];

% AMSRE BT.
amsreGridFolder = 'AMSRE_QuarterDegreeGrid\';
amsreBtClipFolder = 'AMSRE_3_BT_ClipCN_Tif\';
amsreBtYearFolder = ['AMSRE_', yearString, 'XXXX\'];
amsreBtClipYearPath = [dataRootPath, amsreGridFolder, amsreBtClipFolder, amsreBtYearFolder];

% SRTM slope, CN extent and MODIS pixel area images.
srtmSlpCnPath = [dataRootPath, 'SRTM\1_SRTM_CN_TIF\SRTM_0.01d_CN_Slp.tif'];
cnExtentPath = [dataRootPath, 'Zones\GeographicalZones_merge.tif'];
modisCnAreaPath = [dataRootPath, modisLstFolder, modisLstMosaicFolder, 'MYD11A1.PixelArea.tif'];

% ----------------------------------------- Output paths -------------------------------------------
% MODIS LST.
modisLstMaskFolder = 'MYD11A1_5_LST_MaskCN_TIF\';
modisLstMaskYearPath = [dataRootPath, modisLstFolder, modisLstMaskFolder, modisLstYearFolder];
if ~exist(modisLstMaskYearPath, 'dir')
    mkdir(modisLstMaskYearPath);
end

% AMSRE BT.
amsreBtMaskFolder = 'AMSRE_5_BT_MaskCN_Tif\';
amsreBtMaskYearPath = [dataRootPath, amsreGridFolder, amsreBtMaskFolder, amsreBtYearFolder];
if ~exist(amsreBtMaskYearPath, 'dir')
    mkdir(amsreBtMaskYearPath);
end

% Figure.
outFigureFolder = 'ScatterAmsreModis2\';
outFigureYearFolder = ['Scatter_', yearString, 'XXXX\'];

%% Data lists.
% Read AMSRE BT daily folder list.
amsreBtClipDayFolderList = dir([amsreBtClipYearPath, 'AMSRE_*']);
amsreBtClipDayFolderList = {amsreBtClipDayFolderList.name}';
amsreBtClipDayFolderListN = length(amsreBtClipDayFolderList);

% Read MODIS LST, QC, and Emissivity daily file list.
modisLstMosaicDayList = dir([modisLstMosaicPath, 'MYD*LST*', dayNight{flg1}, '.tif']);
modisLstMosaicDayList = {modisLstMosaicDayList.name}';
modisLstMosaicDayListN = length(modisLstMosaicDayList);

modisQcMosaicDayList = dir([modisLstMosaicPath, 'MYD*QC*', dayNight{flg1}, '.tif']);
modisQcMosaicDayList = {modisQcMosaicDayList.name}';
modisQcMosaicDayListN = length(modisQcMosaicDayList);

modisEmis31MosaicDayList = dir([modisLstMosaicPath, 'MYD*Emis31.tif']);
modisEmis31MosaicDayList = {modisEmis31MosaicDayList.name}';
modisEmis31MosaicDayListN = length(modisEmis31MosaicDayList);

modisEmis32MosaicDayList = dir([modisLstMosaicPath, 'MYD*Emis32.tif']);
modisEmis32MosaicDayList = {modisEmis32MosaicDayList.name}';
modisEmis32MosaicDayListN = length(modisEmis32MosaicDayList);

% Read the date string lists of MODIS LST, QC, Emissivity and AMSRE BT images, and get the
%   intersection of them.
amsreBtDateFolderList = cell(amsreBtClipDayFolderListN, 1);
for i = 1 : amsreBtClipDayFolderListN
    amsreBtDateFolderList{i} = amsreBtClipDayFolderList{i}(7:14);
end

modisLstDateList = cell(modisLstMosaicDayListN, 1);
for i = 1 : modisLstMosaicDayListN
    modisLstDateList{i} = modisLstMosaicDayList{i}(10:17);
end

modisQcDateList = cell(modisQcMosaicDayListN, 1);
for i = 1 : modisQcMosaicDayListN
    modisQcDateList{i} = modisQcMosaicDayList{i}(10:17);
end

modisEmis31DateList = cell(modisEmis31MosaicDayListN, 1);
for i = 1 : modisEmis31MosaicDayListN
    modisEmis31DateList{i} = modisEmis31MosaicDayList{i}(10:17);
end

modisEmis32DateList = cell(modisEmis32MosaicDayListN, 1);
for i = 1 : modisEmis32MosaicDayListN
    modisEmis32DateList{i} = modisEmis32MosaicDayList{i}(10:17);
end

yearDateList = intersect(intersect(intersect(intersect(modisLstDateList, modisQcDateList), ...
    amsreBtDateFolderList), modisEmis31DateList), modisEmis32DateList);
yearDayListN = length(yearDateList);

%% Read references and attributes.
% Get the spatial reference of AMSRE BT images.
[~, amsreBtRef] = geotiffread([amsreBtClipYearPath, amsreBtClipDayFolderList{1},...
    '\AMSRE_D25_', yearString, '0101', orbit{flg1}, '_v03_06H.tif']);

% Get the spatial reference parameters of MODIS images. The spatial reference of MODIS LST, QC,
%   emissivity and pixel area are the same.
[~, modisLstRef] = geotiffread([modisLstMosaicPath, 'MYD11A1.A', yearString, '0101_LST_', ...
    dayNight{flg1}, '.tif']);
modisLstCellsizeX = modisLstRef.CellExtentInLongitude;
modisLstCellsizeY = modisLstRef.CellExtentInLatitude;
modisLstXMin = modisLstRef.LongitudeLimits(1);
modisLstXMax = modisLstRef.LongitudeLimits(2);
modisLstYMin = modisLstRef.LatitudeLimits(1);
modisLstYMax = modisLstRef.LatitudeLimits(2);
modisLstPixelLeftSideXVector = modisLstXMin : modisLstCellsizeX : modisLstXMax - modisLstCellsizeX;
modisLstPixelTopSideYVector = modisLstYMax : - modisLstCellsizeY : modisLstYMin + modisLstCellsizeY;

% Get the spatial reference of SRTM DTM images.
[~, srtmRef] = geotiffread(srtmSlpCnPath);

% Get the spatial reference parameters of cn extent zones image.
[~, cnExtentRef] = geotiffread(cnExtentPath);
cnExtentRowN = cnExtentRef.RasterSize(1);
cnExtentColN = cnExtentRef.RasterSize(2);
cnExtentCellsizeX = cnExtentRef.CellExtentInLongitude;
cnExtentCellsizeY = cnExtentRef.CellExtentInLatitude;
cnExtentXMin = cnExtentRef.LongitudeLimits(1);
cnExtentXMax = cnExtentRef.LongitudeLimits(2);
cnExtentYMin = cnExtentRef.LatitudeLimits(1);
cnExtentYMax = cnExtentRef.LatitudeLimits(2);
cnExtentPixelLeftSideXVector = cnExtentXMin : cnExtentCellsizeX : cnExtentXMax - cnExtentCellsizeX;
cnExtentPixelTopSideYVector = cnExtentYMax : - cnExtentCellsizeY : cnExtentYMin + cnExtentCellsizeY;

% --------------------------------------------------------------------------------------------------
% Get the starting/ending row numbers and column numbers from the paramter images (MODIS, SRTM and
%   AMSRE) that are aligned with the boundary of CN extent image.
% !!! The spatial resolution of CN extent image is the same as AMSRE BT image and the pixel boundary
%   of CN extent image is also aligned with that of AMSRE BT image. So the CN extent image is used
%   here for clipping the MODIS, SRTM and AMSRE images, then the origianl three steps (Clip, Filter
%   and Mask) can be combined !!!
[modisStartRow, modisEndRow, modisStartCol, modisEndCol] = getBdyRowCol(modisLstRef, cnExtentRef);
[srtmStartRow, srtmEndRow, srtmStartCol, srtmEndCol] = getBdyRowCol(srtmRef, cnExtentRef);
[amsreStartRow, amsreEndRow, amsreStartCol, amsreEndCol] = getBdyRowCol(amsreBtRef, cnExtentRef);

% Get the number of rows and columns of the block. 
blockColN = round(cnExtentCellsizeX / modisLstCellsizeX);  % cnExtentCellsizeX = amsreBtCellsizeX
blockRowN = round(cnExtentCellsizeY / modisLstCellsizeY);  % cnExtentCellsizeY = amsreBtCellsizeY
blockN = blockColN * blockRowN;

% Get the left and right row numbers and the top and bottom column numbers from the items product
%   image that are aligned with the boundary of the first pixel at the up-left corner of CN extent
%   image.
itemsPixelLeftSideXVector = modisLstPixelLeftSideXVector(modisStartCol : modisEndCol);
itemsPixelTopSideYVector = modisLstPixelTopSideYVector(modisStartRow : modisEndRow+1);
itemsColN = length(itemsPixelLeftSideXVector);
itemsRowN = length(itemsPixelTopSideYVector);

for i = 1 : cnExtentColN
    leftSideXDistance = abs(itemsPixelLeftSideXVector - cnExtentPixelLeftSideXVector(i));
    minLeftSideDistance = min(leftSideXDistance);
    if minLeftSideDistance <= modisLstCellsizeX
        itemsStartLeftSideCol = find(leftSideXDistance == minLeftSideDistance);
        break;
    end
end
itemsStartRightSideCol = itemsStartLeftSideCol + blockColN - 1;

for i = 1 : cnExtentRowN
    topSideYDistance = abs(itemsPixelTopSideYVector - cnExtentPixelTopSideYVector(i));
    minTopSideYDistance = min(topSideYDistance);
    if minTopSideYDistance <= modisLstCellsizeY
        itemsStartTopSideRow = find(topSideYDistance == minTopSideYDistance);
        break;
    end
end
itemsStartBottomSideCol = itemsStartTopSideRow + blockRowN - 1;

% --------------------------------------------------------------------------------------------------
% Create the arrays restoring the statistics of AMSRE BT and MODIS LST images.
if flg3 == 1
    pixelPercentList = zeros(yearDayListN, 1) * nan;
    rArray = zeros(yearDayListN, 12) * nan;
    rmseArray = zeros(yearDayListN, 12) * nan;
    availablePixelCountArray = zeros(cnExtentRowN, cnExtentColN);
end

% Process the matched daily MODIS and AMSRE images.
[srtmSlpLayer, ~] = geotiffread(srtmSlpCnPath);
[modisAreaLayer, ~] = geotiffread(modisCnAreaPath);
[cnExtentLayer, ~] = geotiffread(cnExtentPath);
modisAreaLayer = double(modisAreaLayer(modisStartRow:modisEndRow+1,modisStartCol:modisEndCol))*0.01;
srtmSlpLayer = double(srtmSlpLayer(srtmStartRow : srtmEndRow+1, srtmStartCol : srtmEndCol));
cosGamaLayer = cosd(srtmSlpLayer);
for i = 1 : yearDayListN
    yearDay = yearDateList{i};
    disp(yearDay);
    
    % Get the matched MODIS LST, QC, emissivity filenames and AMSRE BT folder names.
    amsreBtDayFolderName = amsreBtClipDayFolderList{strcmp(yearDay, amsreBtDateFolderList)};
    modisLstDayName = modisLstMosaicDayList{strcmp(yearDay, modisLstDateList)};
    modisQcDayName = modisQcMosaicDayList{strcmp(yearDay, modisQcDateList)};
    modisEmis31DayName = modisEmis31MosaicDayList{strcmp(yearDay, modisEmis31DateList)};
    modisEmis32DayName = modisEmis32MosaicDayList{strcmp(yearDay, modisEmis32DateList)};
    
    % Check whether the AMSRE BT images and MODIS LST image exist.
    outAmsreBtDayFolderPath = [amsreBtMaskYearPath, amsreBtDayFolderName, '\'];
    outModisLstFilteredPath = [modisLstMaskYearPath, modisLstDayName];
    if exist(outAmsreBtDayFolderPath, 'dir')
        amsreBtN = length(dir([outAmsreBtDayFolderPath, 'AMSRE*', orbit{flg1}, '_v03*.tif']));
    else
        amsreBtN = 0;
    end
    if exist(outModisLstFilteredPath, 'file') && amsreBtN == 12
        continue;
    end
    
    % ----------------------------------------------------------------------------------------------
    % Read AMSRE BT images in 12 channels and control their qualities.
    amsreBtDayFolderPath = [amsreBtClipYearPath, amsreBtDayFolderName, '\'];
    amsreBtList = dir([amsreBtDayFolderPath, 'AMSRE*', orbit{flg1}, '_v03*.tif']);
    amsreBtList = {amsreBtList.name}';
    for ii = 1 : length(amsreBtList)
        if contains(amsreBtList{ii}, 'TIM')
            timIdx = ii;
            break;
        end
    end
    amsreBtList(ii) = [];
    amsreBtListN = length(amsreBtList);
    
    amsreBtArray = zeros(cnExtentRowN, cnExtentColN, amsreBtListN) * nan;
    for ii = 1 : amsreBtListN
        [amsreBtLayer, ~] = geotiffread([amsreBtDayFolderPath, amsreBtList{ii}]);
        amsreBtLayer = amsreBtLayer(amsreStartRow : amsreEndRow, amsreStartCol : amsreEndCol);
        amsreBtArray(:, :, ii) = double(amsreBtLayer) * 0.1;
    end
    amsreBtArray(amsreBtArray == 0) = nan;
    
    % Exclude the pixel with PR higher than 1.
    amsrePrArray = zeros(cnExtentRowN, cnExtentColN, amsreBtListN / 2) * nan;
    for ii = 1 : amsreBtListN / 2
        amsrePrArray(:, :, ii) = amsreBtArray(:, :, ii*2-1) ./ amsreBtArray(:, :, ii*2);
    end
    amsrePrIndexArray = amsrePrArray > 1;
    amsrePrIndexLayer = logical(sum(amsrePrIndexArray, 3));
    amsrePrIndexArray = repmat(amsrePrIndexLayer, 1, 1, amsreBtListN);
    amsreBtArray(amsrePrIndexArray) = nan;
    
    % Read MODIS LST, QC and emissivity images.
    [modisLstLayer, ~] = geotiffread([modisLstMosaicPath, modisLstDayName]);
    [modisQcLayer, ~] = geotiffread([modisLstMosaicPath, modisQcDayName]);
    [modisEmis31Layer, ~] = geotiffread([modisLstMosaicPath, modisEmis31DayName]);
    [modisEmis32Layer, ~] = geotiffread([modisLstMosaicPath, modisEmis32DayName]);
    
    % Clip MODIS LST, QC and emissivity images to the extent of AMSRE BT image.
    modisLstLayer = modisLstLayer(modisStartRow : modisEndRow+1, modisStartCol : modisEndCol);
    modisQcLayer = modisQcLayer(modisStartRow : modisEndRow+1, modisStartCol : modisEndCol);
    modisEmis31Layer = modisEmis31Layer(modisStartRow : modisEndRow+1, modisStartCol : modisEndCol);
    modisEmis32Layer = modisEmis32Layer(modisStartRow : modisEndRow+1, modisStartCol : modisEndCol);
    
    % Convert the value of MODIS LST and emissivity image.
    modisLstLayer = double(modisLstLayer) * 0.02;
    modisLstLayer(~ismember(modisQcLayer, qcControl * 0.02) | (modisLstLayer == 65535*0.02)) = nan;
    modisEmis31Layer = double(modisEmis31Layer) * 0.002 + 0.49;
    modisEmis31Layer(modisEmis31Layer == 0 * 0.002 + 0.49) = nan;
    modisEmis32Layer = double(modisEmis32Layer) * 0.002 + 0.49;
    modisEmis32Layer(modisEmis32Layer == 0 * 0.002 + 0.49) = nan;
    
    % ----------------------------------------------------------------------------------------------
    % Upscaling the MODIS LST images.
    modisBbeLayer = 0.273 + 1.778 * modisEmis31Layer - 1.807 * modisEmis31Layer .* ...
        modisEmis32Layer - 1.037 * modisEmis32Layer + 1.774 * modisEmis32Layer .^ 2;
    numeratorLayer = modisAreaLayer .* modisBbeLayer .* modisLstLayer .^ 4 .* ...
        secd(srtmSlpLayer) ./ cosGamaLayer;
    numeratorLayer(numeratorLayer < 0) = nan;
    denominatorLayer = modisBbeLayer .* modisAreaLayer .* secd(srtmSlpLayer);
    
    [numeratorSumLayer, denominatorSumLayer] = deal(zeros(cnExtentRowN, cnExtentColN) * nan);
    for ii = 1 : cnExtentRowN
        for jj = 1 : cnExtentColN
            % Locate each items product block.
            itemsBlockTopRow = itemsStartTopSideRow + blockRowN * (ii - 1);
            itemsBlockBottomRow = itemsStartBottomSideCol + blockRowN * (ii - 1);
            itemsBlockLeftCol = itemsStartLeftSideCol + blockColN * (jj - 1);
            itemsBlockRightCol = itemsStartRightSideCol + blockColN * (jj - 1);
            if itemsBlockRightCol > itemsColN || itemsBlockBottomRow > itemsRowN
                continue;
            end
            % Block calculation.
            numeratorBlock = numeratorLayer(itemsBlockTopRow : itemsBlockBottomRow, ...
                itemsBlockLeftCol: itemsBlockRightCol);
            denominatorBlock = denominatorLayer(itemsBlockTopRow : itemsBlockBottomRow, ...
                itemsBlockLeftCol: itemsBlockRightCol);
            blockAvailableIndex = ~isnan(numeratorBlock) & ~isnan(denominatorBlock);
            availablePixelRatio = sum(blockAvailableIndex(:)) / blockN;
            % Skip the block when the ratio of avaliable MODIS LST pixels is less than 60%.
            if availablePixelRatio < modisPixelPercentage
                continue;
            end
            numeratorBlock(~blockAvailableIndex) = nan;
            denominatorBlock(~blockAvailableIndex) = nan;
            numeratorSumLayer(ii, jj) = nansum(numeratorBlock(:));
            denominatorSumLayer(ii, jj) = nansum(denominatorBlock(:));
        end
    end
    N = ones(cnExtentRowN, cnExtentColN) * 4;
    upScaledModisLstLayer = nthroot(numeratorSumLayer ./ denominatorSumLayer, N);
    
    % Get the common available pixels in AMSRE BT and resampled MODIS LST images.
    amsreNanPixelIndex = isnan(amsreBtArray(:, :, 1));
    modisNanPixelIndex = isnan(upScaledModisLstLayer);
    outCnExtentPixelIndex = (cnExtentLayer == -128);
    uselessPixelIndex = amsreNanPixelIndex | modisNanPixelIndex | outCnExtentPixelIndex;
    upScaledModisLstLayer(uselessPixelIndex) = nan;
    for ii = 1 : amsreBtListN
        amsreBtLayer = amsreBtArray(:, :, ii);
        amsreBtLayer(uselessPixelIndex) = nan;
        amsreBtArray(:, :, ii) = amsreBtLayer;
    end
    
    % Export the QC controlled MODIS LST and AMSRE BT images.
    if ~exist(outAmsreBtDayFolderPath, 'dir')
        mkdir(outAmsreBtDayFolderPath);
    end
    for ii = 1 : amsreBtListN
        outAmsreBtFilteredPath = [outAmsreBtDayFolderPath, amsreBtList{ii}];
        geotiffwrite(outAmsreBtFilteredPath, single(amsreBtArray(:, :, ii)), cnExtentRef);
    end
    geotiffwrite(outModisLstFilteredPath, single(upScaledModisLstLayer), cnExtentRef);
    
    % ----------------------------------------------------------------------------------------------
    % Calculate the statistics and export the scatter plot of MODIS LST and AMSRE BT image in 12
    %   channels.
    if flg3 == 1
        % Calculate the CC of upscaled MODIS LST and AMSRE BT in 12 channels.
%         for ii = 1 : amsreBtListN
%             r = corrcoef(upScaledModisLstLayer, amsreBtArray(:, :, ii), 'rows', 'complete');
%             rArray(i, ii) = r(2);
%             rmseArray(i, ii) = rmse(upScaledModisLstLayer, amsreBtArray(:, :, ii));
%         end
        % Calculate the percentage of usable pixels in the extent of experimental area.
        pixelPercentList(i) = 1 - sum(uselessPixelIndex(:)) / numel(uselessPixelIndex);
        % Calculate the count of available days in each pixel of the extent of China in one year.
        availablePixelCountArray = availablePixelCountArray + ~uselessPixelIndex;

        % Export scatters.
        outFigurePath = [dataRootPath, outFigureFolder];
        if ~exist(outFigurePath, 'dir')
            mkdir(outFigurePath);
        end
        outFigureYearPath = [outFigurePath, outFigureYearFolder];
        if ~exist(outFigureYearPath, 'dir')
            mkdir(outFigureYearPath);
        end
        outFigureDayFolderPath = [outFigureYearPath, amsreBtDayFolderName, '\'];
        if ~exist(outFigureDayFolderPath, 'dir')
            mkdir(outFigureDayFolderPath);
        end
        for ii = 1 : amsreBtListN
            outScatterPath = [outFigureDayFolderPath, replace(amsreBtList{ii}, '.tif', '.png')];
            if exist(outScatterPath, 'file')
                continue;
            end
            
            f1 = figure; f1.Visible = 'off';
            plot(modisMeanLstInAmsrePixel, amsreBtArray(:, :, ii), '.k', ...
                [240, 320], [240, 320], 'r');
            xlabel('MODIS LST'); ylabel('AMSRE BT');
            title(replace(amsreBtList{ii}, '_', '\_'));
            
            txt1 = ['Pixel Percentage: ', num2str(pixelPercentList(i))];
            txt2 = ['R: ', num2str(rArray(i, ii))];
            txt3 = ['RMSE: ', num2str(rmseArray(i, ii))];
            text(0.05, 0.95, txt1, 'Units', 'normalized', 'FontSize', 12);
            text(0.05, 0.89, txt2, 'Units', 'normalized', 'FontSize', 12);
            text(0.05, 0.83, txt3, 'Units', 'normalized', 'FontSize', 12);
            print(f1, outScatterPath, '-dpng');
        end
        close all;
    end
    
end

% --------------------------------------------------------------------------------------------------
% Save the variables to matlab mat files.
if flg4 == 1
    geotiffwrite([dataRootPath, 'Data\availablePixelCount_', num2str(operationalYear), '.tif'], ...
        availablePixelCountArray, amsreBtRef);
    save([dataRootPath, 'pixelPercentList.mat'], 'pixelPercentList');
    save([dataRootPath, 'rArray.mat'], 'rArray');
    save([dataRootPath, 'rmseArray.mat'], 'rmseArray');
end
