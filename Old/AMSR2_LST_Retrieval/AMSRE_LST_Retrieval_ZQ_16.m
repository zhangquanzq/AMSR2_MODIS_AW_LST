%% AMSRE_LST_Retrival_ZhangQuan

%% Flags and Preseted parameters.
% Flag assign the daytime or night-time. 1 means Daytime, 2 means Night-time.
flg1 = 2;
% Flag assign the Zone type. 1 means Geographical, 2 means Land cover by Zhou.
flg2 = 1;
% Flag assign the linear or quadratic equation that the regression model apply. 1 means Linear, 2
%   means Quadratic.
flg3 = 2;
% Flag control whether export the zonal AMSRE BT and MODIS LST iamges and their scatterplots.
%   1 means Yes, 0 means No.
flg4 = 0;
% Flag assign the regression strategy. 1 means Year, 2 means Season, 3 means Month, 4 means
%   optimized model.
flg5 = 2;
% Flag assgin only get the regression coefficients or also retrieve the masked and full extent AMSRE
%   LST using single Year, Season, or Month models. 1 means only coefficients, 2 means also get LST.
% (If the flg5 was set as 4, the flg6 should only set as 1.)
flg6 = 1;
% Flag assgin the testing year. 1 means 2010, 2 means 2005. The data in 2010 is used for
%   constructing the regression model, while the date in 2005 is used for examing the effect of
%   regression model in other years.
flg7 = 1;

% Flag checking.
if flg5 == 4 && flg6 == 2
    message = ['The retrieving of AMSRE LST using single model is not supported when optimized', ...
        ' model was appointed!'];
   error(message); 
end

% The operational year.
testYears = {'2010', '2005', '2011'};
% Regression Model.
regressModels = {'linear_', 'quadratic_'};
% The orbit of AMSRE satellite, also indicates the day and night. A means day, D means night.
orbits = {'A', 'D'};
dayNight = {'Day', 'Night'};
% The regression strategy.
strategy = {'Year', 'Season', 'Month'};
% The percentage threshold of avaliable Desert, Water, Glacier, Buiding, or Snow pxin the
%   extent of one AMSRE pixel.
snowThreshold = 0.6;
lcThreshold = 0.6;
% The channel of AMSRE BT that is combined with MODIS LST to plot scatter.
channels = {'06h', '06v', '10h', '10v', '18h', '18v', '23h', '23v', '36h', '36v', '89h', '89v'};
channelsN = length(channels);
% The variable number involved in the regerssion model.
variablesN = channelsN + 3; % (12 channels) + (2 quadratic term) + (1 intercept).
% Season months.
seasonMonths = {[3, 4, 5], [6, 7, 8], [9, 10, 11], [12, 1, 2]};
% Regions dividing the landcover and zones.
%   [Step1_EasternDongbei: 1: 1 ~ 4, Step1_Huabei: 2: 5 ~ 8, Step1_Hua'nan: 3: 9 ~ 13]
%   [Step2_Xi'nan: 4: 14 ~ 17, Step2_EasternXibei: 5: 18 ~ 27, Step2_WesternDongbei: 6: 28 ~ 30]
%   [Step2_WesternXibei: 7: 31 ~ 53, Step3_QingzangPlateau: 8: 54 ~ 71].
regionNodes = [1, 5, 9, 14, 18, 28, 31, 54, 71+1];
regionsN = length(regionNodes) - 1;

%% Paths.
% Root Paths.
addpath('E:\PaperFusionLST\Code\Matlab\Functions\');
dataRootPath = 'E:\PaperFusionLST\Data\';
figureRootPath = 'E:\PaperFusionLST\Figure\';

% ------------------------------------------ Input paths -------------------------------------------
% AMSRE BT images.
amsreGridFolder = 'AMSRE_QuarterDegreeGrid\';
amsreBtCnFolder = 'AMSRE_3_BT_ClipCN_Tif\';
amsreBtMaskFolder = 'AMSRE_5_BT_MaskCN_Tif\';
amsreYearFolder = ['AMSRE_', testYears{flg7}, 'XXXX\'];
amsreBtCnYearPath = [dataRootPath, amsreGridFolder, amsreBtCnFolder, amsreYearFolder];
amsreBtMaskYearPath = [dataRootPath, amsreGridFolder, amsreBtMaskFolder, amsreYearFolder];

% MODIS LST images.
modisLstFolder = 'MYD11A1\';
modisLstCnFolder = 'MYD11A1_3_LST_ClipCN_TIF\';
modisLstMaskFolder = 'MYD11A1_5_LST_MaskCN_TIF\';
modisLstYearFolder = ['MYD11A1_', testYears{flg7}, 'XXXX\'];
modisLstCnYearPath = [dataRootPath, modisLstFolder, modisLstCnFolder, modisLstYearFolder];
modisLstMaskYearPath = [dataRootPath, modisLstFolder,modisLstMaskFolder,modisLstYearFolder];

% MODIS Snow images.
modisSnowFolder = 'MYD10C1\';
modisSnowCnFolder = 'MYD10C1_2_Snow_CN_TIF\';
modisSnowYearFolder = ['MYD10C1_', testYears{flg7}, 'XXXX\'];
modisSnowYearPath = [dataRootPath, modisSnowFolder, modisSnowCnFolder, modisSnowYearFolder];

% MODIS Land cover images.
modisLcFolder = 'Landcover\';
modisLcPath = [dataRootPath, modisLcFolder, 'MCD12Q1.A2010001.051_gcs_reclsfy_resmpl.tif'];
desertPath = [dataRootPath, modisLcFolder, 'Desert_CN_gcs.tif'];

% Zone image.
zonesFolder = 'Zones\';
zonesFileNameList = {'GeographicalZones_merge.tif', 'LandCoverZones.tif'};
zonesFilePath = [dataRootPath, zonesFolder, zonesFileNameList{flg2}];

% ----------------------------------------- Output paths -------------------------------------------
% Matlab variables.
matlabVariableFolder = 'Matlab\AMSRE_LST_Retrieval\';
matlabVariablePath = [dataRootPath, matlabVariableFolder];
if ~exist(matlabVariablePath, 'dir')
    mkdir(matlabVariablePath);
end

% Zonal AMSRE BT image.
outAmsreBtTifFolder = 'AMSRE_6_BT_Zonal_TIF\';
outAsmreBtTifPath = [dataRootPath, amsreGridFolder, outAmsreBtTifFolder];
outAmsreBtYearPath = [outAsmreBtTifPath, amsreYearFolder];
if ~exist(outAsmreBtTifPath, 'dir')
    mkdir(outAsmreBtTifPath);
end
if ~exist(outAmsreBtYearPath, 'dir')
    mkdir(outAmsreBtYearPath);
end

% Zonal MODIS LST image.
outModisLstTifFolder = 'MYD11A1_6_LST_Zonal_TIF\';
outModisLstTifPath = [dataRootPath, modisLstFolder, outModisLstTifFolder];
outModisLstYearPath = [outModisLstTifPath, modisLstYearFolder];
if ~exist(outModisLstTifPath, 'dir')
    mkdir(outModisLstTifPath);
end
if ~exist(outModisLstYearPath, 'dir')
    mkdir(outModisLstYearPath);
end

% Retreived AMSRE LST image.
outAmsreLstTifFolder = 'AMSRE_7_LST_TIF\';
outAsmreLstTifPath = [dataRootPath, amsreGridFolder, outAmsreLstTifFolder];
outAmsreLstYearPath = [outAsmreLstTifPath, amsreYearFolder];
if ~exist(outAsmreLstTifPath, 'dir')
    mkdir(outAsmreLstTifPath);
end
if ~exist(outAmsreLstYearPath, 'dir')
    mkdir(outAmsreLstYearPath);
end

% Figures.
outFigureFolder = 'ScatterAmsreModisZones\';
outFigureYearFolder = ['Scatter_', testYears{flg7}, 'XXXX\'];
outFigurePath = [figureRootPath, outFigureFolder];
outFigureYearPath = [outFigurePath, outFigureYearFolder];
if ~exist(outFigurePath, 'dir')
    mkdir(outFigurePath);
end
if ~exist(outFigureYearPath, 'dir')
    mkdir(outFigureYearPath);
end

%% Data Lists.
% Read the complete and masked AMSRE daily folder list.
amsreBtCnDayFolderList = dir([amsreBtCnYearPath, 'AMSRE*']);
amsreBtCnDayFolderList = {amsreBtCnDayFolderList.name}';
amsreBtCnDayFolderListN = length(amsreBtCnDayFolderList);

amsreBtMaskDayFolderList = dir([amsreBtMaskYearPath, 'AMSRE*']);
amsreBtMaskDayFolderList = {amsreBtMaskDayFolderList.name}';
amsreBtMaskDayFolderListN = length(amsreBtMaskDayFolderList);

% Read the complete and masked MODIS LST daily image list.
modisLstCnDayList = dir([modisLstCnYearPath, 'MYD11A1*LST*', dayNight{flg1}, '.tif']);
modisLstCnDayList = {modisLstCnDayList.name}';
modisLstCnDayListN = length(modisLstCnDayList);

modisLstMaskDayList = dir([modisLstMaskYearPath, 'MYD11A1*', dayNight{flg1}, '.tif']);
modisLstMaskDayList = {modisLstMaskDayList.name}';
modisLstMaskDayListN = length(modisLstMaskDayList);

% The number of days of the complete and masked data are the same in one year, but may not be the
%   same in different years. Pay attention to this difference if the regression model was created
%   using the data in one year but was applied on the data in another year.

% Check whether the number of masked AMSRE BT and MODIS LST images are the same.
if (amsreBtMaskDayFolderListN ~= modisLstMaskDayListN)
    error('The number of masked AMSRE BT and MODIS LST images are different!');
end

% Read the MODIS Snow image list.
modisSnowList = dir([modisSnowYearPath, 'MYD10C1*.tif']);
modisSnowList = {modisSnowList.name}';
modisSnowListN = length(modisSnowList);

% --------------------------------------------------------------------------------------------------
% Create the date list of the whole year.
firstDayOfYear = datenum([testYears{flg7}, '0101'], 'yyyymmdd');
endDayOfYear = datenum([testYears{flg7}, '1231'], 'yyyymmdd');
tempDateList = datevec(firstDayOfYear : endDayOfYear);
tempDateList = mat2cell(tempDateList(:, 1:3), ones(size(tempDateList, 1), 1), 3);
yearDateList = cellfun(@(t){sprintf('%4i%02i%02i', t)}, tempDateList);
yearDateListN = length(yearDateList);
yearMonthList = zeros(yearDateListN, 1);
for i = 1 : yearDateListN
    yearMonthList(i) = str2double(yearDateList{i}(5:6));
end

% Get the date list of each data category.
amsreBtCnDateList =  cell(amsreBtCnDayFolderListN, 1);
for i = 1 : amsreBtCnDayFolderListN
    amsreBtCnDateList{i} = amsreBtCnDayFolderList{i}(7:end);
end

modisLstCnDateList = cell(modisLstCnDayListN, 1);
modisLstCnMonthList = zeros(modisLstCnDayListN, 1);
for i = 1 : modisLstCnDayListN
    modisLstCnDateList{i} = modisLstCnDayList{i}(10:17);
    modisLstCnMonthList(i) = str2double(modisLstCnDayList{i}(14:15));
end

modisLstMaskDateList = cell(modisLstMaskDayListN, 1);
modisLstMaskMonthList = zeros(modisLstMaskDayListN, 1);
for i = 1 : modisLstMaskDayListN
    modisLstMaskDateList{i} = modisLstMaskDayList{i}(10:17);
    modisLstMaskMonthList(i) = str2double(modisLstMaskDayList{i}(14:15));
end

modisSnowDateList = cell(modisSnowListN, 1);
for i = 1 : modisSnowListN
    modisSnowDateList{i} = modisSnowList{i}(10:17);
end

%% Read references and attributes.
% Read the spatial reference of complete AMSRE BT images. The spatial references of all AMSRE BT
%   images are the same.
[~, btReference] = geotiffread([amsreBtCnYearPath, amsreBtCnDayFolderList{1}, ...
    '\AMSRE_D25_', testYears{flg7}, '0101', orbits{flg1}, '_v03_06H.tif']);
btRowN = btReference.RasterSize(1);
btColN = btReference.RasterSize(2);
btCellsizeX = btReference.CellExtentInLongitude;
btCellsizeY = btReference.CellExtentInLatitude;
btXMin = btReference.LongitudeLimits(1);
btXMax = btReference.LongitudeLimits(2);
btYMin = btReference.LatitudeLimits(1);
btYMax = btReference.LatitudeLimits(2);
btPixelLeftSideXVector = btXMin : btCellsizeX : btXMax - btCellsizeX;
btPixelTopSideYVector = btYMax : - btCellsizeY : btYMin + btCellsizeY;

% Read the spatial reference of masked AMSRE BT images and MODIS LST images. The spatial references
%   of all the masked AMSRE BT images and MODIS LST iamges are the same.
[~, btLstReference] = geotiffread([amsreBtMaskYearPath, amsreBtMaskDayFolderList{1}, ...
    '\AMSRE_D25_', testYears{flg7}, '0101', orbits{flg1}, '_v03_06H.tif']);
btLstRowN = btLstReference.RasterSize(1);
btLstColN = btLstReference.RasterSize(2);
btLstCellsizeX = btLstReference.CellExtentInLongitude;
btLstCellsizeY = btLstReference.CellExtentInLatitude;
btLstXMin = btLstReference.LongitudeLimits(1);
btLstXMax = btLstReference.LongitudeLimits(2);
btLstYMin = btLstReference.LatitudeLimits(1);
btLstYMax = btLstReference.LatitudeLimits(2);
btLstPixelLeftSideXVector = btLstXMin : btLstCellsizeX : btLstXMax - btLstCellsizeX;
btLstPixelTopSideYVector = btLstYMax : - btLstCellsizeY : btLstYMin + btLstCellsizeY;

% Read the spatial reference of MODIS land cover image. The spatial reference of desert image is the
% same as land cover image.
[~, lcReference] = geotiffread(modisLcPath);
lcRowN = lcReference.RasterSize(1);
lcColN = lcReference.RasterSize(2);
lcCellsizeX = lcReference.CellExtentInLongitude;
lcCellsizeY = lcReference.CellExtentInLatitude;
lcXMin = lcReference.LongitudeLimits(1);
lcXMax = lcReference.LongitudeLimits(2);
lcYMin = lcReference.LatitudeLimits(1);
lcYMax = lcReference.LatitudeLimits(2);
lcPixelLeftSideXVector = lcXMin : lcCellsizeX : lcXMax - lcCellsizeX;
lcPixelTopSideYVector = lcYMax : - lcCellsizeY : lcYMin + lcCellsizeY;

% Read the spatial reference of MODIS snow images. The spatial references of all the MODIS snow
%  images are the same.
[~, modisSnowReference] = geotiffread([modisSnowYearPath, modisSnowList{1}]);
snowRowN = modisSnowReference.RasterSize(1);
snowColN = modisSnowReference.RasterSize(2);
snowCellsizeX = modisSnowReference.CellExtentInLongitude;
snowCellsizeY = modisSnowReference.CellExtentInLatitude;
snowXMin = modisSnowReference.LongitudeLimits(1);
snowXMax = modisSnowReference.LongitudeLimits(2);
snowYMin = modisSnowReference.LatitudeLimits(1);
snowYMax = modisSnowReference.LatitudeLimits(2);
snowPixelLeftSideXVector = snowXMin : snowCellsizeX : snowXMax - snowCellsizeX;
snowPixelTopSideYVector = snowYMax : - snowCellsizeY : snowYMin + snowCellsizeY;

% --------------------------------------------------------------------------------------------------
% Get the starting/ending row numbers and column numbers from the complete AMSRE BT images that are
%   aligned with the boundary of AMSRE BT and MODIS LST images.
if btLstXMin > btXMin
    btStartCol = 1 + ceil((btLstXMin - btXMin) / btCellsizeX);
else
    btStartCol = 1;
end

if btYMax > btLstYMax
    btStartRow = 1 + ceil((btYMax - btLstYMax) / btCellsizeY);
else
    btStartRow = 1;
end

if btLstXMax < btXMax
    btEndCol = btColN - ceil((btXMax - btLstXMax) / btCellsizeX);
else
    btEndCol = btColN;
end

if btYMin < btLstYMin
    btEndRow = btRowN - ceil((btLstYMin - btYMin) / btCellsizeY);
else
    btEndRow = btRowN;
end

% --------------------------------------------------------------------------------------------------
% Get the starting/ending row numbers and column numbers from the fused land cover image that are
%   aligned with the boundary of AMSRE BT and MODIS LST images.
if btLstXMin > lcXMin
    lcStartCol = 1 + ceil((btLstXMin - lcXMin) / lcCellsizeX);
else
    lcStartCol = 1;
end

if lcYMax > btLstYMax
    lcStartRow = 1 + ceil((lcYMax - btLstYMax) / lcCellsizeY);
else
    lcStartRow = 1;
end

if btLstXMax < lcXMax
    lcEndCol = lcColN - ceil((lcXMax - btLstXMax) / lcCellsizeX);
else
    lcEndCol = lcColN;
end

if lcYMin < btLstYMin
    lcEndRow = lcRowN - ceil((btLstYMin - lcYMin) / lcCellsizeY);
else
    lcEndRow = lcRowN;
end

% Get the number of rows and columns of the block in land cover image.
lcBlockColN = round(btLstCellsizeX / lcCellsizeX);
lcBlockRowN = round(btLstCellsizeY / lcCellsizeY);
lcBlockN = lcBlockColN * lcBlockRowN;

% Get the left and right row numbers and the top and bottom column numbers from the land cover image
%   that are aligned with the boundary of the first pixel at the up-left corner of AMSRE BT and
%   MODIS LST images.
for ii = 1 : btLstColN
    leftSideXDistance = abs(lcPixelLeftSideXVector - btLstPixelLeftSideXVector(ii));
    minLeftSideDistance = min(leftSideXDistance);
    if minLeftSideDistance <= lcCellsizeX
        lcStartLeftSideCol = find(leftSideXDistance == minLeftSideDistance);
        break;
    end
end
lcStartRightSideCol = lcStartLeftSideCol + lcBlockColN - 1;

for ii = 1 : btLstRowN
    topSideYDistance = abs(lcPixelTopSideYVector - btLstPixelTopSideYVector(ii));
    minTopSideYDistance = min(topSideYDistance);
    if minTopSideYDistance <= lcCellsizeY
        lcStartTopSideRow = find(topSideYDistance == minTopSideYDistance);
        break;
    end
end
lcStartBottomSideCol = lcStartTopSideRow + lcBlockRowN - 1;

% --------------------------------------------------------------------------------------------------
% Get the starting/ending row numbers and column numbers from the MODIS snow image that are aligned
%   with the boundary of AMSRE BT and MODIS LST images.
if btLstXMin > snowXMin
    snowStartCol = 1 + ceil((btLstXMin - snowXMin) / snowCellsizeX);
else
    snowStartCol = 1;
end

if snowYMax > btLstYMax
    snowStartRow = 1 + ceil((snowYMax - btLstYMax) / snowCellsizeY);
else
    snowStartRow = 1;
end

if btLstXMax < snowXMax
    snowEndCol = snowColN - ceil((snowXMax - btLstXMax) / snowCellsizeX);
else
    snowEndCol = snowColN;
end

if snowYMin < btLstYMin
    snowEndRow = snowRowN - ceil((btLstYMin - snowYMin) / snowCellsizeY);
else
    snowEndRow = snowRowN;
end

% Get the number of rows and columns of the block in MODIS snow image.
snowBlockColN = round(btLstCellsizeX / snowCellsizeX);
snowBlockRowN = round(btLstCellsizeY / snowCellsizeY);
snowBlockN = snowBlockColN * snowBlockRowN;

% Get the left and right row numbers and the top and bottom column numbers from the MODIS snow image
%   that are aligned with the boundary of the first pixel at the up-left corner of AMSRE BT and
%   MODIS LST images.
for ii = 1 : btLstColN
    leftSideXDistance = abs(snowPixelLeftSideXVector - btLstPixelLeftSideXVector(ii));
    minLeftSideDistance = min(leftSideXDistance);
    if minLeftSideDistance <= snowCellsizeX
        snowStartLeftSideCol = find(leftSideXDistance == minLeftSideDistance);
        break;
    end
end
snowStartRightSideCol = snowStartLeftSideCol + snowBlockColN - 1;

for ii = 1 : btLstRowN
    topSideYDistance = abs(snowPixelTopSideYVector - btLstPixelTopSideYVector(ii));
    minTopSideYDistance = min(topSideYDistance);
    if minTopSideYDistance <= snowCellsizeY
        snowStartTopSideRow = find(topSideYDistance == minTopSideYDistance);
        break;
    end
end
snowStartBottomSideCol = snowStartTopSideRow + snowBlockRowN - 1;

%% Read Data.
% Create 3D arrays to restore the AMSRE BT, MODIS LST, and MODIS Snow images in the year, skip the
%   corresponding layer if there is no image in certain day.

% Read the complete AMSRE BT images from the 12 channels in the operational year.
outAmsreBtCnYearArrayMatPath = [matlabVariablePath, '1_AMSRE_BT_CN_', dayNight{flg1}, '_', ...
    testYears{flg7}, '.mat'];
if ~exist(outAmsreBtCnYearArrayMatPath, 'file')
    [amsreBt06hCnYearArray, amsreBt06vCnYearArray, amsreBt10hCnYearArray, ...
        amsreBt10vCnYearArray, amsreBt18hCnYearArray, amsreBt18vCnYearArray, ...
        amsreBt23hCnYearArray, amsreBt23vCnYearArray, amsreBt36hCnYearArray, ...
        amsreBt36vCnYearArray, amsreBt89hCnYearArray, amsreBt89vCnYearArray] ...
        = deal(zeros(btLstRowN, btLstColN, yearDateListN) * nan);
    parfor i = 1 : yearDateListN
        % Techniques handeling the issues in parallel computing.
        orbit = orbits;
        tempAmsreBtCnDayFolderList = amsreBtCnDayFolderList;
        
        % Read the AMSRE BT images in the day that the AMSRE BT images are avaliable.
        [~, amsreDayIndex] = ismember(yearDateList{i}, amsreBtCnDateList);
        if amsreDayIndex == 0, continue; end
        amsreFolderName = tempAmsreBtCnDayFolderList{amsreDayIndex};
        disp(['Read complete AMSRE BT image: ', amsreFolderName]);
        amsreBtCnDayPath = [amsreBtCnYearPath, amsreFolderName, '\'];
        amsreBtCnChannelList = dir([amsreBtCnDayPath, 'AMSRE_D*', orbit{flg1},'*.tif']);
        amsreBtCnChannelList = {amsreBtCnChannelList.name}';
        amsreBtCnChannelListN = length(amsreBtCnChannelList) - 1; % Remove AMSRE TIM image.
        amsreBtCnArray = zeros(btLstRowN, btLstColN, amsreBtCnChannelListN);
        for j = 1 : amsreBtCnChannelListN
            amsreBtCnLayer = geotiffread([amsreBtCnDayPath, amsreBtCnChannelList{j}]);
            amsreBtCnLayer = amsreBtCnLayer(btStartRow : btEndRow, btStartCol : btEndCol);
            amsreBtCnArray(:, :, j) = double(amsreBtCnLayer) * 0.1;
        end
        amsreBtCnArray(amsreBtCnArray <= 0) = nan;
        
        amsreBt06hCnYearArray(:, :, i) = amsreBtCnArray(:, :, 1);
        amsreBt06vCnYearArray(:, :, i) = amsreBtCnArray(:, :, 2);
        amsreBt10hCnYearArray(:, :, i) = amsreBtCnArray(:, :, 3);
        amsreBt10vCnYearArray(:, :, i) = amsreBtCnArray(:, :, 4);
        amsreBt18hCnYearArray(:, :, i) = amsreBtCnArray(:, :, 5);
        amsreBt18vCnYearArray(:, :, i) = amsreBtCnArray(:, :, 6);
        amsreBt23hCnYearArray(:, :, i) = amsreBtCnArray(:, :, 7);
        amsreBt23vCnYearArray(:, :, i) = amsreBtCnArray(:, :, 8);
        amsreBt36hCnYearArray(:, :, i) = amsreBtCnArray(:, :, 9);
        amsreBt36vCnYearArray(:, :, i) = amsreBtCnArray(:, :, 10);
        amsreBt89hCnYearArray(:, :, i) = amsreBtCnArray(:, :, 11);
        amsreBt89vCnYearArray(:, :, i) = amsreBtCnArray(:, :, 12);
    end
    save(outAmsreBtCnYearArrayMatPath, 'amsreBt06hCnYearArray', 'amsreBt06vCnYearArray', ...
        'amsreBt10hCnYearArray', 'amsreBt10vCnYearArray', 'amsreBt18hCnYearArray', ...
        'amsreBt18vCnYearArray', 'amsreBt23hCnYearArray', 'amsreBt23vCnYearArray', ...
        'amsreBt36hCnYearArray', 'amsreBt36vCnYearArray', 'amsreBt89hCnYearArray', ...
        'amsreBt89vCnYearArray');
else
    lia = ismember({'amsreBt06hCnYearArray', 'amsreBt06vCnYearArray', 'amsreBt10hCnYearArray',...
        'amsreBt10vCnYearArray', 'amsreBt18hCnYearArray', 'amsreBt18vCnYearArray', ...
        'amsreBt23hCnYearArray', 'amsreBt23vCnYearArray', 'amsreBt36hCnYearArray', ...
        'amsreBt36vCnYearArray', 'amsreBt89hCnYearArray', 'amsreBt89vCnYearArray'}, who);
    if sum(lia) < 12
        load(outAmsreBtCnYearArrayMatPath);
    end
end

% Read the masked AMSRE BT and MODIS LST images in the operational year.
outBtLstMaskYearArrayMatPath = [matlabVariablePath, '2_AMSRE_BT_MODIS_LST_Mask_', ...
    dayNight{flg1}, '_', testYears{flg7}, '.mat'];
if ~exist(outBtLstMaskYearArrayMatPath, 'file')
    [amsreBt06hMaskYearArray, amsreBt06vMaskYearArray, amsreBt10hMaskYearArray, ...
        amsreBt10vMaskYearArray, amsreBt18hMaskYearArray, amsreBt18vMaskYearArray, ...
        amsreBt23hMaskYearArray, amsreBt23vMaskYearArray, amsreBt36hMaskYearArray, ...
        amsreBt36vMaskYearArray, amsreBt89hMaskYearArray, amsreBt89vMaskYearArray, ...
        modisLstMaskYearArray] = deal(zeros(btLstRowN, btLstColN, yearDateListN) * nan);
    % The number of days of the maksed AMSRE BT and MODIS LST images are the same.
    parfor i = 1 : yearDateListN
        % Techniques handeling the issues in parallel computing.
        orbit = orbits;
        tempModisLstMaskDayList = modisLstMaskDayList;
        tempAmsreBtMaskDayFolderList = amsreBtMaskDayFolderList;
        
        % Get index of the day in which the AMSRE BT and MODIS LST images are avaliable.
        [~, btLstDayIndex] = ismember(yearDateList{i}, modisLstMaskDateList);
        if btLstDayIndex == 0, continue; end
        disp(['Read masked AMSRE BT and MODIS LST images: ', modisLstMaskDateList{btLstDayIndex}]);
        
        % Read the masked MODIS LST iamge.
        modisLstMaskName = tempModisLstMaskDayList{btLstDayIndex};
        modisLstMaskLayer = geotiffread([modisLstMaskYearPath, modisLstMaskName]);
        modisLstMaskLayer(modisLstMaskLayer < 0) = nan;
        
        % Read the masked AMSRE BT images in 12 channels.
        amsreBtMaskDayFolderName = tempAmsreBtMaskDayFolderList{btLstDayIndex};
        amsreBtMaskDayPath = [amsreBtMaskYearPath, amsreBtMaskDayFolderName, '\'];
        amsreBtMaskChannelList = dir([amsreBtMaskDayPath, 'AMSRE_D*', orbit{flg1},'*.tif']);
        amsreBtMaskChannelList = {amsreBtMaskChannelList.name}';
        amsreBtMaskChannelListN = length(amsreBtMaskChannelList);
        amsreBtMaskArray = zeros(btLstRowN, btLstColN, amsreBtMaskChannelListN);
        for j = 1 : amsreBtMaskChannelListN
            amsreBtMaskArray(:, :, j) = ...
                geotiffread([amsreBtMaskDayPath, amsreBtMaskChannelList{j}]);
        end
        amsreBtMaskArray(amsreBtMaskArray < 0) = nan;
        
        % Get the index layer and array used for masking the overlapped region of AMSRE BT and MODIS
        %   LST images. The masked AMSRE BT and MODIS LST images should already have the same
        %     coverage, the overlap index layer and array here is used to ensure the same coverage
        %     again in case of any anomaly.
        overlapIndexLayer = ~isnan(modisLstMaskLayer) & ~isnan(sum(amsreBtMaskArray, 3));
        overlapIndexArray = repmat(overlapIndexLayer, 1, 1, amsreBtMaskChannelListN);
        
        % Masking the overlap of MODIS LST and AMSRE BT images.
        modisLstMaskLayer(~overlapIndexLayer) = nan;
        modisLstMaskYearArray(:, :, i) = modisLstMaskLayer;
        
        amsreBtMaskArray(~overlapIndexArray) = nan;
        amsreBt06hMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 1);
        amsreBt06vMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 2);
        amsreBt10hMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 3);
        amsreBt10vMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 4);
        amsreBt18hMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 5);
        amsreBt18vMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 6);
        amsreBt23hMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 7);
        amsreBt23vMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 8);
        amsreBt36hMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 9);
        amsreBt36vMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 10);
        amsreBt89hMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 11);
        amsreBt89vMaskYearArray(:, :, i) = amsreBtMaskArray(:, :, 12);
    end
    save(outBtLstMaskYearArrayMatPath, 'modisLstMaskYearArray', ...
        'amsreBt06hMaskYearArray', 'amsreBt06vMaskYearArray', 'amsreBt10hMaskYearArray', ...
        'amsreBt10vMaskYearArray', 'amsreBt18hMaskYearArray', 'amsreBt18vMaskYearArray', ...
        'amsreBt23hMaskYearArray', 'amsreBt23vMaskYearArray', 'amsreBt36hMaskYearArray', ...
        'amsreBt36vMaskYearArray', 'amsreBt89hMaskYearArray', 'amsreBt89vMaskYearArray');
else
    lia = ismember({'modisLstMaskYearArray', ...
        'amsreBt06hMaskYearArray', 'amsreBt06vMaskYearArray', 'amsreBt10hMaskYearArray', ...
        'amsreBt10vMaskYearArray', 'amsreBt18hMaskYearArray', 'amsreBt18vMaskYearArray', ...
        'amsreBt23hMaskYearArray', 'amsreBt23vMaskYearArray', 'amsreBt36hMaskYearArray', ...
        'amsreBt36vMaskYearArray', 'amsreBt89hMaskYearArray', 'amsreBt89vMaskYearArray'}, who);
    if sum(lia) < 13
        load(outBtLstMaskYearArrayMatPath);
    end
end

% --------------------------------------------------------------------------------------------------
% Blend land cover image, desert image, zonal image and snow imges.
% Original land cover codes:
%   [BaredLand: 1000, VegetationLand: 2000, Water: 3000, Glacier: 4000, Building: 5000].
% Blended and recoded land cover codes:
%   [Others: 0, Water: 3000, Glacier: 4000, Building: 5000, Snow: 6000, Mixture: 10000].
% Desert code: [1000, 1100, ..., 2000];
% Zones Code: [1, 2, ..., 71];

% Get desert image and codes.
[desertLayer, ~] = geotiffread(desertPath);
desertLayer = single(desertLayer);
desertCodeList = unique(desertLayer); % [1000, 1100, ..., 2000, 65535]
desertNodataCode = desertCodeList(end);
desertCodeList = desertCodeList(1:end-1);
desertCodeListN = length(desertCodeList);

% Get Land cover iamge and codes.
[lcLayer, ~] = geotiffread(modisLcPath);
lcLayer = single(lcLayer) * 1000;
lcCodeList = unique(lcLayer);  % [1000, 2000, 3000, 4000, 5000]
uselessLcCodes = lcCodeList(1:2);
remainedLcCodes = lcCodeList(3:5);  % Water, Glacier, Building.
lcCodeList = [0; remainedLcCodes; 6000; 10000];
lcCodeListN = length(lcCodeList);

% Remove the original BaredLand and VegetationLand covers in Lc image, and then blend the rest with
%   desert. Keep the Water, Glacier, Building unchanged when they overlap with the desert.
lcLayer(ismember(lcLayer, uselessLcCodes)) = 0;
desertIndex = (desertLayer ~= desertNodataCode) & ~ismember(lcLayer, remainedLcCodes);
lcLayer(desertIndex) = desertLayer(desertIndex);

% Combine the Lc image and Zonal image. The land covers except snow are added into the zonal image
%   here. The area of snow cover varies day by day, so the snow cover is handled in the next step.
%   The mixed Lc pixel was defined as the pixel with no Lc type holding the area over 60% of the
%   pixel. The combined land cover code is determined by: zoneCode + lcCode. For example: 3010
%   indicates the Water in zone 10.
[otherPercentLayer, desertPercentLayer, waterPercentLayer, glacierPercentLayer, ...
    buildingPercentLayer, desertCodeLayer] = deal(zeros(btLstRowN, btLstColN) * nan);
[zonesLayer, ~] = geotiffread(zonesFilePath);
zonesLayer = single(zonesLayer);
zonesLayer(ismember(zonesLayer, [-128, 128])) = nan;
zonesLcLayer = zonesLayer;
for ii = 1 : btLstRowN
    for jj = 1 : btLstColN
        % Skip the Nodata pixels in zonal iamge.
        if isnan(zonesLcLayer(ii, jj))
            continue;
        end
        % Locate each Lc block.
        lcBlockTopRow = lcStartTopSideRow + lcBlockRowN * (ii - 1);
        lcBlockBottomRow = lcStartBottomSideCol + lcBlockRowN * (ii - 1);
        lcBlockLeftCol = lcStartLeftSideCol + lcBlockColN * (jj - 1);
        lcBlockRightCol = lcStartRightSideCol + lcBlockColN * (jj - 1);
        if lcBlockRightCol > lcColN || lcBlockBottomRow > lcRowN
            continue;
        end
        % Lc block calculation.
        lcBlock = lcLayer(lcBlockTopRow : lcBlockBottomRow, lcBlockLeftCol : lcBlockRightCol);
        otherPercentLayer(ii, jj) = sum(lcBlock(:) == lcCodeList(1)) / lcBlockN;
        desertPercentLayer(ii, jj) = sum(ismember(lcBlock(:), desertCodeList)) / lcBlockN;
        waterPercentLayer(ii, jj) = sum(lcBlock(:) == lcCodeList(2)) / lcBlockN;
        glacierPercentLayer(ii, jj) = sum(lcBlock(:) == lcCodeList(3)) / lcBlockN;
        buildingPercentLayer(ii, jj) = sum(lcBlock(:) == lcCodeList(4)) / lcBlockN;
        lcPercents = [otherPercentLayer(ii, jj); desertPercentLayer(ii, jj); ...
            waterPercentLayer(ii, jj); glacierPercentLayer(ii, jj); buildingPercentLayer(ii, jj)];
        
        % Treat the pixel as pure if it holds one of the land cover more than 60% of the pixel area.
        % If any of the land cover doesn't exceed 60% in one AMSRE pixel, treat it as a mixed pixel.
        majorLc = max(lcPercents);
        if length(majorLc) == 1 && majorLc > lcThreshold  % one land cover exceed 60% of the area.
            maxIndex = find(lcPercents == majorLc);
            if maxIndex == 2  % desert.
                zonesLcLayer(ii, jj) = mode(lcBlock(:));  % The most frequent values.
            elseif maxIndex >= 3 % water, glacier, building.
                zonesLcLayer(ii, jj) = zonesLcLayer(ii, jj) + mode(lcBlock(:));
            end
        else  % mixed pixel.
            zonesLcLayer(ii, jj) = zonesLcLayer(ii, jj) + lcCodeList(6);
        end
        
        % Get desert code Layer with the resolution same as the AMSER BT image.
        if desertPercentLayer(ii, jj) > 0
            desertCodePixelCount = zeros(desertCodeListN, 1);
            for n = 1 : desertCodeListN
                desertCodePixelCount(n) = sum(lcBlock(:) == desertCodeList(n));
            end
            maxCountDesertCode = max(desertCodePixelCount);
            maxCountDesertCode = maxCountDesertCode(1);
            desertCodeLayer(ii, jj) = desertCodeList(desertCodePixelCount == maxCountDesertCode);
        end
    end
end

% Read the daily MODIS snow image and blend it into zonal landcover image, and fix the land cover
%   percent layers.
outLcPercentArrayMatPath = [matlabVariablePath, '3_LandCoverZones_', dayNight{flg1}, '_', ...
    testYears{flg7}, '.mat'];
if ~exist(outLcPercentArrayMatPath, 'file')
    [zonesLcArray, otherPercentArray, desertPercentArray, glacierPercentArray, waterPercentArray,...
        snowPercentArray, buildingPercentArray] = deal(zeros(btLstRowN, btLstColN, yearDateListN));
    % snowZonesLayer indicating the zone that snow is located in.
    snowZonesLayer = zonesLayer + lcCodeList(5);
    parfor i = 1 : yearDateListN
        % Techniques handeling the issues in parallel computing.
        tempZonesLcLayer = zonesLcLayer;
        tempSnowZonesLayer = snowZonesLayer;
        tempModisSnowList = modisSnowList;
        
        % Get index of the day in which the MODIS Snow images is avaliable.
        [~, snowDayIndex] = ismember(yearDateList{i}, modisSnowDateList);
        if snowDayIndex == 0, continue; end
        disp(['Read Snow image: ', modisSnowDateList{snowDayIndex}]);
        
        % Read MODIS Snow image.
        [modisSnowLayer, ~] = geotiffread([modisSnowYearPath, tempModisSnowList{snowDayIndex}]);
        modisSnowLayer = single(modisSnowLayer);
        % Get the snow percent layer.
        snowPercentLayer = zeros(btLstRowN, btLstColN);
        for ii = 1 : btLstRowN
            for jj = 1 : btLstColN
                % Skip the Nodata pixel in zonesLcLayer, that is the region out of China.
                if isnan(tempZonesLcLayer(ii, jj))
                    continue;
                end
                % Locate each Snow block.
                snowBlockTopRow = snowStartTopSideRow + snowBlockRowN * (ii - 1);
                snowBlockBottomRow = snowStartBottomSideCol + snowBlockRowN * (ii - 1);
                snowBlockLeftCol = snowStartLeftSideCol + snowBlockColN * (jj - 1);
                snowBlockRightCol = snowStartRightSideCol + snowBlockColN * (jj - 1);
                if snowBlockRightCol > snowColN || snowBlockBottomRow > snowRowN
                    continue;
                end
                % Snow block calculation.
                modisSnowBlock = modisSnowLayer(snowBlockTopRow : snowBlockBottomRow, ...
                    snowBlockLeftCol : snowBlockRightCol);
                modisSnowPixelIndex = (modisSnowBlock >= 0) & (modisSnowBlock <= 100);
                modisSnowPixelN = sum(modisSnowPixelIndex(:));
                snowPercentLayer(ii, jj) = sum(modisSnowBlock(modisSnowPixelIndex)) ./ ...
                    (100 * modisSnowPixelN);
            end
        end
        % Suppose there is no snow cover in nodata pixel.
        snowPercentLayer(isnan(snowPercentLayer)) = 0; 
        snowPercentArray(:, :, i) = snowPercentLayer;
        
        % Blend the snow cover into the zone lc image only in China.
        pureSnowIndexLayer = (snowPercentLayer > snowThreshold) & ~isnan(tempZonesLcLayer);
        tempZonesLcLayer(pureSnowIndexLayer) = tempSnowZonesLayer(pureSnowIndexLayer);
        zonesLcArray(:, :, i) = tempZonesLcLayer;
        
        % Fix the land cover percent layers.
        otherPercentArray(:, :, i) = otherPercentLayer .* (1 - snowPercentLayer);
        desertPercentArray(:, :, i) = desertPercentLayer .* (1 - snowPercentLayer);
        waterPercentArray(:, :, i) = waterPercentLayer .* (1 - snowPercentLayer);
        glacierPercentArray(:, :, i) = glacierPercentLayer .* (1 - snowPercentLayer);
        buildingPercentArray(:, :, i) = buildingPercentLayer .* (1 - snowPercentLayer);
    end
    save(outLcPercentArrayMatPath, 'otherPercentArray', 'desertPercentArray',...
        'waterPercentArray', 'glacierPercentArray', 'buildingPercentArray',...
        'snowPercentArray', 'zonesLcArray');
else
    lia = ismember({'otherPercentArray', 'desertPercentArray','waterPercentArray', ...
        'glacierPercentArray', 'buildingPercentArray', 'snowPercentArray', 'zonesLcArray'}, who);
    if sum(lia) < 7
        load(outLcPercentArrayMatPath);
    end
end

% Reclassify the land cover and zonal image code.
%   [Step1_EasternDongbei: 1: 1 ~ 4, Step1_Huabei: 2: 5 ~ 8, Step1_Hua'nan: 3: 9 ~ 13]
%   [Step2_Xi'nan: 4: 14 ~ 17, Step2_EasternXibei: 5: 18 ~ 27, Step2_WesternDongbei: 6: 28 ~ 30]
%   [Step2_WesternXibei: 7: 31 ~ 53, Step3_QingzangPlateau: 8: 54 ~ 71].
zonesLcArray(zonesLcArray == 0) = nan;
zonesLcIdList = unique(zonesLcArray);
zonesLcIdList(isnan(zonesLcIdList)) = [];
fixedZonesLcArray = zonesLcArray;
for i = 1 : regionsN
    for j = 2 : lcCodeListN - 1
        lcIndexVector = (zonesLcIdList >= lcCodeList(j) + regionNodes(i)) & ...
            (zonesLcIdList < lcCodeList(j) + regionNodes(i+1));
        zonesLcIdList(lcIndexVector) = lcCodeList(j) + i;
        lcIndexArray = (zonesLcArray >= lcCodeList(j) + regionNodes(i)) & ...
            (zonesLcArray < lcCodeList(j) + regionNodes(i+1));
        fixedZonesLcArray(lcIndexArray) = lcCodeList(j) + i;
    end
end
zonesLcIdList = unique(zonesLcIdList);
zonesLcIdListN = numel(zonesLcIdList);
mixedLcN = sum(zonesLcIdList >= lcCodeList(end));
pureLcIdN = zonesLcIdListN - mixedLcN;

%% Regression and Export.
% Regression, statistics and export of images and figures by Zones.
outRegressionParametersMatPath = [matlabVariablePath, '4_Regression_Parameters_', ...
    regressModels{flg3}, dayNight{flg1}, '_', testYears{flg7}, '.mat'];
if ~exist(outRegressionParametersMatPath, 'file')    
    [rmseYearVector, nYearVector, r2YearVector] = deal(zeros(pureLcIdN, 1) * nan);
    [rmseSeasonArray, nSeasonArray, r2SeasonArray] = deal(zeros(pureLcIdN, 4) * nan);
    [rmseMonthArray, nMonthArray, r2MonthArray] = deal(zeros(pureLcIdN, 12) * nan);
    
    coefficientYearArray = zeros(pureLcIdN, variablesN);
    coefficientSeasonArray = zeros(pureLcIdN, variablesN, 4);
    coefficientMonthArray = zeros(pureLcIdN, variablesN, 12);
    if flg6 == 2
        [amsreLstMaskYearArray, amsreLstCnYearArray] = ...
            deal(zeros(btLstRowN, btLstColN, yearDateListN));
    end
    for n = 1 : pureLcIdN
        zoneName = ['Zone', num2str(zonesLcIdList(n))];
        disp(['Regress ', zoneName]);
        
        % Pick out the zonal AMSRE BT and MODIS LST images from the array of the whole year.
        zonesLcIndexArray = (fixedZonesLcArray == zonesLcIdList(n));
        
        modisLstMaskZoneYearArray = setnan(modisLstMaskYearArray, ~zonesLcIndexArray);
        amsreBt06hMaskZoneYearArray = setnan(amsreBt06hMaskYearArray, ~zonesLcIndexArray);
        amsreBt06vMaskZoneYearArray = setnan(amsreBt06vMaskYearArray, ~zonesLcIndexArray);
        amsreBt10hMaskZoneYearArray = setnan(amsreBt10hMaskYearArray, ~zonesLcIndexArray);
        amsreBt10vMaskZoneYearArray = setnan(amsreBt10vMaskYearArray, ~zonesLcIndexArray);
        amsreBt18hMaskZoneYearArray = setnan(amsreBt18hMaskYearArray, ~zonesLcIndexArray);
        amsreBt18vMaskZoneYearArray = setnan(amsreBt18vMaskYearArray, ~zonesLcIndexArray);
        amsreBt23hMaskZoneYearArray = setnan(amsreBt23hMaskYearArray, ~zonesLcIndexArray);
        amsreBt23vMaskZoneYearArray = setnan(amsreBt23vMaskYearArray, ~zonesLcIndexArray);
        amsreBt36hMaskZoneYearArray = setnan(amsreBt36hMaskYearArray, ~zonesLcIndexArray);
        amsreBt36vMaskZoneYearArray = setnan(amsreBt36vMaskYearArray, ~zonesLcIndexArray);
        amsreBt89hMaskZoneYearArray = setnan(amsreBt89hMaskYearArray, ~zonesLcIndexArray);
        amsreBt89vMaskZoneYearArray = setnan(amsreBt89vMaskYearArray, ~zonesLcIndexArray);
        amsreBtQd1MaskZoneYearArray = ...
            (amsreBt36vMaskZoneYearArray - amsreBt18vMaskZoneYearArray) .^ 2;
        amsreBtQd2MaskZoneYearArray = ...
            (amsreBt36vMaskZoneYearArray - amsreBt23vMaskZoneYearArray) .^ 2;
        if flg6 == 2
            amsreBt06hCnZoneYearArray = setnan(amsreBt06hCnYearArray, ~zonesLcIndexArray);
            amsreBt06vCnZoneYearArray = setnan(amsreBt06vCnYearArray, ~zonesLcIndexArray);
            amsreBt10hCnZoneYearArray = setnan(amsreBt10hCnYearArray, ~zonesLcIndexArray);
            amsreBt10vCnZoneYearArray = setnan(amsreBt10vCnYearArray, ~zonesLcIndexArray);
            amsreBt18hCnZoneYearArray = setnan(amsreBt18hCnYearArray, ~zonesLcIndexArray);
            amsreBt18vCnZoneYearArray = setnan(amsreBt18vCnYearArray, ~zonesLcIndexArray);
            amsreBt23hCnZoneYearArray = setnan(amsreBt23hCnYearArray, ~zonesLcIndexArray);
            amsreBt23vCnZoneYearArray = setnan(amsreBt23vCnYearArray, ~zonesLcIndexArray);
            amsreBt36hCnZoneYearArray = setnan(amsreBt36hCnYearArray, ~zonesLcIndexArray);
            amsreBt36vCnZoneYearArray = setnan(amsreBt36vCnYearArray, ~zonesLcIndexArray);
            amsreBt89hCnZoneYearArray = setnan(amsreBt89hCnYearArray, ~zonesLcIndexArray);
            amsreBt89vCnZoneYearArray = setnan(amsreBt89vCnYearArray, ~zonesLcIndexArray);
            amsreBtQd1CnZoneYearArray = (amsreBt36vCnZoneYearArray - amsreBt18vCnZoneYearArray).^2;
            amsreBtQd2CnZoneYearArray = (amsreBt36vCnZoneYearArray - amsreBt23vCnZoneYearArray).^2;
        end
        
        % ------------------------------------------------------------------------------------------
        % Regress the zonal AMSRE BT and MODIS LST by year, and get the retreived AMSRE LST images.
        if flg5 == 1 || flg5 == 4
            valueIndex = find(~isnan(modisLstMaskZoneYearArray));
            modisLstMaskZoneYearVector = modisLstMaskZoneYearArray(valueIndex);
            amsreBt06hMaskZoneYearVector = amsreBt06hMaskZoneYearArray(valueIndex);
            amsreBt06vMaskZoneYearVector = amsreBt06vMaskZoneYearArray(valueIndex);
            amsreBt10hMaskZoneYearVector = amsreBt10hMaskZoneYearArray(valueIndex);
            amsreBt10vMaskZoneYearVector = amsreBt10vMaskZoneYearArray(valueIndex);
            amsreBt18hMaskZoneYearVector = amsreBt18hMaskZoneYearArray(valueIndex);
            amsreBt18vMaskZoneYearVector = amsreBt18vMaskZoneYearArray(valueIndex);
            amsreBt23hMaskZoneYearVector = amsreBt23hMaskZoneYearArray(valueIndex);
            amsreBt23vMaskZoneYearVector = amsreBt23vMaskZoneYearArray(valueIndex);
            amsreBt36hMaskZoneYearVector = amsreBt36hMaskZoneYearArray(valueIndex);
            amsreBt36vMaskZoneYearVector = amsreBt36vMaskZoneYearArray(valueIndex);
            amsreBt89hMaskZoneYearVector = amsreBt89hMaskZoneYearArray(valueIndex);
            amsreBt89vMaskZoneYearVector = amsreBt89vMaskZoneYearArray(valueIndex);
            amsreBtQd1MaskZoneYearVector = amsreBtQd1MaskZoneYearArray(valueIndex);
            amsreBtQd2MaskZoneYearVector = amsreBtQd2MaskZoneYearArray(valueIndex);               
            if flg5 == 1 && flg3 == 2
                amsreBtsMaskZoneYearRecords = [amsreBt06hMaskZoneYearVector, ...
                    amsreBt06vMaskZoneYearVector, amsreBt10hMaskZoneYearVector, ...
                    amsreBt10vMaskZoneYearVector, amsreBt18hMaskZoneYearVector, ...
                    amsreBt18vMaskZoneYearVector, amsreBt23hMaskZoneYearVector, ...
                    amsreBt23vMaskZoneYearVector, amsreBt36hMaskZoneYearVector, ...
                    amsreBt36vMaskZoneYearVector, amsreBt89hMaskZoneYearVector, ...
                    amsreBt89vMaskZoneYearVector, amsreBtQd1MaskZoneYearVector, ...
                    amsreBtQd2MaskZoneYearVector];
            else
                amsreBtsMaskZoneYearRecords = [amsreBt06hMaskZoneYearVector, ...
                    amsreBt06vMaskZoneYearVector, amsreBt10hMaskZoneYearVector, ...
                    amsreBt10vMaskZoneYearVector, amsreBt18hMaskZoneYearVector, ...
                    amsreBt18vMaskZoneYearVector, amsreBt23hMaskZoneYearVector, ...
                    amsreBt23vMaskZoneYearVector, amsreBt36hMaskZoneYearVector, ...
                    amsreBt36vMaskZoneYearVector, amsreBt89hMaskZoneYearVector, ...
                    amsreBt89vMaskZoneYearVector];
            end
            
            modisLstMaskZoneYearVectorN = length(modisLstMaskZoneYearVector);
            if modisLstMaskZoneYearVectorN >= 2
                mdl = stepwiselm(amsreBtsMaskZoneYearRecords, modisLstMaskZoneYearVector, ...
                    'constant', 'Lower', 'constant', 'Upper', 'linear', 'NSteps', 30);
                rmse1 = mdl.RMSE;
				rSquared = mdl.Rsquared.Ordinary;
                amsreLstMaskZoneYearVector = mdl.Fitted;
                inModel = mdl.VariableInfo.InModel;
                variableIndex = find(inModel == 1);
                pYear = mdl.Coefficients.Estimate;
                pYear2 = zeros(variablesN, 1);
                pYear2(variableIndex) = pYear(2:end);
                pYear2(end) = pYear(1);
                coefficientYearArray(n, :) = pYear2;
                if flg6 == 2
                    amsreLstMaskZoneYearVector2 = zeros(btLstRowN * btLstColN * yearDateListN, 1);
                    amsreLstMaskZoneYearVector2(valueIndex) = amsreLstMaskZoneYearVector;
                    amsreLstMaskZoneYearArray = reshape(amsreLstMaskZoneYearVector2, ...
                        btLstRowN, btLstColN, yearDateListN);
                    if flg5 == 1 && flg3 == 2
                        amsreLstCnZoneYearArray = pYear2(15) + ...
                            pYear2(1) .* amsreBt06hCnZoneYearArray + ...
                            pYear2(2) .* amsreBt06vCnZoneYearArray + ...
                            pYear2(3) .* amsreBt10hCnZoneYearArray + ...
                            pYear2(4) .* amsreBt10vCnZoneYearArray + ...
                            pYear2(5) .* amsreBt18hCnZoneYearArray + ...
                            pYear2(6) .* amsreBt18vCnZoneYearArray + ...
                            pYear2(7) .* amsreBt23hCnZoneYearArray + ...
                            pYear2(8) .* amsreBt23vCnZoneYearArray + ...
                            pYear2(9) .* amsreBt36hCnZoneYearArray + ...
                            pYear2(10) .* amsreBt36vCnZoneYearArray + ...
                            pYear2(11) .* amsreBt89hCnZoneYearArray + ...
                            pYear2(12) .* amsreBt89vCnZoneYearArray + ...
                            pYear2(13) .* amsreBtQd1CnZoneYearArray + ...
                            pYear2(14) .* amsreBtQd2CnZoneYearArray;
                    else
                        amsreLstCnZoneYearArray = pYear2(13) + ...
                            pYear2(1) .* amsreBt06hCnZoneYearArray + ...
                            pYear2(2) .* amsreBt06vCnZoneYearArray + ...
                            pYear2(3) .* amsreBt10hCnZoneYearArray + ...
                            pYear2(4) .* amsreBt10vCnZoneYearArray + ...
                            pYear2(5) .* amsreBt18hCnZoneYearArray + ...
                            pYear2(6) .* amsreBt18vCnZoneYearArray + ...
                            pYear2(7) .* amsreBt23hCnZoneYearArray + ...
                            pYear2(8) .* amsreBt23vCnZoneYearArray + ...
                            pYear2(9) .* amsreBt36hCnZoneYearArray + ...
                            pYear2(10) .* amsreBt36vCnZoneYearArray + ...
                            pYear2(11) .* amsreBt89hCnZoneYearArray + ...
                            pYear2(12) .* amsreBt89vCnZoneYearArray;
                    end
                    amsreLstCnZoneYearArray(isnan(amsreLstCnZoneYearArray)) = 0;
                end
            else
                rmse1 = nan;
                amsreLstMaskZoneYearVector = zeros(modisLstMaskZoneYearVectorN, 1) * nan;
                if flg6 == 2
                    amsreLstMaskZoneYearArray = zeros(btLstRowN, btLstColN, yearDateListN);
                    amsreLstCnZoneYearArray = amsreBt89vCnZoneYearArray * 0;
                end
            end
            if flg6 == 2
                amsreLstMaskYearArray = amsreLstMaskYearArray + amsreLstMaskZoneYearArray;
                amsreLstCnYearArray = amsreLstCnYearArray + amsreLstCnZoneYearArray;
            end
            rmseYearVector(n) = rmse1;
            nYearVector(n) = modisLstMaskZoneYearVectorN;
            r2YearVector(n) = rSquared;
            
            % Export the scatterplot of AMSRE BT(LST) vs MODIS LST by year.
            if flg4 == 1
                % Create the path used for restoring the zonal scatterplot.
                outFigureYearZonePath = [outFigureYearPath, zoneName, '\',];
                if ~exist(outFigureYearZonePath, 'dir')
                    mkdir(outFigureYearZonePath);
                end
                
                outZoneOriginalScatterYearPath = [outFigureYearZonePath, zoneName, ...
                    '_Original_', testYears{1}, '.png'];
                if ~exist(outZoneOriginalScatterYearPath, 'file')
                    f1 = btLstScatter(amsreBt89vMaskZoneYearVector, modisLstMaskZoneYearVector, ...
                        zoneName, testYears{1});
                    print(f1, outZoneOriginalScatterYearPath, '-dpng');
                end
                outZoneRegressedScatterYearPath = [outFigureYearZonePath, zoneName, ...
                    '_Regressed_', testYears{1}, '.png'];
                if ~exist(outZoneRegressedScatterYearPath, 'file')
                    f2 = lstScatter(amsreLstMaskZoneYearVector, modisLstMaskZoneYearVector, ...
                        zoneName, testYears{1}, rmse1);
                    print(f2, outZoneRegressedScatterYearPath, '-dpng');
                end
                close all;
            end
        end
        
        % ------------------------------------------------------------------------------------------
        % Regress and scatter the LST and BT by season, and get the retreived AMSRE LST images.
        if flg5 == 2 || flg5 == 4
            for i = 1 : 4
                % Pick out zonal LST & BT images of the season from the 3D array of the whole year.
                seasonIndex = ismember(yearMonthList, seasonMonths{i});
                seasonIndexN = sum(seasonIndex);
                
                modisLstMaskZoneSeasonArray = modisLstMaskZoneYearArray(:, :, seasonIndex);
                amsreBt06hMaskZoneSeasonArray = amsreBt06hMaskZoneYearArray(:, :, seasonIndex);
                amsreBt06vMaskZoneSeasonArray = amsreBt06vMaskZoneYearArray(:, :, seasonIndex);
                amsreBt10hMaskZoneSeasonArray = amsreBt10hMaskZoneYearArray(:, :, seasonIndex);
                amsreBt10vMaskZoneSeasonArray = amsreBt10vMaskZoneYearArray(:, :, seasonIndex);
                amsreBt18hMaskZoneSeasonArray = amsreBt18hMaskZoneYearArray(:, :, seasonIndex);
                amsreBt18vMaskZoneSeasonArray = amsreBt18vMaskZoneYearArray(:, :, seasonIndex);
                amsreBt23hMaskZoneSeasonArray = amsreBt23hMaskZoneYearArray(:, :, seasonIndex);
                amsreBt23vMaskZoneSeasonArray = amsreBt23vMaskZoneYearArray(:, :, seasonIndex);
                amsreBt36hMaskZoneSeasonArray = amsreBt36hMaskZoneYearArray(:, :, seasonIndex);
                amsreBt36vMaskZoneSeasonArray = amsreBt36vMaskZoneYearArray(:, :, seasonIndex);
                amsreBt89hMaskZoneSeasonArray = amsreBt89hMaskZoneYearArray(:, :, seasonIndex);
                amsreBt89vMaskZoneSeasonArray = amsreBt89vMaskZoneYearArray(:, :, seasonIndex);
                amsreBtQd1MaskZoneSeasonArray = amsreBtQd1MaskZoneYearArray(:, :, seasonIndex);
                amsreBtQd2MaskZoneSeasonArray = amsreBtQd2MaskZoneYearArray(:, :, seasonIndex);
                if flg6 == 2
                    amsreBt06hCnZoneSeasonArray = amsreBt06hCnZoneYearArray(:, :, seasonIndex);
                    amsreBt06vCnZoneSeasonArray = amsreBt06vCnZoneYearArray(:, :, seasonIndex);
                    amsreBt10hCnZoneSeasonArray = amsreBt10hCnZoneYearArray(:, :, seasonIndex);
                    amsreBt10vCnZoneSeasonArray = amsreBt10vCnZoneYearArray(:, :, seasonIndex);
                    amsreBt18hCnZoneSeasonArray = amsreBt18hCnZoneYearArray(:, :, seasonIndex);
                    amsreBt18vCnZoneSeasonArray = amsreBt18vCnZoneYearArray(:, :, seasonIndex);
                    amsreBt23hCnZoneSeasonArray = amsreBt23hCnZoneYearArray(:, :, seasonIndex);
                    amsreBt23vCnZoneSeasonArray = amsreBt23vCnZoneYearArray(:, :, seasonIndex);
                    amsreBt36hCnZoneSeasonArray = amsreBt36hCnZoneYearArray(:, :, seasonIndex);
                    amsreBt36vCnZoneSeasonArray = amsreBt36vCnZoneYearArray(:, :, seasonIndex);
                    amsreBt89hCnZoneSeasonArray = amsreBt89hCnZoneYearArray(:, :, seasonIndex);
                    amsreBt89vCnZoneSeasonArray = amsreBt89vCnZoneYearArray(:, :, seasonIndex);
                    amsreBtQd1CnZoneSeasonArray = amsreBtQd1CnZoneYearArray(:, :, seasonIndex);
                    amsreBtQd2CnZoneSeasonArray = amsreBtQd2CnZoneYearArray(:, :, seasonIndex);
                end
                
                % Regress the zonal BT & LST by season, and get the retreived AMSRE LST images.
                valueIndex = find(~isnan(modisLstMaskZoneSeasonArray));
                modisLstMaskZoneSeasonVector = modisLstMaskZoneSeasonArray(valueIndex);
                amsreBt06hMaskZoneSeasonVector = amsreBt06hMaskZoneSeasonArray(valueIndex);
                amsreBt06vMaskZoneSeasonVector = amsreBt06vMaskZoneSeasonArray(valueIndex);
                amsreBt10hMaskZoneSeasonVector = amsreBt10hMaskZoneSeasonArray(valueIndex);
                amsreBt10vMaskZoneSeasonVector = amsreBt10vMaskZoneSeasonArray(valueIndex);
                amsreBt18hMaskZoneSeasonVector = amsreBt18hMaskZoneSeasonArray(valueIndex);
                amsreBt18vMaskZoneSeasonVector = amsreBt18vMaskZoneSeasonArray(valueIndex);
                amsreBt23hMaskZoneSeasonVector = amsreBt23hMaskZoneSeasonArray(valueIndex);
                amsreBt23vMaskZoneSeasonVector = amsreBt23vMaskZoneSeasonArray(valueIndex);
                amsreBt36hMaskZoneSeasonVector = amsreBt36hMaskZoneSeasonArray(valueIndex);
                amsreBt36vMaskZoneSeasonVector = amsreBt36vMaskZoneSeasonArray(valueIndex);
                amsreBt89hMaskZoneSeasonVector = amsreBt89hMaskZoneSeasonArray(valueIndex);
                amsreBt89vMaskZoneSeasonVector = amsreBt89vMaskZoneSeasonArray(valueIndex);
                amsreBtQd1MaskZoneSeasonVector = amsreBtQd1MaskZoneSeasonArray(valueIndex);
                amsreBtQd2MaskZoneSeasonVector = amsreBtQd2MaskZoneSeasonArray(valueIndex);                 
                if flg5 == 2 && flg3 == 2
                    amsreBtsMaskZoneSeasonRecords = [amsreBt06hMaskZoneSeasonVector, ...
                        amsreBt06vMaskZoneSeasonVector, amsreBt10hMaskZoneSeasonVector, ...
                        amsreBt10vMaskZoneSeasonVector, amsreBt18hMaskZoneSeasonVector, ...
                        amsreBt18vMaskZoneSeasonVector, amsreBt23hMaskZoneSeasonVector, ...
                        amsreBt23vMaskZoneSeasonVector, amsreBt36hMaskZoneSeasonVector, ...
                        amsreBt36vMaskZoneSeasonVector, amsreBt89hMaskZoneSeasonVector, ...
                        amsreBt89vMaskZoneSeasonVector, amsreBtQd1MaskZoneSeasonVector, ...
                        amsreBtQd2MaskZoneSeasonVector];
                else
                    amsreBtsMaskZoneSeasonRecords = [amsreBt06hMaskZoneSeasonVector, ...
                        amsreBt06vMaskZoneSeasonVector, amsreBt10hMaskZoneSeasonVector, ...
                        amsreBt10vMaskZoneSeasonVector, amsreBt18hMaskZoneSeasonVector, ...
                        amsreBt18vMaskZoneSeasonVector, amsreBt23hMaskZoneSeasonVector, ...
                        amsreBt23vMaskZoneSeasonVector, amsreBt36hMaskZoneSeasonVector, ...
                        amsreBt36vMaskZoneSeasonVector, amsreBt89hMaskZoneSeasonVector, ...
                        amsreBt89vMaskZoneSeasonVector];
                end
                
                modisLstMaskZoneSeasonVectorN = length(modisLstMaskZoneSeasonVector);
                if modisLstMaskZoneSeasonVectorN >= 2
                    mdl = stepwiselm(amsreBtsMaskZoneSeasonRecords, modisLstMaskZoneSeasonVector,...
                        'constant', 'Lower', 'constant', 'Upper', 'linear', 'NSteps', 30);
                    rmse1 = mdl.RMSE;
                    rSquared = mdl.Rsquared.Ordinary;
                    amsreLstMaskZoneSeasonVector = mdl.Fitted;
                    inModel = mdl.VariableInfo.InModel;
                    variableIndex = find(inModel == 1);
                    pYear = mdl.Coefficients.Estimate;
                    pYear2 = zeros(variablesN, 1);
                    pYear2(variableIndex) = pYear(2:end);
                    pYear2(end) = pYear(1);
                    coefficientSeasonArray(n, :, i) = pYear2;
                    if flg6 == 2
                        btLstSeasonArrayN = btLstRowN * btLstColN * seasonIndexN;
                        amsreLstMaskZoneSeasonVector2 = zeros(btLstSeasonArrayN, 1);
                        amsreLstMaskZoneSeasonVector2(valueIndex) = amsreLstMaskZoneSeasonVector;
                        amsreLstMaskZoneSeasonArray = reshape(amsreLstMaskZoneSeasonVector2, ...
                            btLstRowN, btLstColN, seasonIndexN);
                        if flg5 == 2 && flg3 == 2
                            amsreLstCnZoneSeasonArray = pYear2(15) + ...
                                pYear2(1) .* amsreBt06hCnZoneSeasonArray + ...
                                pYear2(2) .* amsreBt06vCnZoneSeasonArray + ...
                                pYear2(3) .* amsreBt10hCnZoneSeasonArray + ...
                                pYear2(4) .* amsreBt10vCnZoneSeasonArray + ...
                                pYear2(5) .* amsreBt18hCnZoneSeasonArray + ...
                                pYear2(6) .* amsreBt18vCnZoneSeasonArray + ...
                                pYear2(7) .* amsreBt23hCnZoneSeasonArray + ...
                                pYear2(8) .* amsreBt23vCnZoneSeasonArray + ...
                                pYear2(9) .* amsreBt36hCnZoneSeasonArray + ...
                                pYear2(10) .* amsreBt36vCnZoneSeasonArray + ...
                                pYear2(11) .* amsreBt89hCnZoneSeasonArray + ...
                                pYear2(12) .* amsreBt89vCnZoneSeasonArray + ...
                                pYear2(13) .* amsreBtQd1CnZoneSeasonArray + ...
                                pYear2(14) .* amsreBtQd2CnZoneSeasonArray;
                        else
                            amsreLstCnZoneSeasonArray = pYear2(13) + ...
                                pYear2(1) .* amsreBt06hCnZoneSeasonArray + ...
                                pYear2(2) .* amsreBt06vCnZoneSeasonArray + ...
                                pYear2(3) .* amsreBt10hCnZoneSeasonArray + ...
                                pYear2(4) .* amsreBt10vCnZoneSeasonArray + ...
                                pYear2(5) .* amsreBt18hCnZoneSeasonArray + ...
                                pYear2(6) .* amsreBt18vCnZoneSeasonArray + ...
                                pYear2(7) .* amsreBt23hCnZoneSeasonArray + ...
                                pYear2(8) .* amsreBt23vCnZoneSeasonArray + ...
                                pYear2(9) .* amsreBt36hCnZoneSeasonArray + ...
                                pYear2(10) .* amsreBt36vCnZoneSeasonArray + ...
                                pYear2(11) .* amsreBt89hCnZoneSeasonArray + ...
                                pYear2(12) .* amsreBt89vCnZoneSeasonArray;
                        end
                        amsreLstCnZoneSeasonArray(isnan(amsreLstCnZoneSeasonArray)) = 0;
                    end
                else
                    rmse1 = nan;
                    amsreLstMaskZoneSeasonVector = zeros(modisLstMaskZoneSeasonVectorN, 1) * nan;
                    if flg6 == 2
                        amsreLstMaskZoneSeasonArray = zeros(btLstRowN, btLstColN, seasonIndexN);
                        amsreLstCnZoneSeasonArray = amsreBt89vCnZoneSeasonArray * 0;
                    end
                end
                if flg6 == 2
                    zonesLcSeasonIndexArray = zonesLcIndexArray(:, :, seasonIndex);
                    zonesLcSeasonIndexArray2 = false(btLstRowN, btLstColN, yearDateListN);
                    zonesLcSeasonIndexArray2(:, :, seasonIndex) = zonesLcSeasonIndexArray;
                    amsreLstMaskYearArray(zonesLcSeasonIndexArray2) = ...
                        amsreLstMaskZoneSeasonArray(zonesLcSeasonIndexArray);
                    amsreLstCnYearArray(zonesLcSeasonIndexArray2) = ...
                        amsreLstCnZoneSeasonArray(zonesLcSeasonIndexArray);
                end
                rmseSeasonArray(n, i) = rmse1;
                nSeasonArray(n, i) = modisLstMaskZoneSeasonVectorN;
                r2SeasonArray(n, i) = rSquared;
                
                % Export the scatterplot of AMSRE BT(LST) vs MODIS LST by season.
                if flg4 == 1
                    % Create the path used for restoring the zonal scatterplot.
                    outFigureYearZonePath = [outFigureYearPath, zoneName, '\',];
                    if ~exist(outFigureYearZonePath, 'dir')
                        mkdir(outFigureYearZonePath);
                    end
                
                    yearSeason = [testYears{1}, ' Season', num2str(i)];
                    outZoneOriginalScatterSeasonPath = [outFigureYearZonePath, zoneName, ...
                        '_Original_', yearSeason, '.png'];
                    if ~exist(outZoneOriginalScatterSeasonPath, 'file')
                        f1 = btLstScatter(amsreBt89vMaskZoneSeasonVector,  ...
                            modisLstMaskZoneSeasonVector, zoneName, yearSeason);
                        print(f1, outZoneOriginalScatterSeasonPath, '-dpng');
                    end
                    outZoneRegressedScatterSeasonPath = [outFigureYearZonePath, zoneName, ...
                        '_Regressed_', yearSeason, '.png'];
                    if ~exist(outZoneRegressedScatterSeasonPath, 'file')
                        f2 = lstScatter(amsreLstMaskZoneSeasonVector, ...
                            modisLstMaskZoneSeasonVector, zoneName, yearSeason, rmse1);
                        print(f2, outZoneRegressedScatterSeasonPath, '-dpng');
                    end
                    close all;
                end
            end
        end
        
        % ------------------------------------------------------------------------------------------
        % Regress and scatter the LST and BT by month, and get the retreived AMSRE LST images.
        if flg5 == 3 || flg5 == 4
            for i = 1 : 12
                % Pick out the zonal LST & BT images of the month from the 3D array of the year.
                monthIndex = ismember(yearMonthList, i);
                monthIndexN = sum(monthIndex);
                
                modisLstMaskZoneMonthArray = modisLstMaskZoneYearArray(:, :, monthIndex);
                amsreBt06hMaskZoneMonthArray = amsreBt06hMaskZoneYearArray(:, :, monthIndex);
                amsreBt06vMaskZoneMonthArray = amsreBt06vMaskZoneYearArray(:, :, monthIndex);
                amsreBt10hMaskZoneMonthArray = amsreBt10hMaskZoneYearArray(:, :, monthIndex);
                amsreBt10vMaskZoneMonthArray = amsreBt10vMaskZoneYearArray(:, :, monthIndex);
                amsreBt18hMaskZoneMonthArray = amsreBt18hMaskZoneYearArray(:, :, monthIndex);
                amsreBt18vMaskZoneMonthArray = amsreBt18vMaskZoneYearArray(:, :, monthIndex);
                amsreBt23hMaskZoneMonthArray = amsreBt23hMaskZoneYearArray(:, :, monthIndex);
                amsreBt23vMaskZoneMonthArray = amsreBt23vMaskZoneYearArray(:, :, monthIndex);
                amsreBt36hMaskZoneMonthArray = amsreBt36hMaskZoneYearArray(:, :, monthIndex);
                amsreBt36vMaskZoneMonthArray = amsreBt36vMaskZoneYearArray(:, :, monthIndex);
                amsreBt89hMaskZoneMonthArray = amsreBt89hMaskZoneYearArray(:, :, monthIndex);
                amsreBt89vMaskZoneMonthArray = amsreBt89vMaskZoneYearArray(:, :, monthIndex);
                amsreBtQd1MaskZoneMonthArray = amsreBtQd1MaskZoneYearArray(:, :, monthIndex);
                amsreBtQd2MaskZoneMonthArray = amsreBtQd2MaskZoneYearArray(:, :, monthIndex);
                if flg6 == 2
                    amsreBt06hCnZoneMonthArray = amsreBt06hCnZoneYearArray(:, :, monthIndex);
                    amsreBt06vCnZoneMonthArray = amsreBt06vCnZoneYearArray(:, :, monthIndex);
                    amsreBt10hCnZoneMonthArray = amsreBt10hCnZoneYearArray(:, :, monthIndex);
                    amsreBt10vCnZoneMonthArray = amsreBt10vCnZoneYearArray(:, :, monthIndex);
                    amsreBt18hCnZoneMonthArray = amsreBt18hCnZoneYearArray(:, :, monthIndex);
                    amsreBt18vCnZoneMonthArray = amsreBt18vCnZoneYearArray(:, :, monthIndex);
                    amsreBt23hCnZoneMonthArray = amsreBt23hCnZoneYearArray(:, :, monthIndex);
                    amsreBt23vCnZoneMonthArray = amsreBt23vCnZoneYearArray(:, :, monthIndex);
                    amsreBt36hCnZoneMonthArray = amsreBt36hCnZoneYearArray(:, :, monthIndex);
                    amsreBt36vCnZoneMonthArray = amsreBt36vCnZoneYearArray(:, :, monthIndex);
                    amsreBt89hCnZoneMonthArray = amsreBt89hCnZoneYearArray(:, :, monthIndex);
                    amsreBt89vCnZoneMonthArray = amsreBt89vCnZoneYearArray(:, :, monthIndex);
                    amsreBtQd1CnZoneMonthArray = amsreBtQd1CnZoneYearArray(:, :, monthIndex);
                    amsreBtQd2CnZoneMonthArray = amsreBtQd2CnZoneYearArray(:, :, monthIndex);
                end
                
                % Regress the zonal BT & LST by month, and get the retreived AMSRE LST images.
                valueIndex = find(~isnan(modisLstMaskZoneMonthArray));
                modisLstMaskZoneMonthVector = modisLstMaskZoneMonthArray(valueIndex);
                amsreBt06hMaskZoneMonthVector = amsreBt06hMaskZoneMonthArray(valueIndex);
                amsreBt06vMaskZoneMonthVector = amsreBt06vMaskZoneMonthArray(valueIndex);
                amsreBt10hMaskZoneMonthVector = amsreBt10hMaskZoneMonthArray(valueIndex);
                amsreBt10vMaskZoneMonthVector = amsreBt10vMaskZoneMonthArray(valueIndex);
                amsreBt18hMaskZoneMonthVector = amsreBt18hMaskZoneMonthArray(valueIndex);
                amsreBt18vMaskZoneMonthVector = amsreBt18vMaskZoneMonthArray(valueIndex);
                amsreBt23hMaskZoneMonthVector = amsreBt23hMaskZoneMonthArray(valueIndex);
                amsreBt23vMaskZoneMonthVector = amsreBt23vMaskZoneMonthArray(valueIndex);
                amsreBt36hMaskZoneMonthVector = amsreBt36hMaskZoneMonthArray(valueIndex);
                amsreBt36vMaskZoneMonthVector = amsreBt36vMaskZoneMonthArray(valueIndex);
                amsreBt89hMaskZoneMonthVector = amsreBt89hMaskZoneMonthArray(valueIndex);
                amsreBt89vMaskZoneMonthVector = amsreBt89vMaskZoneMonthArray(valueIndex);
                amsreBtQd1MaskZoneMonthVector = amsreBtQd1MaskZoneMonthArray(valueIndex);
                amsreBtQd2MaskZoneMonthVector = amsreBtQd2MaskZoneMonthArray(valueIndex);
                if flg3 == 2
                    amsreBtsMaskZoneMonthRecords = [amsreBt06hMaskZoneMonthVector, ...
                        amsreBt06vMaskZoneMonthVector, amsreBt10hMaskZoneMonthVector, ...
                        amsreBt10vMaskZoneMonthVector, amsreBt18hMaskZoneMonthVector, ...
                        amsreBt18vMaskZoneMonthVector, amsreBt23hMaskZoneMonthVector, ...
                        amsreBt23vMaskZoneMonthVector, amsreBt36hMaskZoneMonthVector, ...
                        amsreBt36vMaskZoneMonthVector, amsreBt89hMaskZoneMonthVector, ...
                        amsreBt89vMaskZoneMonthVector, amsreBtQd1MaskZoneMonthVector, ...
                        amsreBtQd2MaskZoneMonthVector];
                else
                    amsreBtsMaskZoneMonthRecords = [amsreBt06hMaskZoneMonthVector, ...
                        amsreBt06vMaskZoneMonthVector, amsreBt10hMaskZoneMonthVector, ...
                        amsreBt10vMaskZoneMonthVector, amsreBt18hMaskZoneMonthVector, ...
                        amsreBt18vMaskZoneMonthVector, amsreBt23hMaskZoneMonthVector, ...
                        amsreBt23vMaskZoneMonthVector, amsreBt36hMaskZoneMonthVector, ...
                        amsreBt36vMaskZoneMonthVector, amsreBt89hMaskZoneMonthVector, ...
                        amsreBt89vMaskZoneMonthVector];
                end
                
                modisLstMaskZoneMonthVectorN = length(modisLstMaskZoneMonthVector);
                if modisLstMaskZoneMonthVectorN >= 2
                    mdl = stepwiselm(amsreBtsMaskZoneMonthRecords, modisLstMaskZoneMonthVector, ...
                        'constant', 'Lower', 'constant', 'Upper', 'linear', 'NSteps', 30);
                    rmse1 = mdl.RMSE;
                    rSquared = mdl.Rsquared.Ordinary;
                    amsreLstMaskZoneMonthVector = mdl.Fitted;
                    inModel = mdl.VariableInfo.InModel;
                    variableIndex = find(inModel == 1);
                    pYear = mdl.Coefficients.Estimate;
                    pYear2 = zeros(variablesN, 1);
                    pYear2(variableIndex) = pYear(2:end);
                    pYear2(end) = pYear(1);
                    coefficientMonthArray(n, :, i) = pYear2;
                    if flg6 == 2
                        btLstMonthArrayN = btLstRowN * btLstColN * monthIndexN;
                        amsreLstMaskZoneMonthVector2 = zeros(btLstMonthArrayN, 1);
                        amsreLstMaskZoneMonthVector2(valueIndex) = amsreLstMaskZoneMonthVector;
                        amsreLstMaskZoneMonthArray = reshape(amsreLstMaskZoneMonthVector2, ...
                            btLstRowN, btLstColN, monthIndexN);
                        if flg3 == 1
                            amsreLstCnZoneMonthArray = pYear2(13) + ...
                                pYear2(1) .* amsreBt06hCnZoneMonthArray + ...
                                pYear2(2) .* amsreBt06vCnZoneMonthArray + ...
                                pYear2(3) .* amsreBt10hCnZoneMonthArray + ...
                                pYear2(4) .* amsreBt10vCnZoneMonthArray + ...
                                pYear2(5) .* amsreBt18hCnZoneMonthArray + ...
                                pYear2(6) .* amsreBt18vCnZoneMonthArray + ...
                                pYear2(7) .* amsreBt23hCnZoneMonthArray + ...
                                pYear2(8) .* amsreBt23vCnZoneMonthArray + ...
                                pYear2(9) .* amsreBt36hCnZoneMonthArray + ...
                                pYear2(10) .* amsreBt36vCnZoneMonthArray + ...
                                pYear2(11) .* amsreBt89hCnZoneMonthArray + ...
                                pYear2(12) .* amsreBt89vCnZoneMonthArray;
                        elseif flg3 == 2
                            amsreLstCnZoneMonthArray = pYear2(15) + ...
                                pYear2(1) .* amsreBt06hCnZoneMonthArray + ...
                                pYear2(2) .* amsreBt06vCnZoneMonthArray + ...
                                pYear2(3) .* amsreBt10hCnZoneMonthArray + ...
                                pYear2(4) .* amsreBt10vCnZoneMonthArray + ...
                                pYear2(5) .* amsreBt18hCnZoneMonthArray + ...
                                pYear2(6) .* amsreBt18vCnZoneMonthArray + ...
                                pYear2(7) .* amsreBt23hCnZoneMonthArray + ...
                                pYear2(8) .* amsreBt23vCnZoneMonthArray + ...
                                pYear2(9) .* amsreBt36hCnZoneMonthArray + ...
                                pYear2(10) .* amsreBt36vCnZoneMonthArray + ...
                                pYear2(11) .* amsreBt89hCnZoneMonthArray + ...
                                pYear2(12) .* amsreBt89vCnZoneMonthArray + ...
                                pYear2(13) .* amsreBtQd1CnZoneMonthArray + ...
                                pYear2(14) .* amsreBtQd2CnZoneMonthArray;
                        end
                        amsreLstCnZoneMonthArray(isnan(amsreLstCnZoneMonthArray)) = 0;
                    end
                else
                    rmse1 = nan;
                    amsreLstMaskZoneMonthVector = zeros(modisLstMaskZoneMonthVectorN, 1) * nan;
                    if flg6 == 2
                        amsreLstMaskZoneMonthArray = zeros(btLstRowN, btLstColN, monthIndexN);
                        amsreLstCnZoneMonthArray = amsreBt89vCnZoneMonthArray * 0;
                    end
                end
                if flg6 == 2
                    zonesLcMonthIndexArray = zonesLcIndexArray(:, :, monthIndex);
                    zonesLcMonthIndexArray2 = false(btLstRowN, btLstColN, yearDateListN);
                    zonesLcMonthIndexArray2(:, :, monthIndex) = zonesLcMonthIndexArray;
                    amsreLstMaskYearArray(zonesLcMonthIndexArray2) = ...
                        amsreLstMaskZoneMonthArray(zonesLcMonthIndexArray);
                    amsreLstCnYearArray(zonesLcMonthIndexArray2) = ...
                        amsreLstCnZoneMonthArray(zonesLcMonthIndexArray);
                end
                rmseMonthArray(n, i) = rmse1;
                nMonthArray(n, i) = modisLstMaskZoneMonthVectorN;
                r2MonthArray(n, i) = rSquared;
                
                % Export the scatterplot of AMSRE BT(LST) vs MODIS LST by month.
                if flg4 == 1
                    % Create the path used for restoring the zonal scatterplot.
                    outFigureYearZonePath = [outFigureYearPath, zoneName, '\',];
                    if ~exist(outFigureYearZonePath, 'dir')
                        mkdir(outFigureYearZonePath);
                    end
                    
                    yearMonth = [testYears{1}, num2str(i, '%02d')];
                    outZoneOriginalScatterMonthPath = [outFigureYearZonePath, zoneName, ...
                        '_Original_', yearMonth, '.png'];
                    if ~exist(outZoneOriginalScatterMonthPath, 'file')
                        f1 = btLstScatter(amsreBt89vMaskZoneMonthVector, ...
                            modisLstMaskZoneMonthVector, zoneName, yearMonth);
                        print(f1, outZoneOriginalScatterMonthPath, '-dpng');
                    end
                    outZoneRegressedScatterMonthPath = [outFigureYearZonePath, zoneName, ...
                        '_Regressed_', yearMonth, '.png'];
                    if ~exist(outZoneRegressedScatterMonthPath, 'file')
                        f2 = lstScatter(amsreLstMaskZoneMonthVector, ...
                            modisLstMaskZoneMonthVector, zoneName, yearMonth, rmse1);
                        print(f2, outZoneRegressedScatterMonthPath, '-dpng');
                    end
                    close all;
                end
            end
        end
        
        % ------------------------------------------------------------------------------------------
        % Export the zonal daily AMMSRE BT and MODIS LST images, and scatterplots.
        if flg4 == 2  % flg4 should be set 0 and 1. 2 means doesn't export figures in any cases.
            % Create the path used for restoring the zonal scatterplot.
            outFigureYearZonePath = [outFigureYearPath, zoneName, '\',];
            if ~exist(outFigureYearZonePath, 'dir')
                mkdir(outFigureYearZonePath);
            end
            
            for i = 1 : yearDateListN
                lstBtDate = modisLstMaskDateList{i};
                disp(lstBtDate);
                
                modisLstMaskZoneLayer = modisLstMaskZoneYearArray(:, :, i);
                amsreBt89vMaskZoneLayer = amsreBt89vMaskZoneYearArray(:, :, i);
                
                % AMSRE BT images.
                amsreBtMaskDayPath = [amsreBtMaskYearPath,amsreBtMaskDayFolderList{i},'\'];
                amsreBtMaskChannelList = dir([amsreBtMaskDayPath,'AMSRE_D*',orbits{flg1},'*.tif']);
                amsreBtChannelName2List = {amsreBtMaskChannelList.name}';
                
                outAmsreBtDayFolderPath = [outAmsreBtYearPath, amsreBtMaskDayFolderList{i}];
                if ~exist(outAmsreBtDayFolderPath, 'dir')
                    mkdir(outAmsreBtDayFolderPath);
                end
                outAmsreZoneImagePath = [outAmsreBtDayFolderPath, '\', zoneName, '_', ...
                    amsreBtChannelName2List{12}];
                if ~exist(outAmsreZoneImagePath, 'file')
                    geotiffwrite(outAmsreZoneImagePath, amsreBt89vMaskZoneLayer, btLstReference);
                end
                
                % MODIS LST images.
                outModisZoneFolderPath = [outModisLstYearPath, zoneName];
                if ~exist(outModisZoneFolderPath, 'dir')
                    mkdir(outModisZoneFolderPath);
                end
                outModisZoneImagePath = [outModisZoneFolderPath, '\', zoneName, '_', ...
                    modisLstMaskDayList{i}];
                if ~exist(outModisZoneImagePath, 'file')
                    geotiffwrite(outModisZoneImagePath, modisLstMaskZoneLayer, btLstReference);
                end
                
                % Export scatterplots.
                outZoneScatterPath = [outFigureYearZonePath, zoneName, '_', ...
                    replace(amsreBtChannelName2List{12}, '.tif', '.png')];
                if ~exist(outZoneScatterPath, 'file')
                    f1 = btLstScatter(amsreBt89vMaskZoneLayer, modisLstMaskZoneLayer, zoneName, ...
                        lstBtDate);
                    print(f1, outZoneScatterPath, '-dpng');
                end
                close all;
            end
        end
    end
    save(outRegressionParametersMatPath, 'rmseYearVector', 'rmseSeasonArray', 'rmseMonthArray', ...
        'nYearVector',  'nSeasonArray', 'nMonthArray', 'r2YearVector', 'r2SeasonArray', ...
        'r2MonthArray', 'coefficientYearArray', 'coefficientSeasonArray', 'coefficientMonthArray');
else
    lia = ismember({'rmseYearVector', 'rmseSeasonArray', 'rmseMonthArray', 'nYearVector',...
        'nSeasonArray', 'nMonthArray', 'r2YearVector', 'r2SeasonArray', 'r2MonthArray',...
        'coefficientYearArray', 'coefficientSeasonArray', 'coefficientMonthArray'}, who);
    if sum(lia) < 12
        load(outRegressionParametersMatPath);
    end
end

%% Optimized Model.
% The Optimized model is based on the Month model, in which the coefficients of a specific zone and
%   month are replaced by the coefficients of corresponding season model or year model when the
%   sample is too small to represent the trait of the specific zone and month.

% !!! The strategy here was verified to be not suitable, the replacement may cause huge error for
%   the retrieved AMSRE LST !!!

if flg5 == 4
    outfixedCoefficientMatPath = [matlabVariablePath, '5_Fixed_Coefficients_', ...
        regressModels{flg3}, dayNight{flg1}, '_', testYears{flg7}, '.mat'];
    if ~exist(outfixedCoefficientMatPath, 'file')
        % Average pixel number in the extent of a zone. The areas of ordinary zones(1 ~ 71) vary
        %   every day throughout the year due to the impact of snow cover, so the averaged pixel
        %   number in a zone is used to rule the model.
        zonesLcIdAreaCountVector = zeros(zonesLcIdListN, 1);
        for i = 1 : zonesLcIdListN
            zoneLcIdConut = sum(fixedZonesLcArray(:) == zonesLcIdList(i));
            zonesLcIdAreaCountVector(i) = zoneLcIdConut / modisSnowListN;
        end
        
        % Replace the coefficients of the zone in the month that has less than 1 time of pixel
        %   samples to the full extent pixel number by the coefficients of the season. If the pixel
        %   number of the season is less than 3 times to the full extent pixel number, apply the
        %   coefficients of the year. If the pixel number in a year is still rare, borrow the
        %   coefficients of similar land cover in other zones.
        fixedCoefficientMonthArray = coefficientMonthArray;
        fixedRmseMonthArray = rmseMonthArray;
        fixedR2MonthArray = r2MonthArray;
        fixedNMonthArray = nMonthArray;
        for i = 1 : pureLcIdN  % zones ID.
            for j = 1 : 12  % Month.
                monthSampleRatio = nMonthArray(i, j) / zonesLcIdAreaCountVector(i);
                if monthSampleRatio < 1
                    for n = 1 : 4  % Season.
                        if ismember(j, seasonMonths{n})
                            break;
                        end
                    end
                    seasonSampleRatio = nSeasonArray(i, n) / zonesLcIdAreaCountVector(i);
                    if seasonSampleRatio < 3
                        fixedCoefficientMonthArray(i, :, j) = coefficientYearArray(i, :);
                        fixedRmseMonthArray(i, j) = rmseYearVector(i);
                        fixedR2MonthArray(i, j) = r2YearVector(i);
                        fixedNMonthArray(i, j) = r2YearVector(i);
                    else
                        fixedCoefficientMonthArray(i, :, j) = coefficientSeasonArray(i, :, n);
                        fixedRmseMonthArray(i, j) = rmseSeasonArray(i, n);
                        fixedR2MonthArray(i, j) = r2SeasonArray(i, n);
                        fixedNMonthArray(i, j) = nSeasonArray(i, n);
                    end
                end
            end
        end

        % Handle the coefficients of rest several zones that has no pixels throughout the year.
        % Copy the coefficients in lc code 3006 to 3001, and 4008 to 4007.
        waterRegion1Index = find(zonesLcIdList == 3001);
        waterRegion6Index = find(zonesLcIdList == 3006);
        glacierRegion7Index = find(zonesLcIdList == 4007);
        glacierRegion8Index = find(zonesLcIdList == 4008);
        fixedCoefficientMonthArray(waterRegion1Index, :, :) = ...
            coefficientMonthArray(waterRegion6Index, :, :);
        fixedCoefficientMonthArray(glacierRegion7Index, :, :) = ...
            coefficientMonthArray(glacierRegion8Index, :, :);
        fixedRmseMonthArray(waterRegion1Index, :) = fixedRmseMonthArray(waterRegion6Index, :);
        fixedRmseMonthArray(glacierRegion7Index, :) = fixedRmseMonthArray(glacierRegion8Index, :);
        
        % Fix the coefficients in the zones 1 ~ 6008.
        % Zone 28, 29.
        fixedCoefficientMonthArray(zonesLcIdList == 1, :, [1, 2, 3, 11, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 6001, :, [1, 2, 3, 11, 12]);
        fixedCoefficientMonthArray(zonesLcIdList == 2, :, [1, 2, 3, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 6001, :, [1, 2, 3, 12]);
        fixedCoefficientMonthArray(zonesLcIdList == 3, :, [1, 2, 3, 4, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 6001, :, [1, 2, 3, 4, 12]);
        fixedCoefficientMonthArray(zonesLcIdList == 13, :, [4, 9]) = ...
            coefficientMonthArray(zonesLcIdList == 13, :, [3, 10]);
        fixedCoefficientMonthArray(zonesLcIdList == 17, :, 9) = ...
            coefficientMonthArray(zonesLcIdList == 17, :, 8);
        fixedCoefficientMonthArray(zonesLcIdList == 18, :, 1) = ...
            coefficientMonthArray(zonesLcIdList == 18, :, 2);
        fixedCoefficientMonthArray(zonesLcIdList == 27, :, 1) = ...
            coefficientMonthArray(zonesLcIdList == 27, :, 1);
        fixedCoefficientMonthArray(zonesLcIdList == 28, :, [1, 2, 3, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 6006, :, [1, 2, 3, 12]);
        fixedCoefficientMonthArray(zonesLcIdList == 29, :, [1, 2, 3, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 6006, :, [1, 2, 3, 12]);
        fixedCoefficientMonthArray(zonesLcIdList == 31, :, 2) = ...
            coefficientMonthArray(zonesLcIdList == 31, :, 2);
        fixedCoefficientMonthArray(zonesLcIdList == 33, :, [1, 2, 3, 4, 11, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 6007, :, [1, 2, 3, 4, 11, 12]);
        fixedCoefficientMonthArray(zonesLcIdList == 34, :, [1, 2, 3, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 6007, :, [1, 2, 3, 12]);
        fixedCoefficientMonthArray(zonesLcIdList == 37, :, [1, 2, 3, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 6007, :, [1, 2, 3, 12]);
        fixedCoefficientMonthArray(zonesLcIdList == 43, :, [1, 2, 3, 11, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 6007, :, [1, 2, 3, 11, 12]);
        fixedCoefficientMonthArray(zonesLcIdList == 45, :, [1, 2]) = ...
            coefficientMonthArray(zonesLcIdList == 44, :, [1, 2]);
        fixedCoefficientMonthArray(zonesLcIdList == 46, :, [1, 2]) = ...
            coefficientMonthArray(zonesLcIdList == 44, :, [1, 2]);
        fixedCoefficientMonthArray(zonesLcIdList == 47, :, [3, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 44, :, [4, 11]);
        fixedCoefficientMonthArray(zonesLcIdList == 47, :, [1, 2]) = ...
            coefficientMonthArray(zonesLcIdList == 6007, :, [1, 2]);
        fixedCoefficientMonthArray(zonesLcIdList == 50, :, [1, 2, 3, 11, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 32, :, [1, 2, 3, 11, 12]);
        fixedCoefficientMonthArray(zonesLcIdList == 65, :, 6) = ...
            coefficientMonthArray(zonesLcIdList == 65, :, 6);
        fixedCoefficientMonthArray(zonesLcIdList == 2000, :, [1, 2, 3, 12]) = ...
            coefficientMonthArray(zonesLcIdList == 6007, :, [1, 2, 3, 12]);
        
        save(outfixedCoefficientMatPath, 'fixedCoefficientMonthArray');
    else
        lia = ismember('fixedCoefficientMonthArray', who);
        if sum(lia) < 1
            load(outfixedCoefficientMatPath);
        end
    end
    
    % ----------------------------------------------------------------------------------------------
    % Retrieve AMSRE LST in pure zones(1:1:71, 1000:100:2000, 3000:1000:6000) using optimized model.
    outAmsreLstPureMatPath = [matlabVariablePath, '6_AMSRE_LST_PureLc_', regressModels{flg3}, ...
        dayNight{flg1}, '_', testYears{flg7}, '.mat'];
    if ~exist(outAmsreLstPureMatPath, 'file')
        amsreLstCnPureLcYearArray = zeros(btLstRowN, btLstColN, yearDateListN);
        % ------------------
%         fixedCoefficientMonthArray = coefficientMonthArray;
        % ------------------
        for n = 1 : pureLcIdN
            zoneName = ['Zone', num2str(zonesLcIdList(n))];
            disp(['Retrieve AMSRE LST in pure pixels: ', zoneName]);
            
            % Pick out the zonal AMSRE BT and MODIS LST images from the array of the whole year.
            zonesLcIndexArray = (fixedZonesLcArray == zonesLcIdList(n));
            
            amsreBt06hCnZoneYearArray = setnan(amsreBt06hCnYearArray, ~zonesLcIndexArray);
            amsreBt06vCnZoneYearArray = setnan(amsreBt06vCnYearArray, ~zonesLcIndexArray);
            amsreBt10hCnZoneYearArray = setnan(amsreBt10hCnYearArray, ~zonesLcIndexArray);
            amsreBt10vCnZoneYearArray = setnan(amsreBt10vCnYearArray, ~zonesLcIndexArray);
            amsreBt18hCnZoneYearArray = setnan(amsreBt18hCnYearArray, ~zonesLcIndexArray);
            amsreBt18vCnZoneYearArray = setnan(amsreBt18vCnYearArray, ~zonesLcIndexArray);
            amsreBt23hCnZoneYearArray = setnan(amsreBt23hCnYearArray, ~zonesLcIndexArray);
            amsreBt23vCnZoneYearArray = setnan(amsreBt23vCnYearArray, ~zonesLcIndexArray);
            amsreBt36hCnZoneYearArray = setnan(amsreBt36hCnYearArray, ~zonesLcIndexArray);
            amsreBt36vCnZoneYearArray = setnan(amsreBt36vCnYearArray, ~zonesLcIndexArray);
            amsreBt89hCnZoneYearArray = setnan(amsreBt89hCnYearArray, ~zonesLcIndexArray);
            amsreBt89vCnZoneYearArray = setnan(amsreBt89vCnYearArray, ~zonesLcIndexArray);
            amsreBtQd1CnZoneYearArray = (amsreBt36vCnZoneYearArray - amsreBt18vCnZoneYearArray).^2;
            amsreBtQd2CnZoneYearArray = (amsreBt36vCnZoneYearArray - amsreBt23vCnZoneYearArray).^2;

            for i = 1 : 12  % Month.
                % Pick out the zonal LST & BT images of the month from the 3D array of the year.
                monthIndex = ismember(modisLstMaskMonthList, i);
                monthIndexN = sum(monthIndex);
                
                amsreBt06hCnZoneMonthArray = amsreBt06hCnZoneYearArray(:, :, monthIndex);
                amsreBt06vCnZoneMonthArray = amsreBt06vCnZoneYearArray(:, :, monthIndex);
                amsreBt10hCnZoneMonthArray = amsreBt10hCnZoneYearArray(:, :, monthIndex);
                amsreBt10vCnZoneMonthArray = amsreBt10vCnZoneYearArray(:, :, monthIndex);
                amsreBt18hCnZoneMonthArray = amsreBt18hCnZoneYearArray(:, :, monthIndex);
                amsreBt18vCnZoneMonthArray = amsreBt18vCnZoneYearArray(:, :, monthIndex);
                amsreBt23hCnZoneMonthArray = amsreBt23hCnZoneYearArray(:, :, monthIndex);
                amsreBt23vCnZoneMonthArray = amsreBt23vCnZoneYearArray(:, :, monthIndex);
                amsreBt36hCnZoneMonthArray = amsreBt36hCnZoneYearArray(:, :, monthIndex);
                amsreBt36vCnZoneMonthArray = amsreBt36vCnZoneYearArray(:, :, monthIndex);
                amsreBt89hCnZoneMonthArray = amsreBt89hCnZoneYearArray(:, :, monthIndex);
                amsreBt89vCnZoneMonthArray = amsreBt89vCnZoneYearArray(:, :, monthIndex);
                amsreBtQd1CnZoneMonthArray = amsreBtQd1CnZoneYearArray(:, :, monthIndex);
                amsreBtQd2CnZoneMonthArray = amsreBtQd2CnZoneYearArray(:, :, monthIndex);
                
                pOptimized = fixedCoefficientMonthArray(n, :, i);
                if flg3 == 1
                    amsreLstCnZoneMonthArray = pOptimized(end) + ...
                        pOptimized(1) .* amsreBt06hCnZoneMonthArray + ...
                        pOptimized(2) .* amsreBt06vCnZoneMonthArray + ...
                        pOptimized(3) .* amsreBt10hCnZoneMonthArray + ...
                        pOptimized(4) .* amsreBt10vCnZoneMonthArray + ...
                        pOptimized(5) .* amsreBt18hCnZoneMonthArray + ...
                        pOptimized(6) .* amsreBt18vCnZoneMonthArray + ...
                        pOptimized(7) .* amsreBt23hCnZoneMonthArray + ...
                        pOptimized(8) .* amsreBt23vCnZoneMonthArray + ...
                        pOptimized(9) .* amsreBt36hCnZoneMonthArray + ...
                        pOptimized(10) .* amsreBt36vCnZoneMonthArray + ...
                        pOptimized(11) .* amsreBt89hCnZoneMonthArray + ...
                        pOptimized(12) .* amsreBt89vCnZoneMonthArray;
                elseif flg3 == 2
                    amsreLstCnZoneMonthArray = pOptimized(end) + ...
                        pOptimized(1) .* amsreBt06hCnZoneMonthArray + ...
                        pOptimized(2) .* amsreBt06vCnZoneMonthArray + ...
                        pOptimized(3) .* amsreBt10hCnZoneMonthArray + ...
                        pOptimized(4) .* amsreBt10vCnZoneMonthArray + ...
                        pOptimized(5) .* amsreBt18hCnZoneMonthArray + ...
                        pOptimized(6) .* amsreBt18vCnZoneMonthArray + ...
                        pOptimized(7) .* amsreBt23hCnZoneMonthArray + ...
                        pOptimized(8) .* amsreBt23vCnZoneMonthArray + ...
                        pOptimized(9) .* amsreBt36hCnZoneMonthArray + ...
                        pOptimized(10) .* amsreBt36vCnZoneMonthArray + ...
                        pOptimized(11) .* amsreBt89hCnZoneMonthArray + ...
                        pOptimized(12) .* amsreBt89vCnZoneMonthArray + ...
                        pOptimized(13) .* amsreBtQd1CnZoneMonthArray + ...
                        pOptimized(14) .* amsreBtQd2CnZoneMonthArray;
                end
                amsreLstCnZoneMonthArray(isnan(amsreLstCnZoneMonthArray)) = 0;
                
                zonesLcMonthIndexArray = zonesLcIndexArray(:, :, monthIndex);
                zonesLcMonthIndexArray2 = false(btLstRowN, btLstColN, yearDateListN);
                zonesLcMonthIndexArray2(:, :, monthIndex) = zonesLcMonthIndexArray;
                amsreLstCnPureLcYearArray(zonesLcMonthIndexArray2) = ...
                    amsreLstCnZoneMonthArray(zonesLcMonthIndexArray);
            end
        end
        save(outAmsreLstPureMatPath, 'amsreLstCnPureLcYearArray');
    else
        lia = ismember('amsreLstCnPureLcYearArray', who);
        if sum(lia) < 1
            load(outAmsreLstPureMatPath);
        end
    end

    % Retrieve the AMSRE LST in mixed Lc pixels.
    % In the first loop, the partial lst consisting of the land cover of other, water, glacier,
    %   building, and snow was combined. In the second loop the rest part of desert LST was
    %   appended. The partition model of desert is different from that of other land covers, they
    %   cannot easily handled in a single loop.
    outAmsreLstMixedMatPath = [matlabVariablePath, '7_AMSRE_LST_MixedLc_', regressModels{flg3}, ...
        dayNight{flg1}, '_', testYears{flg7}, '.mat'];
    if ~exist(outAmsreLstMixedMatPath, 'file')
        % First loop. Traverse (bian li) each ordinary zone (1 ~ 71) and find mixed pixels in this
        %   zone, then derive the partial compound temperature of the pixel by weighted averaging
        %   the other, water, glacier, building, and snow cover temperature.
        zonesIdList = unique(zonesLayer);
        zonesIdList(isnan(zonesIdList)) = [];
        zonesIdListN = length(zonesIdList);
        amsreLstCnMixPixelYearArray = zeros(btLstRowN, btLstColN, yearDateListN);
        for n = 1 : zonesIdListN
            % Find the mixed pixels throughout the entire year in each zone. For example: 1 + 10000
            %   means mixed pixels in zone 1.
            zoneId = zonesIdList(n);
            mixLcId = zoneId + lcCodeList(end);
            mixPixelIndexArray = (fixedZonesLcArray == mixLcId);
            disp(['Retrieve AMSRE LST in mixed pixels: ', num2str(mixLcId)]);
            
            % Skip the zone that has no mixed pixels.
            if sum(mixPixelIndexArray(:)) == 0
                continue;
            end
            
            % Get the complete AMSRE BT images of the zone in 12 channels throughout the year.
            amsreBt06hCnMixPixelYearArray = setnan(amsreBt06hCnYearArray, ~mixPixelIndexArray);
            amsreBt06vCnMixPixelYearArray = setnan(amsreBt06vCnYearArray, ~mixPixelIndexArray);
            amsreBt10hCnMixPixelYearArray = setnan(amsreBt10hCnYearArray, ~mixPixelIndexArray);
            amsreBt10vCnMixPixelYearArray = setnan(amsreBt10vCnYearArray, ~mixPixelIndexArray);
            amsreBt18hCnMixPixelYearArray = setnan(amsreBt18hCnYearArray, ~mixPixelIndexArray);
            amsreBt18vCnMixPixelYearArray = setnan(amsreBt18vCnYearArray, ~mixPixelIndexArray);
            amsreBt23hCnMixPixelYearArray = setnan(amsreBt23hCnYearArray, ~mixPixelIndexArray);
            amsreBt23vCnMixPixelYearArray = setnan(amsreBt23vCnYearArray, ~mixPixelIndexArray);
            amsreBt36hCnMixPixelYearArray = setnan(amsreBt36hCnYearArray, ~mixPixelIndexArray);
            amsreBt36vCnMixPixelYearArray = setnan(amsreBt36vCnYearArray, ~mixPixelIndexArray);
            amsreBt89hCnMixPixelYearArray = setnan(amsreBt89hCnYearArray, ~mixPixelIndexArray);
            amsreBt89vCnMixPixelYearArray = setnan(amsreBt89vCnYearArray, ~mixPixelIndexArray);
            amsreBtQd1CnMixPixelYearArray = ...
                (amsreBt36vCnMixPixelYearArray - amsreBt18vCnMixPixelYearArray) .^ 2;
            amsreBtQd2CnMixPixelYearArray = ...
                (amsreBt36vCnMixPixelYearArray - amsreBt23vCnMixPixelYearArray) .^ 2;
            
            % Get the percentage of each land cover in the zone throughout the year.
            otherPercentInMixPixelYearArray = setnan(otherPercentArray, ~mixPixelIndexArray);
            waterPercentInMixPixelYearArray = setnan(waterPercentArray, ~mixPixelIndexArray);
            glacierPercentInMixPixelYearArray = setnan(glacierPercentArray, ~mixPixelIndexArray);
            buildingPercentInMixPixelYearArray = setnan(buildingPercentArray, ~mixPixelIndexArray);
            snowPercentInMixPixelYearArray = setnan(snowPercentArray, ~mixPixelIndexArray);
            
            % Derive AMSRE BT in mixed pixels of each zone by month.
            for i = 1 : 12
                monthIndex = ismember(modisLstMaskMonthList, i);
                monthIndexN = sum(monthIndex);
                mixPixelMonthIndexArray = mixPixelIndexArray(:, :, monthIndex);
                % Skip the zone that has no mixed pixels.
                if sum(mixPixelMonthIndexArray(:)) == 0
                    continue;
                end
                % Get the complete AMSRE BT images of the zone in the month.
                amsreBt06hCnMixPixelMonthArray = amsreBt06hCnMixPixelYearArray(:, :, monthIndex);
                amsreBt06vCnMixPixelMonthArray = amsreBt06vCnMixPixelYearArray(:, :, monthIndex);
                amsreBt10hCnMixPixelMonthArray = amsreBt10hCnMixPixelYearArray(:, :, monthIndex);
                amsreBt10vCnMixPixelMonthArray = amsreBt10vCnMixPixelYearArray(:, :, monthIndex);
                amsreBt18hCnMixPixelMonthArray = amsreBt18hCnMixPixelYearArray(:, :, monthIndex);
                amsreBt18vCnMixPixelMonthArray = amsreBt18vCnMixPixelYearArray(:, :, monthIndex);
                amsreBt23hCnMixPixelMonthArray = amsreBt23hCnMixPixelYearArray(:, :, monthIndex);
                amsreBt23vCnMixPixelMonthArray = amsreBt23vCnMixPixelYearArray(:, :, monthIndex);
                amsreBt36hCnMixPixelMonthArray = amsreBt36hCnMixPixelYearArray(:, :, monthIndex);
                amsreBt36vCnMixPixelMonthArray = amsreBt36vCnMixPixelYearArray(:, :, monthIndex);
                amsreBt89hCnMixPixelMonthArray = amsreBt89hCnMixPixelYearArray(:, :, monthIndex);
                amsreBt89vCnMixPixelMonthArray = amsreBt89vCnMixPixelYearArray(:, :, monthIndex);
                amsreBtQd1CnMixPixelMonthArray = amsreBtQd1CnMixPixelYearArray(:, :, monthIndex);
                amsreBtQd2CnMixPixelMonthArray = amsreBtQd2CnMixPixelYearArray(:, :, monthIndex);
                
                % Get the percentage arrays of different land covers in the month.
                otherPercentInMixPixelMonthArray = ...
                    otherPercentInMixPixelYearArray(:, :, monthIndex);
                waterPercentInMixPixelMonthArray = ...
                    waterPercentInMixPixelYearArray(:, :, monthIndex);
                glacierPercentInMixPixelMonthArray = ...
                    glacierPercentInMixPixelYearArray(:, :, monthIndex);
                buildingPercentInMixPixelMonthArray = ...
                    buildingPercentInMixPixelYearArray(:, :, monthIndex);
                snowPercentMixedLcMonthArray = ...
                    snowPercentInMixPixelYearArray(:, :, monthIndex);
                
                % Find the regression coefficients of AMSRE LST in certain land cover.
                % regionNodes: [1, 5, 9, 14, 18, 28, 31, 57, 71+1]
                % regionIDs: [1, 2, 3, 4, 5, 6, 7, 8]
                % Convert zoneId to regionID and find the index of landcover in zonesLcIdList.
                for j = 1 : regionsN
                    if (regionNodes(j) <= zoneId) && (zoneId < regionNodes(j+1))
                        regionId = j;
                        break;
                    end
                end
                
                % If there exist the certain land cover in the zone, get its coefficients, otherwise
                %   set the coefficients as zeros.
                waterLcIdIndex = (zonesLcIdList == (lcCodeList(2) + regionId));
                glacierLcIdIndex = (zonesLcIdList == (lcCodeList(3) + regionId));
                buildingLcIdIndex = (zonesLcIdList == (lcCodeList(4) + regionId));
                snowLcIdIndex = (zonesLcIdList == (lcCodeList(5) + regionId));
                
                pOtherMonth = fixedCoefficientMonthArray(n, :, i);
                [pWaterMonth, pSnowMonth, pGlacierMonth] = deal(zeros(1, variablesN));
                if sum(waterLcIdIndex) == 1
                    pWaterMonth = fixedCoefficientMonthArray(waterLcIdIndex, :, i);
                end
                if sum(snowLcIdIndex) == 1
                    pSnowMonth = fixedCoefficientMonthArray(snowLcIdIndex, :, i);
                end
                if sum(glacierLcIdIndex) == 1
                    pGlacierMonth = fixedCoefficientMonthArray(glacierLcIdIndex, :, i);
                end
                % Fix the coefficients of land cover having few or no samples in the regions.
                pBuildingMonth = fixedCoefficientMonthArray(zonesLcIdList == 10, :, i);
                
                % Regression.
                % First combine the coefficients, then multiply it by AMSRE BTs in 12 channels.
                coefficientsAmsreBt = zeros(btLstRowN, btLstColN, monthIndexN, variablesN);
                for j = 1 : variablesN
                    coefficientsAmsreBt(:, :, :, j) = ...
                        pOtherMonth(j) .* otherPercentInMixPixelMonthArray + ...
                        pWaterMonth(j) .* waterPercentInMixPixelMonthArray + ...
                        pGlacierMonth(j) .* glacierPercentInMixPixelMonthArray + ...
                        pBuildingMonth(j) .* buildingPercentInMixPixelMonthArray + ...
                        pSnowMonth(j) .* snowPercentMixedLcMonthArray;
                end
                if flg3 == 1
                    amsrePartialLstCnMixedLcMonthArray = coefficientsAmsreBt(:, :, :, end) + ...
                        coefficientsAmsreBt(:, :, :, 1) .* amsreBt06hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 2) .* amsreBt06vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 3) .* amsreBt10hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 4) .* amsreBt10vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 5) .* amsreBt18hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 6) .* amsreBt18vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 7) .* amsreBt23hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 8) .* amsreBt23vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 9) .* amsreBt36hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 10) .* amsreBt36vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 11) .* amsreBt89hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 12) .* amsreBt89vCnMixPixelMonthArray;
                elseif flg3 == 2
                    amsrePartialLstCnMixedLcMonthArray = coefficientsAmsreBt(:, :, :, end) + ...
                        coefficientsAmsreBt(:, :, :, 1) .* amsreBt06hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 2) .* amsreBt06vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 3) .* amsreBt10hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 4) .* amsreBt10vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 5) .* amsreBt18hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 6) .* amsreBt18vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 7) .* amsreBt23hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 8) .* amsreBt23vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 9) .* amsreBt36hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 10) .* amsreBt36vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 11) .* amsreBt89hCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 12) .* amsreBt89vCnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 13) .* amsreBtQd1CnMixPixelMonthArray + ...
                        coefficientsAmsreBt(:, :, :, 14) .* amsreBtQd2CnMixPixelMonthArray;
                end
                amsrePartialLstCnMixedLcMonthArray(isnan(amsrePartialLstCnMixedLcMonthArray)) = 0;
                
                mixPixelMonthIndexArray2 = false(btLstRowN, btLstColN, yearDateListN);
                mixPixelMonthIndexArray2(:, :, monthIndex) = mixPixelMonthIndexArray;
                amsreLstCnMixPixelYearArray(mixPixelMonthIndexArray2) = ...
                    amsrePartialLstCnMixedLcMonthArray(mixPixelMonthIndexArray);
            end
        end
        save(outAmsreLstMixedMatPath, 'amsreLstCnMixPixelYearArray');
    else
        lia = ismember('amsreLstCnMixPixelYearArray', who);
        if sum(lia) < 1
            load(outAmsreLstMixedMatPath);
        end
    end

    % Second loop. Add the rest partical LST of desert into the partical LST derived from the first
    %   loop. Traverse each desert region (1000, 1100, ..., 2000)
    outAmsreLstDesertMatPath = [matlabVariablePath, '8_AMSRE_LST_Desert_', regressModels{flg3}, ...
        dayNight{flg1}, '_', testYears{flg7}, '.mat'];
    if ~exist(outAmsreLstDesertMatPath, 'file')
        % Get desertCodeArray used for indicating the desert area that the mixed pixels belong to. 
        desertCodeArray = repmat(desertCodeLayer, 1, 1, yearDateListN);
        mixedDesertIndexArray = (0 < desertPercentArray) & (desertPercentArray <= lcThreshold);
        mixedLcsIndexArray = (fixedZonesLcArray > lcCodeList(end));
        for n = 1 : desertCodeListN
            % Find the mixed pixels in each desert area.
            desertId = desertCodeList(n);
            disp(['Retrieve AMSRE LST in rest desert part: ', num2str(desertId)]);
            
            desertCodeIndexArray = (desertCodeArray == desertId);
            mixPixelInDesertIndexArray = ...
                mixedDesertIndexArray & mixedLcsIndexArray & desertCodeIndexArray;
            
            desertPercentDesertYearArray = setnan(desertPercentArray, ~mixPixelInDesertIndexArray);
            amsreBt06hCnDesertYearArray = setnan(amsreBt06hCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt06vCnDesertYearArray = setnan(amsreBt06vCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt10hCnDesertYearArray = setnan(amsreBt10hCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt10vCnDesertYearArray = setnan(amsreBt10vCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt18hCnDesertYearArray = setnan(amsreBt18hCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt18vCnDesertYearArray = setnan(amsreBt18vCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt23hCnDesertYearArray = setnan(amsreBt23hCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt23vCnDesertYearArray = setnan(amsreBt23vCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt36hCnDesertYearArray = setnan(amsreBt36hCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt36vCnDesertYearArray = setnan(amsreBt36vCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt89hCnDesertYearArray = setnan(amsreBt89hCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBt89vCnDesertYearArray = setnan(amsreBt89vCnYearArray, ~mixPixelInDesertIndexArray);
            amsreBtQd1CnDesertYearArray = ...
                (amsreBt36vCnDesertYearArray - amsreBt18vCnDesertYearArray) .^ 2;
            amsreBtQd2CnDesertYearArray = ...
                (amsreBt36vCnDesertYearArray - amsreBt23vCnDesertYearArray) .^ 2;
            
            for i = 1 : 12  % Month.
                monthIndex = ismember(yearMonthList, i);
                monthIndexN = sum(monthIndex);
                desertInMixPixelMonthIndexArray = mixPixelInDesertIndexArray(:, :, monthIndex);
                % Skip the month that has no mixed pixels in the desert.
                if sum(desertInMixPixelMonthIndexArray(:)) == 0
                    continue;
                end
                
                desertPercentDesertMonthArray = desertPercentDesertYearArray(:, :, monthIndex);
                amsreBt06hCnDesertMonthArray = amsreBt06hCnDesertYearArray(:, :, monthIndex);
                amsreBt06vCnDesertMonthArray = amsreBt06vCnDesertYearArray(:, :, monthIndex);
                amsreBt10hCnDesertMonthArray = amsreBt10hCnDesertYearArray(:, :, monthIndex);
                amsreBt10vCnDesertMonthArray = amsreBt10vCnDesertYearArray(:, :, monthIndex);
                amsreBt18hCnDesertMonthArray = amsreBt18hCnDesertYearArray(:, :, monthIndex);
                amsreBt18vCnDesertMonthArray = amsreBt18vCnDesertYearArray(:, :, monthIndex);
                amsreBt23hCnDesertMonthArray = amsreBt23hCnDesertYearArray(:, :, monthIndex);
                amsreBt23vCnDesertMonthArray = amsreBt23vCnDesertYearArray(:, :, monthIndex);
                amsreBt36hCnDesertMonthArray = amsreBt36hCnDesertYearArray(:, :, monthIndex);
                amsreBt36vCnDesertMonthArray = amsreBt36vCnDesertYearArray(:, :, monthIndex);
                amsreBt89hCnDesertMonthArray = amsreBt89hCnDesertYearArray(:, :, monthIndex);
                amsreBt89vCnDesertMonthArray = amsreBt89vCnDesertYearArray(:, :, monthIndex);
                amsreBtQd1CnDesertMonthArray = amsreBtQd1CnDesertYearArray(:, :, monthIndex);
                amsreBtQd2CnDesertMonthArray = amsreBtQd2CnDesertYearArray(:, :, monthIndex);
                
                % Find the regression coefficients of AMSRE LST in desert. Desert 1200 has no 
                %   coefficients, appoint its coefficients same as Desert 1100.
                pDesertMonth = zeros(1, variablesN);
                desertLcIdIndex = find(zonesLcIdList == desertId);
                if desertId == 1200
                    desertLcIdIndex = find(zonesLcIdList == 1100);
                end
                if ~isempty(desertLcIdIndex)
                    pDesertMonth = fixedCoefficientMonthArray(desertLcIdIndex, :, i);
                end
                
                % Regression.
                if flg3 == 1
                    amsreLstCnPartial2MonthArray = ...
                        desertPercentDesertMonthArray .* (pDesertMonth(end) + ...
                        pDesertMonth(1) .* amsreBt06hCnDesertMonthArray + ...
                        pDesertMonth(2) .* amsreBt06vCnDesertMonthArray + ...
                        pDesertMonth(3) .* amsreBt10hCnDesertMonthArray + ...
                        pDesertMonth(4) .* amsreBt10vCnDesertMonthArray + ...
                        pDesertMonth(5) .* amsreBt18hCnDesertMonthArray + ...
                        pDesertMonth(6) .* amsreBt18vCnDesertMonthArray + ...
                        pDesertMonth(7) .* amsreBt23hCnDesertMonthArray + ...
                        pDesertMonth(8) .* amsreBt23vCnDesertMonthArray + ...
                        pDesertMonth(9) .* amsreBt36hCnDesertMonthArray + ...
                        pDesertMonth(10) .* amsreBt36vCnDesertMonthArray + ...
                        pDesertMonth(11) .* amsreBt89hCnDesertMonthArray + ...
                        pDesertMonth(12) .* amsreBt89vCnDesertMonthArray);
                elseif flg3 == 2
                    amsreLstCnPartial2MonthArray = ...
                        desertPercentDesertMonthArray .* (pDesertMonth(end) + ...
                        pDesertMonth(1) .* amsreBt06hCnDesertMonthArray + ...
                        pDesertMonth(2) .* amsreBt06vCnDesertMonthArray + ...
                        pDesertMonth(3) .* amsreBt10hCnDesertMonthArray + ...
                        pDesertMonth(4) .* amsreBt10vCnDesertMonthArray + ...
                        pDesertMonth(5) .* amsreBt18hCnDesertMonthArray + ...
                        pDesertMonth(6) .* amsreBt18vCnDesertMonthArray + ...
                        pDesertMonth(7) .* amsreBt23hCnDesertMonthArray + ...
                        pDesertMonth(8) .* amsreBt23vCnDesertMonthArray + ...
                        pDesertMonth(9) .* amsreBt36hCnDesertMonthArray + ...
                        pDesertMonth(10) .* amsreBt36vCnDesertMonthArray + ...
                        pDesertMonth(11) .* amsreBt89hCnDesertMonthArray + ...
                        pDesertMonth(12) .* amsreBt89vCnDesertMonthArray + ...
                        pDesertMonth(13) .* amsreBtQd1CnDesertMonthArray + ...
                        pDesertMonth(14) .* amsreBtQd2CnDesertMonthArray);
                end
                amsreLstCnPartial2MonthArray(isnan(amsreLstCnPartial2MonthArray)) = 0;
                
                desertInMixPixelMonthIndexArray2 = false(btLstRowN, btLstColN, yearDateListN);
                desertInMixPixelMonthIndexArray2(:, :, monthIndex) = desertInMixPixelMonthIndexArray;
                amsreLstCnMixPixelYearArray(desertInMixPixelMonthIndexArray2) = ...
                    amsreLstCnMixPixelYearArray(desertInMixPixelMonthIndexArray2) + ...
                    amsreLstCnPartial2MonthArray(desertInMixPixelMonthIndexArray);
            end
        end
        save(outAmsreLstDesertMatPath, 'amsreLstCnMixPixelYearArray');
    else
        lia = ismember('amsreLstCnMixPixelYearArray', who);
        if sum(lia) < 1
            load(outAmsreLstDesertMatPath);
        end
    end
end

% Merge the pure pixels and mixed pixels.
outAmsreLstCnMatPath = [matlabVariablePath, '9_AMSRE_LST_CN_', regressModels{flg3}, ...
    dayNight{flg1}, '_', testYears{flg7}, '.mat'];
if ~exist(outAmsreLstCnMatPath, 'file')
    amsreLstCnYearArray = amsreLstCnPureLcYearArray + amsreLstCnMixPixelYearArray;
    amsreLstCnYearArray(amsreLstCnYearArray <= 0) = nan;
    save(outAmsreLstCnMatPath, 'amsreLstCnYearArray');
else
    lia = ismember('amsreLstCnYearArray', who);
    if sum(lia) < 1
        load(outAmsreLstCnMatPath);
    end
end

% Remove the outlier of retrieved AMSRE LST.


%% Patch the gaps in AMSRE LST images.
outAmsreLstFillMatPath = ...
    [matlabVariablePath, '10_AMSRE_LST_Filled_', regressModels{flg3}, testYears{flg7}, '.mat'];
if ~exist(outAmsreLstFillMatPath, 'file')
    % The days interval for interpolate the AMSRE LST in orbit gaps.
    interpDays = 9;
    halfInterval = floor(interpDays / 2);
    
    % Add the head and tail layers to the AMSRE LST array.
    amsreLstCnHeadArray = amsreLstCnYearArray(:, :, end-halfInterval+1 : end);
    amsreLstCnTailArray = amsreLstCnYearArray(:, :, 1 : halfInterval);
    amsreLstCnYearArray2 = cat(3, amsreLstCnHeadArray, amsreLstCnYearArray, amsreLstCnTailArray);
    
    % Add the head and tail layers to the date list.
    headYearDateList = yearDateList(end-halfInterval+1 : end);
    tailYearDateList = yearDateList(1 : halfInterval);
    yearDateList2 = [headYearDateList; yearDateList; tailYearDateList];
    
    % Fill the gaps in the daily AMSRE LST images by interpolating the timeseries of the pixels
    %   ahead and behind them.
    for i = (1 + halfInterval) : (size(amsreLstCnYearArray2, 3) - halfInterval)
        disp(yearDateList2{i});
        
        % Prepare the data for interpolation.
        amsreLstCnDayLayer = amsreLstCnYearArray2(:, :, i);
        amsreLstCnInterpDaysArray = amsreLstCnYearArray2(:, :, i - halfInterval : i + halfInterval);
        amsreLstCnDayGapIndexList = find(isnan(amsreLstCnDayLayer) & ~isnan(zonesLayer));
        amsreLstCnDayGapIndexN = length(amsreLstCnDayGapIndexList);
        
        % Interpolation.
        for j = 1 : amsreLstCnDayGapIndexN
            gapIndex = amsreLstCnDayGapIndexList(j);
            yVector = zeros(interpDays, 1) * nan;
            for n = 1 : interpDays
                amsreLstCnInterpDayLayer = amsreLstCnInterpDaysArray(:, :, n);
                yVector(n) = amsreLstCnInterpDayLayer(gapIndex);
            end
            xVector = 1 : interpDays;
            xVector(isnan(yVector)) = [];
            yVector(isnan(yVector)) = [];
            if length(xVector) >= 2
                amsreLstGapValue = interp1(xVector, yVector, ceil(interpDays / 2));
            elseif isempty(xVector)
                amsreLstGapValue = nan;
            else
                amsreLstGapValue = yVector;
            end
            amsreLstCnDayLayer(gapIndex) = amsreLstGapValue;
        end
        amsreLstCnYearArray2(:, :, i) = amsreLstCnDayLayer;
    end
    amsreLstCnFillYearArray = ...
        amsreLstCnYearArray2(:, :, 1 + halfInterval : size(amsreLstCnYearArray2, 3) - halfInterval);
    
    save(outAmsreLstFillMatPath, 'amsreLstCnFillYearArray');
    
else
    lia = ismember('amsreLstCnFillYearArray', who);
    if sum(lia) < 1
        load(outAmsreLstFillMatPath);
    end
end

% Export the merged and gap-filled AMSRE LST iamges and scatterplot between them and MODIS LSTs.
% amsreLstCnFillYearArray(amsreLstCnFillYearArray < 240) = nan;
% modisLstMaskYearArray(amsreLstCnFillYearArray < 240) = nan;
if flg4 == 1
    for i = 1 : yearDateListN
        % Daily AMSRE and MODIS LST scatterplots for the extent of China.
        testDate = yearDateList{i};
        outLstScatterYearPath = [outFigureYearPath, 'China_', testDate, '_Regressed_', ...
            testYears{flg7}, '.png'];
        if ~exist(outLstScatterYearPath, 'file')
            amsreLstCnLayer = amsreLstCnYearArray(:, :, i);
            modisLstMaskLayer = modisLstMaskYearArray(:, :, i);
            rmse1 = rmse(amsreLstCnLayer, modisLstMaskLayer);
            disp([testDate, ': ', num2str(rmse1)]);
            f = lstScatter(amsreLstCnLayer, modisLstMaskLayer, 'China', testDate, rmse1);
            print(f, outLstScatterYearPath, '-dpng');
            close all;
        end
        
        % Merged AMSRE LST for China.
        outAmsreLstImagePath = [outAmsreLstYearPath, 'AMSRE_LST_merged_', testDate, '_Day.tif'];
        if ~exist(outAmsreLstImagePath, 'file')
            amsreLstCnLayer = amsreLstCnYearArray(:, :, i);
            geotiffwrite(outAmsreLstImagePath, amsreLstCnLayer, btLstReference);
        end
        
%         % Gap-filled AMSRE LST for China.
%         outAmsreLstImagePath = [outAmsreLstYearPath, 'AMSRE_LST_filled_', testDate, '_Day.tif'];
%         if ~exist(outAmsreLstImagePath, 'file')
%             amsreLstCnFillLayer = amsreLstCnFillYearArray(:, :, i);
%             geotiffwrite(outAmsreLstImagePath, amsreLstCnFillLayer, btLstReference);
%         end
    end
    
    % AMSRE and MODIS LST scatterplots for the extent of China in the whole year.
    outLstScatterYearPath = [outFigureYearPath, 'China', '_Regressed_', testYears{flg7}, '.png'];
    if ~exist(outLstScatterYearPath, 'file')
        rmse1 = rmse(amsreLstCnYearArray, modisLstMaskYearArray);
        f = lstScatter(amsreLstCnYearArray, modisLstMaskYearArray, 'China', testYears{flg7}, rmse1);
        print(f, outLstScatterYearPath, '-dpng');
    end
end

% Export Scatterplot between AMSRE LST and MODIS LST based on the zonal strategy of Zhou.
zonesZhouFilePath = [dataRootPath, zonesFileNameList{2}];
[zonesZhouLayer, ~] = geotiffread(zonesZhouFilePath);
zonesZhouLayer = single(zonesZhouLayer);
zonesZhouLayer(ismember(zonesZhouLayer, [-128, 128])) = nan;
zonesZhouIdList = unique(zonesZhouLayer);
zonesZhouIdList(isnan(zonesZhouIdList)) = [];
zonesZhouIdListN = length(zonesZhouIdList);
for i = 1 : zonesZhouIdListN
    zonesIdIndexLayer = (zonesZhouLayer == zonesZhouIdList(i));
    zonesIdIndexArray = repmat(zonesIdIndexLayer, 1, 1, yearDateListN);
    
    amsreLstZoneFillYearArray = setnan(amsreLstCnFillYearArray, ~zonesIdIndexArray);
    modisLstMaskZoneYearArray = setnan(modisLstMaskYearArray, ~zonesIdIndexArray);
    
    for n = 1 : length(seasonMonths)
        seasonIndex = ismember(yearMonthList, seasonMonths{n});
        seasonIndexN = sum(seasonIndex);
        outLstScatterYearPath = [outFigureYearPath, ['Landcover', num2str(zonesZhouIdList(i))], ...
            ['_Season', num2str(n)], '_Regressed_', testYears{flg7}, '.png'];
        if ~exist(outLstScatterYearPath, 'file')
            amsreLstZoneFillSeasonArray = amsreLstZoneFillYearArray(:, :, seasonIndex);
            modisLstMaskZoneSeasonArray = modisLstMaskZoneYearArray(:, :, seasonIndex);
            rmseSeason = rmse(amsreLstZoneFillSeasonArray, modisLstMaskZoneSeasonArray);
            f = lstScatter(amsreLstZoneFillSeasonArray, modisLstMaskZoneSeasonArray, ...
                ['Landcover', num2str(zonesZhouIdList(i))], ...
                [testYears{flg7}, ' Season', num2str(n)], rmseSeason);
            print(f, outLstScatterYearPath, '-dpng');
            close all;
        end
    end
end

%% Histogram.
% rmseYearVector(rmseYearVector == 0) = nan;
% rmse2SeasonArr(rmse2SeasonArr == 0) = nan;
% rmse2MonthArr(rmse2MonthArr == 0) = nan;
% 
h1 = histRmse(rmseYearVector, 0:1:10, 'Year');
% 
% h2 = histRmse(rmse2SeasonArr(:, 1), 0:1:10, 'Spring');
% h3 = histRmse(rmse2SeasonArr(:, 2), 0:1:10, 'Summer');
% h4 = histRmse(rmse2SeasonArr(:, 3), 0:1:10, 'Altumn');
% h5 = histRmse(rmse2SeasonArr(:, 4), 0:1:10, 'Winter');

% h6 = histRmse(rmse2MonthArr(:, 1), 0:1:10, 'Jan');
% h7 = histRmse(rmse2MonthArr(:, 2), 0:1:10, 'Feb');
% h8 = histRmse(rmse2MonthArr(:, 3), 0:1:10, 'Mar');
% h9 = histRmse(rmse2MonthArr(:, 4), 0:1:10, 'Apr');
% h10 = histRmse(rmse2MonthArr(:, 5), 0:1:10, 'May');
% h11 = histRmse(rmse2MonthArr(:, 6), 0:1:10, 'Jun');
% h12 = histRmse(rmse2MonthArr(:, 7), 0:1:10, 'Jul');
% h13 = histRmse(rmse2MonthArr(:, 8), 0:1:10, 'Aug');
% h14 = histRmse(rmse2MonthArr(:, 9), 0:1:10, 'Sep');
% h15 = histRmse(rmse2MonthArr(:, 10), 0:1:10, 'Obt');
% h16 = histRmse(rmse2MonthArr(:, 11), 0:1:10, 'Nov');
% h17 = histRmse(rmse2MonthArr(:, 12), 0:1:10, 'Dec');


%% Functions.
% Set nan.
function array2 = setnan(array1, nanarray)
array2 = array1;
array2(nanarray) = nan;
end

% Scatter AMSRE BT and MODIS LST.
function f = btLstScatter(amsreBt, modisLst, zoneName, dateString)
amsreBt = amsreBt(~isnan(amsreBt));
modisLst = modisLst(~isnan(modisLst));

p = [nan, nan];
if ~isempty(amsreBt)
    p = polyfit(amsreBt, modisLst, 1);
end
btLstRmse = rmse(amsreBt, modisLst);
btLstR = corrcoef(amsreBt, modisLst, 'rows', 'complete');
if numel(btLstR) > 1
    btLstR = btLstR(1, 2);
end

x1 = min(amsreBt) - 8;
x2 = max(amsreBt) + 5;
y1 = x1*p(1) + p(2);
y2 = x2*p(1) + p(2);

f = figure; f.Visible = 'off';
plot(amsreBt, modisLst, '.k', [220, 360], [220, 360], 'r', [x1, x2], [y1, y2], 'b');
xlabel('AMSRE BT'); ylabel('MODIS LST');
title(['AMSRE 89V BT vs MODIS LST in ', zoneName, ' ', dateString]);

txt1 = ['Pixel Number: ', num2str(sum(~isnan(modisLst)))];
txt2 = ['R: ', num2str(btLstR, '%.3f')];
txt3 = ['RMSE: ', num2str(btLstRmse, '%.3f')];
txt4 = ['Formula: ', 'y = ', num2str(p(1), '%.3f'), 'x + ', num2str(p(2), '%.3f')];
text(0.6, 0.17, txt1, 'Units', 'normalized', 'FontSize', 12);
text(0.6, 0.11, txt2, 'Units', 'normalized', 'FontSize', 12);
text(0.6, 0.05, txt3, 'Units', 'normalized', 'FontSize', 12);
text(0.03, 0.95, txt4, 'Units', 'normalized', 'FontSize', 12);
end

% Scatter AMSRE LST and MODIS LST.
function f = lstScatter(amsreLst, modisLst, zoneName, dateString, rmse1)
valueIndex = ~isnan(amsreLst) & ~isnan(modisLst);
amsreLst = amsreLst(valueIndex);
modisLst = modisLst(valueIndex);

bias = sum(amsreLst-modisLst) / length(amsreLst);
mae1 = sum(abs(amsreLst-modisLst)) / length(amsreLst);

btLstR = corrcoef(amsreLst, modisLst, 'rows', 'complete');
if numel(btLstR) > 1
    btLstR = btLstR(1, 2);
end

f = figure; f.Visible = 'off';
plot(amsreLst, modisLst, '.k', [220, 360], [220, 360], 'r');
xlabel('AMSRE LST'); ylabel('MODIS LST');
title(['AMSRE LST vs MODIS LST in ', zoneName, ' ', dateString]);

txt1 = ['Pixel Number: ', num2str(sum(~isnan(modisLst)))];
txt2 = ['R: ', num2str(btLstR, '%.3f')];
txt3 = ['RMSE: ', num2str(rmse1, '%.3f')];
txt4 = ['Bias: ', num2str(bias, '%.3f')];
txt5 = ['MAE:', num2str(mae1, '%.3f')];
text(0.6, 0.29, txt1, 'Units', 'normalized', 'FontSize', 12);
text(0.6, 0.23, txt2, 'Units', 'normalized', 'FontSize', 12);
text(0.6, 0.17, txt3, 'Units', 'normalized', 'FontSize', 12);
text(0.6, 0.11, txt4, 'Units', 'normalized', 'FontSize', 12);
text(0.6, 0.05, txt5, 'Units', 'normalized', 'FontSize', 12);
end

% Histogram of RMSE.
function h = histRmse(rmseVector, edges, period)
figure;
h = histogram(rmseVector, edges);
ylim([0, 45]);
xlabel('RMSE (K)'); ylabel('Zone Count');
title(['RMSE of retrieved AMSRE LST in 2010 ', period]);
end
