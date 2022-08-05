%% 将AMSR2 H5格式的各通道亮温数据根据指定空间范围输出为TIF格式, 并保存为Mat格式.
% 虽然此脚本输出的影像范围包含了全球, 但实际并不使用.

%% 标记和预设参数.
% 指定输出范围的标记. 1表示全球, 2表示中国(包含高亚洲地区).
flg1 = 1;
% 指定白天和晚上的标记. 1表示白天(升轨), 2表示晚上(降轨).
flg2 = 1;

% 输出AMSR2影像的经纬度范围和空间参考.
latLim = {[-90, 90], [17, 55]}; lonLim = {[0, 360], [66, 136]};  % {[World], [CN]}
latLim = latLim{flg1}; lonLim = lonLim{flg1}; cellsize = 0.1;

startLon = (lonLim(1) - 0) / cellsize + 1; endLon = (lonLim(2) - 0) / cellsize;
startLat = (90 - latLim(2)) / cellsize + 1; endLat = (90 - latLim(1)) / cellsize;
amsr2Ref = georefcells(latLim, lonLim, cellsize, cellsize, 'ColumnsStartFrom', 'north');

% AMSR2数据的通道与极化.
channelList = ["10", "18", "23", "36", "6", "7", "89"];
polarize = ["H", "V"];
[cMatrix, pMatrix] = meshgrid(channelList, polarize);
cpList = cellstr(reshape(cMatrix + pMatrix, [], 1));
cpListN = length(cpList);

% AMSR2数据的范围, 轨道和昼夜.
extent = {'World', 'CN'};
orbit = {'A', 'D'};
daynight = {'Day', 'Night'};

extent = extent{flg1};
orbit = orbit{flg2};
daynight = daynight{flg2};

% 有数据的年份月份列表(时间区间：2012/07/02-2020/12/31).
startDate = [2012, 07, 02]; endDate = [2020, 12, 31];
dateAllList = cellstr(datetime(startDate) : datetime(endDate), 'yyyyMMdd')';
dateAllListN = length(dateAllList);
yearNumList = startDate(1) : endDate(1);
yearNumListN = length(yearNumList);

%% 目录.
% 根目录.
rootPath = 'E:\AMSR_MODIS_Fusion\Data';

% AMSR2 BT数据H5格式总目录.
amsr2H5Path = fullfile(rootPath, 'AMSR2_1_H5');

% AMSR2 BT数据TIF, Mat格式总目录.
amsr2TifPath = fullfile(rootPath, sprintf('AMSR2_2_%s_TIF', extent));
if ~exist(amsr2TifPath, 'dir')
    mkdir(amsr2TifPath)
end

amsr2MatPath = fullfile(rootPath, sprintf('AMSR2_2_%s_Matlab', extent));
if ~exist(amsr2MatPath, 'dir')
    mkdir(amsr2MatPath);
end

%% 将AMSR2 H5格式转为TIF格式.
% 读取AMSR2 H5波段文件夹列表.
disp('将AMSR2 H5格式转为TIF格式.');
amsr2BandFolderList = dir(fullfile(amsr2H5Path, 'L3*'));
amsr2BandFolderList = {amsr2BandFolderList.name}';
for i = 1 : length(amsr2BandFolderList)
    amsr2BandFolderName = amsr2BandFolderList{i};
    amsr2H5BandFolderPath = fullfile(amsr2H5Path, amsr2BandFolderName);
    amsr2TifBandFolderPath = fullfile(amsr2TifPath, amsr2BandFolderName);
    if ~exist(amsr2TifBandFolderPath, 'dir')  % 创建AMSR2 TIF波段文件夹.
        mkdir(amsr2TifBandFolderPath)
    end

    % 读取AMSR2 H5年份文件夹列表.
    amsr2YearFolderList = dir(amsr2H5BandFolderPath);
    amsr2YearFolderList = {amsr2YearFolderList(3:end).name}';
    for j = 1 : length(amsr2YearFolderList)
        amsr2YearFolderName = amsr2YearFolderList{j};
        amsr2H5YearFolderPath = fullfile(amsr2H5BandFolderPath, amsr2YearFolderName);
        amsr2TifYearFolderPath = fullfile(amsr2TifBandFolderPath, amsr2YearFolderName);
        if ~exist(amsr2TifYearFolderPath, 'dir')  % 创建AMSR2 TIF年份文件夹.
            mkdir(amsr2TifYearFolderPath)
        end

        % 读取AMSR2 H5月份文件夹列表.
        amsr2MonthFolderList = dir(amsr2H5YearFolderPath);
        amsr2MonthFolderList = {amsr2MonthFolderList(3:end).name}';
        for k = 1 : length(amsr2MonthFolderList)
            amsr2MonthFolderName = amsr2MonthFolderList{k};
            amsr2H5MonthFolderPath = fullfile(amsr2H5YearFolderPath, amsr2MonthFolderName);
            amsr2TifMonthFolderPath = fullfile(amsr2TifYearFolderPath, amsr2MonthFolderName);
            if ~exist(amsr2TifMonthFolderPath, 'dir')  % 创建AMSR2 TIF月份文件夹.
                mkdir(amsr2TifMonthFolderPath)
            end
            fprintf('转换%s年%s月AMSR2 BT %s数据的格式.\n', amsr2YearFolderName, ...
                amsr2MonthFolderName, amsr2BandFolderName(6:10))

            % 读取AMSR2 H5每日文件列表.
            amsr2H5DailyList = dir(fullfile(amsr2H5MonthFolderPath, '*.h5'));
            amsr2H5DailyList = {amsr2H5DailyList.name}';
            for m = 1 : length(amsr2H5DailyList)
                amsr2H5DailyName = amsr2H5DailyList{m};
                amsr2H5DailyPath = fullfile(amsr2H5MonthFolderPath, amsr2H5DailyName);

                % 保存H5格式中的极化亮温数据为TIF格式.
                for n = 1 : length(polarize)
                    replacedStr = sprintf('_Bt%s.tif', polarize(n));
                    amsr2BtTifDailyName = replace(amsr2H5DailyName, '.h5', replacedStr);
                    amsr2BtTifDailyPath = fullfile(amsr2TifMonthFolderPath, amsr2BtTifDailyName);
                    if ~exist(amsr2BtTifDailyPath, 'file')
                        layerName = sprintf('/Brightness Temperature (%s)', polarize(n));
                        try
                            btArray = h5read(amsr2H5DailyPath, layerName)';
                            btArray = btArray(startLat : endLat, startLon : endLon);
                            geotiffwrite(amsr2BtTifDailyPath,btArray,amsr2Ref,CoordRefSysCode=4326);
                        catch
                            disp(['有问题数据：', amsr2BtTifDailyName]);
                        end
                    end
                end
            end
        end
    end
end

%% 将AMSR2 BT数据从TIF格式转为Mat格式.
% 获取所有年份AMSR2 TIF文件的路径列表.
disp('获取所有年份AMSR2 TIF文件的路径列表.');
amsr2DailyPathAllList = cell(dateAllListN * cpListN, 1); % 每日7个波段, 2各极化.
amsr2TifPath = fullfile(rootPath, 'AMSR2_2_CN_TIF');
amsr2ChannelList = dir(fullfile(amsr2TifPath, 'L3*'));
amsr2ChannelList = {amsr2ChannelList.name}';
startIndex = 1;
for i = 1 : length(amsr2ChannelList)
    amsr2ChannelPath = fullfile(amsr2TifPath, amsr2ChannelList{i});
    amsr2YearList = dir(fullfile(amsr2ChannelPath));
    amsr2YearList = {amsr2YearList(3:end).name}'; % 排除 '.', '..'两个文件夹。
    for j = 1 : length(amsr2YearList)
        amsr2YearPath = fullfile(amsr2ChannelPath, amsr2YearList{j});
        amsr2MonthList = dir(amsr2YearPath);
        amsr2MonthList = {amsr2MonthList(3:end).name}'; % 排除 '.', '..'两个文件夹。
        for k = 1 : length(amsr2MonthList)
            amsr2MonthPath = fullfile(amsr2YearPath, amsr2MonthList{k});
            amsr2DailyList = dir(fullfile(amsr2MonthPath, sprintf('*01D_EQM%s*.tif', orbit)));
            amsr2DailyList = {amsr2DailyList.name}';
            amsr2DailyListN = length(amsr2DailyList);
            amsr2DailyPath = fullfile(amsr2MonthPath, amsr2DailyList);
            endIndex = startIndex + amsr2DailyListN - 1;
            amsr2DailyPathAllList(startIndex : endIndex) = amsr2DailyPath;
            startIndex = endIndex + 1;
        end
    end
end
amsr2DailyPathAllList(cellfun(@isempty, amsr2DailyPathAllList)) = [];

% 读取AMSR2 BT数据, 并存储为Mat格式.
disp('按年份, 通道, 极化和昼夜读取AMSR2亮温数据, 并存储为mat格式.');
amsr2RowN = amsr2Ref.RasterSize(1); amsr2ColN = amsr2Ref.RasterSize(2);
orbitIndex = contains(amsr2DailyPathAllList, ['EQM', orbit]);  % 升降轨.
for i = 1 : yearNumListN
    yearStr = num2str(yearNumList(i));

    % 检查是否已存在当年的AMSR2 BT mat格式文件。
    amsr2YearArrayMatName = sprintf('AMSR2_BT_%s_%s.mat', yearStr, daynight);
    amsr2YearArrayMatPath = fullfile(amsr2MatPath, amsr2YearArrayMatName);
    if exist(amsr2YearArrayMatPath, 'file')
        continue
    end

    % 索引当年的AMSR2亮温数据存储路径, 并分极化方式与波段存储为mat格式文件.
    yearIndex1 = contains(amsr2DailyPathAllList, [yearStr,'\']);
    yearIndex2 = strcmp(extractBefore(dateAllList, 5), yearStr);
    % 极化排序：[H, V]，通道排序：[10, 18, 23, 36, 6, 7, 89]
    amsr2DailyPathYearList = amsr2DailyPathAllList(yearIndex1);
    dateYearList = dateAllList(yearIndex2);
    dateYearListN = length(dateYearList);
    [amsr2H06YearArray, amsr2V06YearArray, amsr2H07YearArray, amsr2V07YearArray, ...
        amsr2H10YearArray, amsr2V10YearArray, amsr2H18YearArray, amsr2V18YearArray, ...
        amsr2H23YearArray, amsr2V23YearArray, amsr2H36YearArray, amsr2V36YearArray, ...
        amsr2H89YearArray, amsr2V89YearArray] = ...
        deal(zeros(amsr2RowN, amsr2ColN, dateYearListN, 'uint16'));
    amsr2PathMatrix = cell(dateYearListN + 1, cpListN);
    amsr2PathMatrix(1, :) = cpList'; % 第一行存极化通道名称, 其他行存各日期文件路径. 
    emptyRowIndex = false(dateYearListN + 1, 1);
    for j = 1 : dateYearListN
        dateYear = dateYearList{j};
        dateYearIndex = contains(amsr2DailyPathYearList, dateYear);
        amsr2DailyPathList = amsr2DailyPathYearList(dateYearIndex);
        amsr2DailyPathListN = length(amsr2DailyPathList);
        if amsr2DailyPathListN ~= cpListN  % cpListN == 14
            emptyRowIndex(j+1) = true;
            fprintf('  %s %s的AMSR2亮温数据有缺失, 请检查.\n', dateYear, daynight);
            continue;
        end

        fprintf('  读取%s %s的AMSR2亮温数据.\n', dateYear, daynight);
        amsr2DailyArray = zeros(amsr2RowN, amsr2ColN, amsr2DailyPathListN, 'uint16');
        for k = 1 : amsr2DailyPathListN
            amsr2DailyArray(:, :, k) = readgeoraster(amsr2DailyPathList{k});
            amsr2PathMatrix(j+1, k) = amsr2DailyPathList(k);
        end
        amsr2H10YearArray(:, :, j) = amsr2DailyArray(:, :, 1);
        amsr2V10YearArray(:, :, j) = amsr2DailyArray(:, :, 2);
        amsr2H18YearArray(:, :, j) = amsr2DailyArray(:, :, 3);
        amsr2V18YearArray(:, :, j) = amsr2DailyArray(:, :, 4);
        amsr2H23YearArray(:, :, j) = amsr2DailyArray(:, :, 5);
        amsr2V23YearArray(:, :, j) = amsr2DailyArray(:, :, 6);
        amsr2H36YearArray(:, :, j) = amsr2DailyArray(:, :, 7);
        amsr2V36YearArray(:, :, j) = amsr2DailyArray(:, :, 8);
        amsr2H06YearArray(:, :, j) = amsr2DailyArray(:, :, 9);
        amsr2V06YearArray(:, :, j) = amsr2DailyArray(:, :, 10);
        amsr2H07YearArray(:, :, j) = amsr2DailyArray(:, :, 11);
        amsr2V07YearArray(:, :, j) = amsr2DailyArray(:, :, 12);
        amsr2H89YearArray(:, :, j) = amsr2DailyArray(:, :, 13);
        amsr2V89YearArray(:, :, j) = amsr2DailyArray(:, :, 14);
    end

    % 删除数据缺失日期行.
    amsr2PathMatrix(emptyRowIndex, :) = []; emptyRowIndex = emptyRowIndex(2:end);
    amsr2H10YearArray(:, :, emptyRowIndex) = []; amsr2V10YearArray(:, :, emptyRowIndex) = [];
    amsr2H18YearArray(:, :, emptyRowIndex) = []; amsr2V18YearArray(:, :, emptyRowIndex) = [];
    amsr2H23YearArray(:, :, emptyRowIndex) = []; amsr2V23YearArray(:, :, emptyRowIndex) = [];
    amsr2H36YearArray(:, :, emptyRowIndex) = []; amsr2V36YearArray(:, :, emptyRowIndex) = [];
    amsr2H06YearArray(:, :, emptyRowIndex) = []; amsr2V06YearArray(:, :, emptyRowIndex) = [];
    amsr2H07YearArray(:, :, emptyRowIndex) = []; amsr2V07YearArray(:, :, emptyRowIndex) = [];
    amsr2H89YearArray(:, :, emptyRowIndex) = []; amsr2V89YearArray(:, :, emptyRowIndex) = [];

    fprintf('  保存%s年%s的AMSR2亮温数据。\n', yearStr, daynight);
    save(amsr2YearArrayMatPath, 'amsr2Ref', 'amsr2*YearArray', 'amsr2PathMatrix');
end
