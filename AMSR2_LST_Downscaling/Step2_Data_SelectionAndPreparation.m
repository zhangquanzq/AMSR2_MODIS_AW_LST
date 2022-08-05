%% 用于AMSR2 LST降尺度试验的数据筛选和准备.

%% 标记和预设参数.
% 指定区域的标记. 1表示内蒙古, 2表示晋南豫西, 3表示云贵高原, 4表示北大河, 5表示那曲.
flg1 = 2;
% 指定昼夜的标记. 1表示白天, 2表示夜晚.
flg2 = 1;
% 指定降尺度策略的标记. 1表示逐级降尺度, 2表示直接降尺度.
flg3 = 1;

region = {'NMG', 'JNYX', 'YGGY', 'BDH', 'Naqu'}; region = region{flg1};
dayNight = {'Day', 'Night'}; dayNight = dayNight{flg2};
cellsizeList = {{'0.01', '0.03', '0.1'}, {'0.01', '0.1'}}; cellsizeList = cellsizeList{flg3};
cellsizeListN = length(cellsizeList);
topoParams = {'Elev', 'Slp'};

yearNum = 2012;
cloudPctEdges = [0, 0.1];  % 针对降尺度试验.
% cloudPctEdges = [0.1, 0.9];  % 针对融合试验.
modisLstQcCode = [0, 17];  % 只保留0和17的话, 满足云比例条件的数据量大幅减小.

dayNightRegion = sprintf('%s_%s', dayNight, region);

%% 路径.
% 根目录.
rootDir1 = 'J:\AMSR2_LST_Retrieval';
rootDir2 = 'J:\AMSR2_LST_Downscaling';
dataDir1 = fullfile(rootDir1, 'Data');
dataDir2 = fullfile(rootDir2, 'Data');
addpath(fullfile(rootDir2, 'Code\functions'));

% 不同类型数据总目录.
amsr2LstDir = fullfile(dataDir1, 'AMSR2_4_LstCn_TIF');
modisLstDir = fullfile(dataDir1, 'MYD11A1_2_PrjCN_TIF');
ndviDir = fullfile(dataDir2, 'MYD13A2_2_PrjCN_TIF');
regionExtentPath = fullfile(dataDir2, 'Feature', sprintf('%s_Extent.shp', region));
regionBufferPath = fullfile(dataDir2, 'Feature', sprintf('%s_Buffer.shp', region));

% 研究区和年份文件夹路径.
regionDir = fullfile(dataDir2, sprintf('Region_%s', region));
if ~exist(regionDir, 'dir')
    mkdir(regionDir)
end
yearDir = fullfile(regionDir, sprintf('%s_%d_%s', region, yearNum, dayNight));
if ~exist(yearDir, 'dir')
    mkdir(yearDir)
end

% 筛选数据指标统计文件路径.
staRegionCsvName = sprintf('DateSelect_%s_%d_%s.csv', region, yearNum, dayNight);
staRegionCsvPath = fullfile(regionDir, staRegionCsvName);

% 读取MODIS地表温度数据文件名和日期列表.
modisLstYearDir = fullfile(modisLstDir, sprintf('MYD11A1_%dXXX_TIF', yearNum));

modisLstNameList = {dir(fullfile(modisLstYearDir, sprintf('MYD11A1*LST_%s.tif', dayNight))).name}';
modisLstNameListN = length(modisLstNameList);
modisDateList = string(cellfun(@(x) yday2ymd(x(10:16)), modisLstNameList, UniformOutput=false));

modisQcNameList = {dir(fullfile(modisLstYearDir, sprintf('MYD11A1*QC_%s.tif', dayNight))).name}';
modisQcNameListN = length(modisQcNameList);
modisQcDateList = string(cellfun(@(x) yday2ymd(x(10:16)), modisQcNameList, UniformOutput=false));

if (modisLstNameListN ~= modisQcNameListN) || sum(modisDateList ~= modisQcDateList) > 0
    error('MODIS LST数据与QC数据不匹配, 请检查.')
end
modisLstPathList = fullfile(modisLstYearDir, modisLstNameList);
modisLstQcPathList = fullfile(modisLstYearDir, modisQcNameList);

% 读取AMSRE/2地表温度数据文件名和日期列表.
amsr2LstYearDir = fullfile(amsr2LstDir, sprintf('AMSR2_LST_%dXXXX_TIF', yearNum));
amsr2LstNameList = {dir(fullfile(amsr2LstYearDir, sprintf('AMSR2*%s*.TIF', dayNight))).name}';
amsr2DateList = string(cellfun(@(x) x(15:22), amsr2LstNameList, UniformOutput=false));
amsr2LstPathList = fullfile(amsr2LstYearDir, amsr2LstNameList);

% 读取MODIS NDVI数据文件名和日期列表.
ndviYearDir = fullfile(ndviDir, sprintf('MYD13A2_%dXXX_TIF', yearNum));
ndviPastYearDir = fullfile(ndviDir, sprintf('MYD13A2_%dXXX_TIF', yearNum - 1));
ndviNextYearDir = fullfile(ndviDir, sprintf('MYD13A2_%dXXX_TIF', yearNum + 1));

ndvi361Path = fullfile(ndviPastYearDir, sprintf('MYD13A2.A%d361.061_NDVI.tif', yearNum - 1));
ndvi009Path = fullfile(ndviNextYearDir, sprintf('MYD13A2.A%d009.061_NDVI.tif', yearNum + 1));
ndviQa361Path = fullfile(ndviPastYearDir, sprintf('MYD13A2.A%d361.061_QA.tif', yearNum - 1));
ndviQa009Path = fullfile(ndviNextYearDir, sprintf('MYD13A2.A%d009.061_QA.tif', yearNum + 1));
if ~exist(ndvi361Path, 'file') || ~exist(ndvi009Path, 'file') || ...
        ~exist(ndviQa361Path, 'file') || ~exist(ndviQa009Path, 'file')
    error('上一年最后一个或下一年第一个MODIS NDVI数据与QA数据不存在, 请检查.')
end

ndviNameList = {dir(fullfile(ndviYearDir, 'MYD13A2*NDVI.tif')).name}';
ndviPathList = fullfile(ndviYearDir, ndviNameList);
ndviPathList = [ndvi361Path; ndviPathList; ndvi009Path];

ndviQaNameList = {dir(fullfile(ndviYearDir, 'MYD13A2*QA.tif')).name}';
ndviQaPathList = fullfile(ndviYearDir, ndviQaNameList);
ndviQaPathList = [ndviQa361Path; ndviQaPathList; ndviQa009Path];

[~, ndviNameList, ~] = cellfun(@(x) fileparts(x), ndviPathList, UniformOutput=false);
ndviDateList = string(cellfun(@(x) yday2ymd(x(10:16)), ndviNameList, UniformOutput=false));

[~, ndviQaNameList, ~] = cellfun(@(x) fileparts(x), ndviQaPathList, UniformOutput=false);
ndviQaDateList = string(cellfun(@(x) yday2ymd(x(10:16)), ndviQaNameList, UniformOutput=false));

if (length(ndviPathList) ~= length(ndviQaPathList)) || sum(ndviDateList ~= ndviQaDateList) > 0
    error('MODIS NDVI数据与QA数据不匹配, 请检查.')
end

%% SRTM数据准备. 裁剪研究区的SRTM高程, 坡度数据, 并生成分辨率序列.
srtmRegionDir = fullfile(regionDir, sprintf('SRTM_%s', region));
if ~exist(srtmRegionDir, 'dir')
    mkdir(srtmRegionDir)
end

snapRasterPath = fullfile(dataDir1, 'AMSR2_2_CN_TIF\L3.TB6GHz_10\2012\07',...
    'GW1AM2_20120703_01D_EQMA_L3SGT06HA2220220_BtH.tif');
for i = 1: length(topoParams)
    topo = topoParams{i};

    % 裁剪研究区SRTM数据.
    srtmTopoPath = fullfile(dataDir1, 'SRTM', sprintf('SRTM_0d01_CN_%s.tif', topo));
    srtmTopoExtentPath = fullfile(srtmRegionDir, sprintf('SRTM_0.01_%s_%s.tif', region, topo));
    if ~exist(srtmTopoExtentPath, 'file')
        clipRaster(srtmTopoPath, regionExtentPath, srtmTopoExtentPath);
    end
    srtmTopoBufferPath = fullfile(srtmRegionDir, sprintf('SRTM_0.01_buffer_%s.tif', topo));
    if ~exist(srtmTopoBufferPath, 'file')
        clipRaster(srtmTopoPath, regionBufferPath, srtmTopoBufferPath);
    end

    % 重采样分辨率序列的SRTM数据.
    for j = 1: cellsizeListN
        cellsize = cellsizeList{j};

        srtmTopoResample1Name = sprintf('SRTM_%s_%s_%s.tif', cellsize, region, topo);
        srtmTopoResample1Path = fullfile(srtmRegionDir, srtmTopoResample1Name);
        if exist(srtmTopoResample1Path, 'file')
            continue
        end

        srtmTopoResample2Name = sprintf('SRTM_%s_buffer_%s.tif', cellsize, topo);
        srtmTopoResample2Path = fullfile(srtmRegionDir, srtmTopoResample2Name);
        if exist(srtmTopoResample2Path, 'file')
            delete(srtmTopoResample2Path)
        end

        resampleRaster(srtmTopoBufferPath, srtmTopoResample2Path, cellsize, snapRasterPath);
        clipRaster(srtmTopoResample2Path, regionExtentPath, srtmTopoResample1Path);
        delete(srtmTopoResample2Path)
    end
    delete(srtmTopoBufferPath)
end

%% 筛选可用日期数据.
% 创建年度日期列表.
yearDateRange = datetime(sprintf('%d-01-01', yearNum)): datetime(sprintf('%d-12-31', yearNum));
yearDateList = string(yearDateRange', 'yyyyMMdd');

% 数据日期筛选.
writelines("YearDate,PixelPercent,R", staRegionCsvPath);
sprintf('筛选%s区%d年%s的数据.\n', region, yearNum, dayNight);
for i = 1: length(yearDateList)
    yearDate = yearDateList(i);

    % 获取均有数据的AMSR2和MODIS地表温度文件路径, 跳过至少有其中一个没有数据的日期.
    [dateIndex1, dateLocate1] = ismember(yearDate, amsr2DateList);
    [dateIndex2, dateLocate2] = ismember(yearDate, modisDateList);
    if dateIndex1 == 1 && dateIndex2 == 1
        amsr2LstPath = amsr2LstPathList{dateLocate1};
        modisLstPath = modisLstPathList{dateLocate2};
        modisQcPath = modisLstQcPathList{dateLocate2};
    else
        continue
    end

    % 创建AMSR2和MODIS都有数据的日期文件夹.
    yearDateDir = fullfile(yearDir, sprintf('%s_%s_%s', region, yearDate, dayNight));
    if ~exist(yearDateDir, 'file')
        mkdir(yearDateDir)
    end

    % 裁剪研究区范围的AMSR2数据.
    [~, amsr2Lst, ext] = fileparts(amsr2LstPath);
    amsr2LstRegionName = replace([amsr2Lst, ext], dayNight, dayNightRegion);
    amsr2LstRegionPath = fullfile(yearDateDir, amsr2LstRegionName);
    if ~exist(amsr2LstRegionPath, 'file')
        clipRaster(amsr2LstPath, regionExtentPath, amsr2LstRegionPath);
    end

    % 如果AMSR2地表温度影像中存在Nodata或0, 则舍弃当日的数据.
    amsr2LstLayer = readgeoraster(amsr2LstRegionPath);
    amsr2LstNodata = georasterinfo(amsr2LstRegionPath).MissingDataIndicator;
    if ~isempty(amsr2LstNodata)
        amsr2LstLayer(amsr2LstLayer == amsr2LstNodata) = 0;
    end
    if sum(amsr2LstLayer == 0, 'all') > 0
        rmdir(yearDateDir, 's')
        continue
    end

    % 裁剪研究区MODIS数据, 计算云覆盖面积, 确定数据可用性. 若不可用, 删除当天AMSR2和MODIS数据.
    [~, modisLst, ext] = fileparts(modisLstPath);
    modisLstRegionName = replace([modisLst, ext], dayNight, [dayNightRegion, '_0.01']);
    modisLstRegionPath = fullfile(yearDateDir, modisLstRegionName);
    if ~exist(modisLstRegionPath, 'file')
        clipRaster(modisLstPath, regionExtentPath, modisLstRegionPath);
    end

    [~, modisQc, ext] = fileparts(modisQcPath);
    modisQcName = [modisQc, ext];
    modisQcRegionName = replace(modisQcName, dayNight, [dayNightRegion, '_0.01']);
    modisQcRegionPath = fullfile(yearDateDir, modisQcRegionName);
    if ~exist(modisQcRegionPath, 'file')
        clipRaster(modisQcPath, regionExtentPath, modisQcRegionPath);
    end

    [modisLstLayer, modisLstRef] = readgeoraster(modisLstRegionPath);
    modisLstNodata = georasterinfo(modisLstRegionPath).MissingDataIndicator;
    if ~isempty(modisLstNodata)
        modisLstLayer(modisLstLayer == modisLstNodata) = 0;
    end
    modisLstLayer = single(modisLstLayer) * 0.02;
    modisLstLayer(~ismember(readgeoraster(modisQcRegionPath), modisLstQcCode)) = 0;
    modisNodateLayer = (modisLstLayer == 0);
    modisNodataPct = sum(modisNodateLayer, 'all') / numel(modisNodateLayer);
    if modisNodataPct < cloudPctEdges(1) || modisNodataPct > cloudPctEdges(2)
        rmdir(yearDateDir, 's')
        continue
    end

    modisLstLayer(modisLstLayer == 0) = nan;
    geotiffwrite(modisLstRegionPath, modisLstLayer, modisLstRef, CoordRefSysCode=4326);
    delete(modisQcRegionPath);

    % 创建存放重采样过程数据的临时文件夹.
    tempDir = fullfile(dataDir2, sprintf('Temp_%s_%s_%s', region, yearDate, dayNight));
    if ~exist(tempDir, 'dir')
        mkdir(tempDir)
    end

    % 若MODIS影像的Nodata像元数在指定范围内, 重采样生成MODIS LST分辨率序列.
    modisLstRegionValidList = false(cellsizeListN, 1);
    for j = 1: cellsizeListN
        modisPath = replace(modisLstRegionPath, '0.01', cellsizeList(j));
        modisLstRegionValidList(j) = exist(modisPath, 'file');
    end
    if ismember(false, modisLstRegionValidList)
        modisLstBufferName = replace([modisLst, ext], dayNight, [dayNight, '_Buffer_0.01']);
        modisLstBufferPath = fullfile(tempDir, modisLstBufferName);
        clipRaster(modisLstPath, regionBufferPath, modisLstBufferPath); % 暂时没加*0.02和SetNull步骤.

        modisQcBufferName = replace([modisLst, ext], dayNight, [dayNight, '_Buffer_0.01']);
        modisQcBufferPath = fullfile(tempDir, modisQcBufferName);
        clipRaster(modisQcPath, regionBufferPath, modisQcBufferPath);

        % 保存质量控制后的MODIS LST.
        [modisLstLayer, modisLstRef] = readgeoraster(modisLstBufferPath);
        modisLstLayer = single(modisLstLayer) * 0.02;
        modisLstLayer(~ismember(readgeoraster(modisQcBufferPath), modisLstQcCode)) = nan;
        geotiffwrite(modisLstBufferPath, modisLstLayer, modisLstRef, CoordRefSysCode=4326);
        delete(modisQcBufferPath);
        
        % 重采样.
        resamplePixelSeries(modisLstBufferPath, regionExtentPath, yearDateDir, cellsizeList, ...
            region, 'initial', amsr2LstRegionPath)
        delete(modisLstBufferPath);
    end

    % 准备NDVI数据.
    % 获取离MODIS LST最近的NDVI数据日期.
    yearDate2 = datetime(yearDate, InputFormat='yyyyMMdd');
    ndviDateDiff = abs(yearDate2 - datetime(ndviDateList, InputFormat='yyyyMMdd'));

    % NDVI数据裁剪与质量控制. 裁剪使用研究区缓冲区范围, 避免重采样时边界像元插值问题.
    ndviPath = ndviPathList{ndviDateDiff == min(ndviDateDiff)};
    [~, ndvi, ndviExt] = fileparts(ndviPath);
    ndviRegionName = replace([ndvi, ndviExt], 'NDVI', sprintf('NDVI_%s_0.01', region));
    ndviRegionPath = fullfile(yearDateDir, ndviRegionName);

    ndviRegionValidList = false(cellsizeListN, 1);
    for j = 1: cellsizeListN
        nPath = replace(ndviRegionPath, '0.01', cellsizeList(j));
        ndviRegionValidList(j) = exist(nPath, 'file');
    end
    if ismember(false, ndviRegionValidList)
        ndviBufferName = replace([ndvi, ndviExt], 'NDVI', 'NDVI_Buffer_0.01');
        ndviBufferPath = fullfile(tempDir, ndviBufferName);
        if ~exist(ndviBufferPath, 'file')
            % 裁剪NDVI数据.
            ndviClipName = replace([ndvi, ndviExt], 'NDVI', 'NDVI_Clip_0.01');
            ndviClipPath = fullfile(tempDir, ndviClipName);
            clipRaster(ndviPath, regionBufferPath, ndviClipPath);

            % 裁剪QA数据.
            ndviQaPath = ndviQaPathList(ndviDateDiff == min(ndviDateDiff));
            [~, qaName, qaExt] = fileparts(ndviQaPath);
            ndviQaClipName = replace([qaName, qaExt], 'QA', 'QA_Clip_0.01');
            ndviQaClipPath = fullfile(tempDir, ndviQaClipName);
            clipRaster(ndviQaPath, regionBufferPath, ndviQaClipPath);

            % 对NDVI进行质量控制, 并导出控制后的NDVI数据.
            [ndviQaLayer, ndviQaRef] = readgeoraster(ndviQaClipPath);
            [qaRowN, qaColN] = size(ndviQaLayer);
            ndviQaBinVector = dec2bin(ndviQaLayer, 16);
            ndviQaBinVectorN = length(ndviQaBinVector);
            [ndviQaBinCutArray1, ndviQaBinCutArray2] = deal(strings(ndviQaBinVectorN, 1));
            for m = 1: ndviQaBinVectorN
                ndviQaBin = ndviQaBinVector(m, :);
                ndviQaBinCutArray1(m) = ndviQaBin(end - 1: end);
                ndviQaBinCutArray2(m) = ndviQaBin(end - 5: end - 4);
            end
            ndviQaBinCutArray1 = reshape(ndviQaBinCutArray1, qaRowN, qaColN);
            ndviQaBinCutArray2 = reshape(ndviQaBinCutArray2, qaRowN, qaColN);
            ndviIndexLayer = ((ndviQaBinCutArray1 == '00') | (ndviQaBinCutArray1 == '01') | ...
                (ndviQaBinCutArray2 ~= '11'));

            ndviClipLayer = readgeoraster(ndviClipPath) .* int16(ndviIndexLayer);
            geotiffwrite(ndviBufferPath, ndviClipLayer, ndviQaRef, CoordRefSysCode=4326);

            % 删除质量控制前的NDVI和QA Clip文件.
            delete(ndviClipPath);
            delete(ndviQaClipPath);
        end

        clipRaster(ndviBufferPath, regionExtentPath, ndviRegionPath);
        resamplePixelSeries(ndviBufferPath, regionExtentPath, yearDateDir, cellsizeList, ...
            region, 'initial', amsr2LstRegionPath)
        delete(ndviBufferPath)
    end
    rmdir(tempDir, 's')

    % 获取AMSR和MODIS地表温度数据的重叠区域, 并计算相关系数.
    modisLstLowestPath = replace(modisLstRegionPath, '0.01', cellsizeList(end));
    [modisLstLayer, amsr2LstLayer] = intersectRaster(modisLstLowestPath, amsr2LstRegionPath);
    
    validIndexLayer = ~isnan(modisLstLayer) & ~isnan(amsr2LstLayer);
    modisLstVector = modisLstLayer(validIndexLayer);
    amsr2LstVector = amsr2LstLayer(validIndexLayer);
    cc = corrcoef(modisLstVector, amsr2LstVector); cc = cc(2);
    fprintf("%s, %s, %s, Area: %.3f, R: %.3f.\n", region, yearDate, dayNight, 1-modisNodataPct, cc);

    % 将统计数据写入CSV文件.
    record = sprintf("%s,%.3f,%.3f", yearDate, 1 - modisNodataPct, cc);
    writelines(record, staRegionCsvPath, WriteMode="append");
end
