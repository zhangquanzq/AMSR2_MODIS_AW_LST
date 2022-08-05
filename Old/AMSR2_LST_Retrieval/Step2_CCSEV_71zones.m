%% AMSR2 LST反演中, 环境变量综合分类体系的构建.

%% 预设参数
% 有数据的年份月份列表(时间区间：2012/07/02-2019/12/31).
yearList = 2012 : 2019;
yearListN = length(yearList);

% 重采样时, 纯像元阈值.
lcThreshold = 0.6;
snowThreshold = 0.6;

% LC和分区影像代码的重分类. 主要针对水域, 冰川, 建筑这三个LC区中LST的反演.
%   [阶梯1_东北的东部: 1: 1-4, 阶梯1_华北: 2: 5-8, 阶梯1_华南: 3: 9-13]
%   [阶梯2_西南: 4: 14-17, 阶梯2_西北的东部: 5: 18-27, 阶梯2_东北的西部: 6: 28-30]
%   [阶梯2_西北的西部： 7： 31-53， 阶梯3_青藏高原: 8: 54-71]
% 重分类时, 使用了右开区间, 导致最后一个分区编码71排除在分类结果中, 71+1, 使其包含在重分类结果中.
regionNodes = [1, 5, 9, 14, 18, 28, 31, 54, 71+1]; 
regionsN = length(regionNodes) - 1;

%% 路径.
% 根目录
rootPath = 'G:\AMSR_MODIS_Fusion';
dataPath = fullfile(rootPath, 'Data');
% funcPath = fullfile(rootPath, 'Code\Functions');
funcPath = fullfile('F:\AMSR_MODIS_Fusion\Code\Functions');
addpath(funcPath);

% 创建保存环境变量综合分类体系数据的文件夹.
ccsevFolder = 'CCSEV_Matlab';
ccsevPath = fullfile(dataPath, ccsevFolder);
if ~exist(ccsevPath, 'dir')
    mkdir(ccsevPath)
end

% 地表覆盖和分区数据路径. 其中,zones和desert长期不变, MODIS地表覆盖每年一幅, MODIS积雪每日一幅.
zonesPath = fullfile(dataPath, 'Zones', 'GeographicalZones_merge.tif');
desertPath = fullfile(dataPath, 'Landcover', 'Desert_CN_gcs.tif');
modisLcPath = fullfile(dataPath, 'MCD12Q1_2_CN_TIF', 'MCD12Q1.A2012001.006_gcs_reclsfy_resmpl.tif');
modisSnowPath = fullfile(dataPath, 'MYD10C1_2_CN_TIF/MYD10C1_2012XXX', ...
    'MYD10C1.A2012002.006.2016081154447.tif');

% 事例AMSR2 BT数据路径.
amsr2BtPath = fullfile(dataPath, 'AMSR2_3_MaskCn_TIF/AMSR2_2012XXXX_TIF/AMSR2_20120703', ...
    'GW1AM2_20120703_01D_EQMA_L3SGT06HA2220220_BtH.tif');

%% 读取数据参考与属性.
% 读取掩膜后的AMSR2 BT图像和MODIS LST图像的空间参考. 所有掩膜AMSR2 BT影像与MODIS LST影像的空间参考都相同.
amsr2Ref = georasterinfo(amsr2BtPath).RasterReference;
amsr2RowN = amsr2Ref.RasterSize(1); amsr2ColN = amsr2Ref.RasterSize(2);

% 读取MODIS LC图像的空间参考. 沙漠图像与土地覆盖图像的空间基准相同.
modisLcRef = georasterinfo(modisLcPath).RasterReference;
lcRowN = modisLcRef.RasterSize(1); lcColN = modisLcRef.RasterSize(2);

% 读取MODIS Snow图像的空间参考.
modisSnowRef = georasterinfo(modisSnowPath).RasterReference;
snowRowN = modisSnowRef.RasterSize(1); snowColN = modisSnowRef.RasterSize(2);

% 获取LC影像左上角第一个滑动窗口的边界行列号, 行列数.
% 滑动窗口的大小为AMSR2影像的像元尺寸. 滑动窗口的行列数是LC影像和AMSR2影像的分辨率倍数. 
% 边界行列号: [topRow, bottomRow, leftCol, rightCol],  行列数: [blockRowN, blockColN].
[lcBlockBdy, lcBlockSize] = getStartBlockRowCol(modisLcRef, amsr2Ref);
lcBlockN = prod(lcBlockSize);

[snowBlockBdy, snowBlockSize] = getStartBlockRowCol(modisSnowRef, amsr2Ref);
snowBlockN = prod(snowBlockSize);

%% 合并地表覆盖、沙漠、分区和积雪分区影像。
% 获取沙漠影像和编码.
desertLayer = single(readgeoraster(desertPath));
desertCodeList = unique(desertLayer); % [1000, 1100, ..., 2000, 65535], 65535是Nodata.
desertNodataCode = desertCodeList(end);
desertCodeList = desertCodeList(1:end-1);
desertCodeListN = length(desertCodeList);

% 获取地形分区影像.
zonesLayer = single(readgeoraster(zonesPath));
zonesLayer(zonesLayer == 128) = nan;

% 获取LC编码信息. MODIS LC编码: [裸地: 1000, 植被: 2000, 水体: 3000, 冰川: 4000, 建筑: 5000]
lcLayer = single(readgeoraster(modisLcPath)) * 1000;
lcCodeList = unique(lcLayer);  % [1, 2, 3, 4, 5, 15] * 1000, 15是Nodata.
lcCodeList = lcCodeList(1:end-1);  % 去掉Nodata.
removedLcCodes = lcCodeList(1:2); % 裸地, 植被.
keptLcCodes = lcCodeList(3:5);  % 水体, 冰川, 建筑.
lcCodeList = [0; keptLcCodes; 6000; 10000]; % [其他:0, 积雪: 6000, 混合: 10000]
lcCodeListN = length(lcCodeList);

% 分年度合并不同LC, 并升尺度到AMSR2数据的分辨率.
for i = 1 : yearListN
    yearStr = num2str(yearList(i));
    fprintf('创建环境变量综合分类体系(CCSEV)年份: %s\n', yearStr)

    % 保存环境变量综合分类体系数据的路径.
    outLcPctArrayMatPath = fullfile(ccsevPath, ['CCSEV_', yearStr, '.mat']);
    if exist(outLcPctArrayMatPath, 'file')
        continue
    end

    % 获取地表覆盖类型影像. 移除原数据中的裸地和植被类型, 然后将剩余的地表覆盖类型与沙漠叠置. 当水体, 冰雪,
    %   建筑与沙漠有重叠时, 保持他们不变.
    modisLcName = sprintf('MCD12Q1.A%s001.006_gcs_reclsfy_resmpl.tif', yearStr);
    modisLcPath = fullfile(dataPath, 'MCD12Q1_2_CN_TIF', modisLcName);
    lcLayer = single(readgeoraster(modisLcPath)) * 1000;
    lcLayer(ismember(lcLayer, removedLcCodes)) = 0;
    desertIndex = (desertLayer ~= desertNodataCode) & ~ismember(lcLayer, keptLcCodes);
    lcLayer(desertIndex) = desertLayer(desertIndex);

    % 叠置地表覆盖影像和分区影像. 除了积雪之外的其他地表覆盖类型与分区影像在此处结合. 因为积雪覆盖的面积每天
    %   都在变化, 所以积雪覆盖在下一步处理. 混合像元定义为像元内没有任何地表覆盖类型的面积超过60%的像元. 合
    %   并的地表覆盖编码定义为: 分区代码 + 地表覆盖代码. 例如: 3010表示在10分区的水域(3000).
    [othersPctLayer, desertPctLayer, waterPctLayer, glacierPctLayer, ...
        buildingPctLayer, desertCodeLayer] = deal(zeros(amsr2RowN, amsr2ColN) * nan);
    zonesLcLayer = zonesLayer;
    for ii = 1 : amsr2RowN
        for jj = 1 : amsr2ColN
            % 跳过Zonal图层的Nodata像元, 既非中国区域.
            if isnan(zonesLcLayer(ii, jj))
                continue
            end
            % 定位每一个LC滑动窗口.
            lcBlockTopRow = lcBlockBdy(1) + lcBlockSize(1) * (ii - 1);
            lcBlockBottomRow = lcBlockBdy(2) + lcBlockSize(1) * (ii - 1);
            lcBlockLeftCol = lcBlockBdy(3) + lcBlockSize(2) * (jj - 1);
            lcBlockRightCol = lcBlockBdy(4) + lcBlockSize(2) * (jj - 1);
            if lcBlockRightCol > lcColN || lcBlockBottomRow > lcRowN
                continue
            end
            % 窗口计算.
            lcBlock = lcLayer(lcBlockTopRow : lcBlockBottomRow, lcBlockLeftCol : lcBlockRightCol);
            othersPctLayer(ii, jj) = sum(lcBlock(:) == lcCodeList(1)) / lcBlockN;
            desertPctLayer(ii, jj) = sum(ismember(lcBlock(:), desertCodeList)) / lcBlockN;
            waterPctLayer(ii, jj) = sum(lcBlock(:) == lcCodeList(2)) / lcBlockN;
            glacierPctLayer(ii, jj) = sum(lcBlock(:) == lcCodeList(3)) / lcBlockN;
            buildingPctLayer(ii, jj) = sum(lcBlock(:) == lcCodeList(4)) / lcBlockN;
            lcPercents = [othersPctLayer(ii, jj), desertPctLayer(ii, jj), waterPctLayer(ii, jj), ...
                glacierPctLayer(ii, jj), buildingPctLayer(ii, jj)]';

            % 如果单个窗口(1个AMSR2像元)内的某土地覆盖类型面积超过窗口面积的60%，则将其视为纯像元, 否则视为
            %   混合像元.
            maxLcPct = max(lcPercents);
            maxIndex = find(lcPercents == maxLcPct);
            if length(maxIndex) == 1 % 最大值数量只有1个时, 单个土地覆盖类型面积比例有超过0.6的可能.
                if maxLcPct >= lcThreshold % 纯像元.
                    % maxIndex == 1 为其他地表覆盖类型(植被, 裸地), 不做处理.
                    if maxIndex == 2 % 沙漠.
                        zonesLcLayer(ii, jj) = mode(lcBlock(:));  % 众数.
                    elseif maxIndex >= 3 % 水域, 冰川, 建筑.
                        zonesLcLayer(ii, jj) = zonesLcLayer(ii, jj) + mode(lcBlock(:));
                    end
                else  % 混合像元.
                    zonesLcLayer(ii, jj) = zonesLcLayer(ii, jj) + lcCodeList(6);
                end
            else % 最大值数量超过1个时, 肯定是混合像元.
                zonesLcLayer(ii, jj) = zonesLcLayer(ii, jj) + lcCodeList(6);
            end

            % 获取与AMSR2影像分辨率相同的沙漠编码数据层.
            if desertPctLayer(ii, jj) > 0
                desertPixelCount = zeros(desertCodeListN, 1);
                for n = 1 : desertCodeListN
                    desertPixelCount(n) = sum(lcBlock(:) == desertCodeList(n));
                end
                maxDesertPixelCount = max(desertPixelCount);
                if sum(desertPixelCount) > lcBlockN / 2
                    majorDesertIndex = find(desertPixelCount == maxDesertPixelCount);
                    majorDesertIndex = majorDesertIndex(1);
                    desertCodeLayer(ii, jj) = desertCodeList(majorDesertIndex);
                end
            end
        end
    end

    % 读取每天的MODIS积雪图像, 将其混合到分区地表覆盖图像中，并修正土地覆盖百分比.
    modisSnowYearFolder = ['MYD10C1_', yearStr, 'XXX'];
    modisSnowYearPath = fullfile(dataPath, 'MYD10C1_2_CN_TIF', modisSnowYearFolder);
    modisSnowDailyFileList = dir(fullfile(modisSnowYearPath, 'MYD10C1*.tif'));
    modisSnowDailyFileList = {modisSnowDailyFileList.name}';
    modisSnowDailyFileListN = length(modisSnowDailyFileList);
    dateList = cell(modisSnowDailyFileListN, 1);
    
    [zonesLcArray, otherPctArray, desertPctArray, glacierPctArray, waterPctArray, snowPctArray, ...
        buildingPctArray] = deal(zeros(amsr2RowN, amsr2ColN, modisSnowDailyFileListN));
    % snowZonesLayer表示积雪所在的分区.
    snowZonesLayer = zonesLayer + lcCodeList(5);
    for j = 1 : modisSnowDailyFileListN
        modisSnowDailyFile = modisSnowDailyFileList{j};
        yearDateStr = split(modisSnowDailyFile, '.');
        dateList{j} = yearDateStr{2}(2:end);

        fprintf('修正积雪覆盖比例: %s\n', modisSnowDailyFile);

        % 处理并行计算中问题的技术.
        tempZonesLcLayer = zonesLcLayer;
        tempSnowZonesLayer = snowZonesLayer;

        % 获取有MODIS积雪覆盖数据的日期.
        modisSnowLayer = readgeoraster(fullfile(modisSnowYearPath, modisSnowDailyFile));
        modisSnowLayer = single(modisSnowLayer);
        % 获取积雪比例数据层.
        snowPctLayer = zeros(amsr2RowN, amsr2ColN);
        for ii = 1 : amsr2RowN
            for jj = 1 : amsr2ColN
                % 跳过zonesLcLayer的Nodata像元, 既非中国区域.
                if isnan(zonesLcLayer(ii, jj))
                    continue
                end

                % 定位每一个snow滑动窗口.
                snowBlockTopRow = snowBlockBdy(1) + snowBlockSize(1) * (ii - 1);
                snowBlockBottomRow = snowBlockBdy(2) + snowBlockSize(1) * (ii - 1);
                snowBlockLeftCol = snowBlockBdy(3) + snowBlockSize(2) * (jj - 1);
                snowBlockRightCol = snowBlockBdy(4) + snowBlockSize(2) * (jj - 1);
                if snowBlockRightCol > snowColN || snowBlockBottomRow > snowRowN
                    continue
                end

                % 积雪滑动窗口计算.
                snowBlock = modisSnowLayer(snowBlockTopRow:snowBlockBottomRow, ...
                    snowBlockLeftCol:snowBlockRightCol);
                snowPixelIndex = (snowBlock >= 0) & (snowBlock <= 100);
                snowPixelN = sum(snowPixelIndex(:));
                snowPctLayer(ii, jj) = sum(snowBlock(snowPixelIndex)) ./ (100 * snowPixelN);

%                 % 如果原沙漠像元有一半以上面积由积雪覆盖, 则将此像元视为非沙漠像元.
%                 if snowPctLayer(ii, jj) > 0.5
%                     desertCodeLayer(ii, jj) = nan;
%                 end
            end
        end
        % 假设在Nodata像元上没有积雪覆盖.
        snowPctLayer(isnan(snowPctLayer)) = 0; 
        snowPctArray(:, :, j) = snowPctLayer;

        % 将积雪覆盖融合到中国区域的Zonal LC图像中.
        pureSnowIndexLayer = snowPctLayer > snowThreshold;
        tempZonesLcLayer(pureSnowIndexLayer) = tempSnowZonesLayer(pureSnowIndexLayer);
        zonesLcArray(:, :, j) = tempZonesLcLayer;

        % 修正地表覆盖比例数据层.
        otherPctArray(:, :, j) = othersPctLayer .* (1 - snowPctLayer);
        desertPctArray(:, :, j) = desertPctLayer .* (1 - snowPctLayer);
        waterPctArray(:, :, j) = waterPctLayer .* (1 - snowPctLayer);
        glacierPctArray(:, :, j) = glacierPctLayer .* (1 - snowPctLayer);
        buildingPctArray(:, :, j) = buildingPctLayer .* (1 - snowPctLayer);
    end

    % LC和分区影像代码的重分类.
    %   [阶梯1_东北的东部: 1: 1-4, 阶梯1_华北: 2: 5-8, 阶梯1_华南: 3: 9-13]
    %   [阶梯2_西南: 4: 14-17, 阶梯2_西北的东部: 5: 18-27, 阶梯2_东北的西部: 6: 28-30]
    %   [阶梯2_西北的西部： 7： 31-53， 阶梯3_青藏高原: 8: 54-71]
    fprintf('中国8大地形分区: %s\n', yearStr);
    zonesLcIdList = unique(zonesLcArray);
    zonesLcIdList(isnan(zonesLcIdList)) = [];
    fixedZonesLcArray = zonesLcArray;
    for j = 1 : regionsN
        for k = 2 : lcCodeListN - 1
            lcIndexVector = (zonesLcIdList >= lcCodeList(k) + regionNodes(j)) & ...
                (zonesLcIdList < lcCodeList(k) + regionNodes(j+1));
            zonesLcIdList(lcIndexVector) = lcCodeList(k) + j;
            lcIndexArray = (zonesLcArray >= lcCodeList(k) + regionNodes(j)) & ...
                (zonesLcArray < lcCodeList(k) + regionNodes(j+1));
            fixedZonesLcArray(lcIndexArray) = lcCodeList(k) + j;
        end
    end
    zonesLcIdList = unique(zonesLcIdList);
    zonesLcIdListN = numel(zonesLcIdList);
    mixedLcN = sum(zonesLcIdList >= lcCodeList(end));
    pureLcIdN = zonesLcIdListN - mixedLcN;

    % 变量 zonesLcArray 可以不保存. 
    % desertCodeLayer 受到积雪覆盖变化影像, 每天也在变化, 可能需要修正后才能用.
    fprintf('保存%s年的CCESV.\n', yearStr)
    save(outLcPctArrayMatPath, 'otherPctArray', 'desertPctArray', 'waterPctArray', ...
        'glacierPctArray', 'buildingPctArray', 'snowPctArray', 'zonesLcArray', ...
        'fixedZonesLcArray', 'desertCodeLayer', 'zonesLcIdList', 'lcCodeList', 'dateList');
end
