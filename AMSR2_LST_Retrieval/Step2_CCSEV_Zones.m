%% 构建AMSR2 LST反演中的环境变量综合分类体系.
% !!! 积雪覆盖目前也按8大综合分区反演地表温度, 但结果精度较低, 尤其东北新疆地区, 需重新分区. !!!

%% 预设参数.
% 有数据的年份月份列表(时间区间：2012/07/02-2019/12/31).
yearList = 2012 : 2019;

% 重采样时, 纯像元阈值.
lcThreshold = 0.6;
snowThreshold = 0.6;

% LC和分区影像代码的重分类. 主要针对水域, 冰川, 建筑这三个LC区中LST的反演.
%   [阶梯1_东北东部(1): 1-3, 阶梯1_华北(2): 4-5, 阶梯1_华南(3): 6-10]
%   [阶梯2_西南(4): 11-14, 阶梯2_西北东部(5): 15-25, 阶梯2_东北西部(6): 26-29]
%   [阶梯2_西北西部(7): 30-46， 阶梯3_青藏高原(8): 47-62]
% 重分类时, 使用了右开区间, 导致最后一个分区编码71排除在分类结果中, 62+1, 使其包含在重分类结果中.
regionNodes = [1, 4, 6, 11, 15, 26, 30, 47, 62+1]; 

%% 路径.
% 根目录.
rootPath = 'E:\AMSR_MODIS_Fusion';
dataPath = fullfile(rootPath, 'Data');
funcPath = fullfile(rootPath, 'Code\Functions');
addpath(funcPath);

% 创建保存环境变量综合分类体系数据的文件夹.
ccsevPath = fullfile(dataPath, 'CCSEV_Matlab');
if ~exist(ccsevPath, 'dir')
    mkdir(ccsevPath)
end

% 分区, 地表覆盖, AMSR2 BT示例数据路径, 用于获取数据的属性信息和空间参考.
%   zones和desert长期不变, MODIS地表覆盖每年一幅, MODIS积雪每日一幅.
zonesPath = fullfile(dataPath, 'Zones', 'GeographicalZones_62_Merged.tif');
desertPath = fullfile(dataPath, 'Landcover', 'Desert_CN_gcs.tif');
modisLcPath = fullfile(dataPath, 'MCD12Q1_2_CN_TIF', 'MCD12Q1.A2012001.006_gcs_reclsfy_resmpl.tif');
modisSnowPath = fullfile(dataPath, 'MYD10C1_2_CN_TIF/MYD10C1_2012XXX', ...
    'MYD10C1.A2012002.006.2016081154447.tif');
amsr2BtPath = fullfile(dataPath, 'AMSR2_3_MaskCn_TIF/AMSR2_2012XXXX_TIF/AMSR2_20120703', ...
    'GW1AM2_20120703_01D_EQMA_L3SGT06HA2220220_BtH.tif');

%% 读取数据参考与属性.
% 读取掩膜后的AMSR2 BT图像的空间参考. 所有掩膜后AMSR2 BT与MODIS LST影像的空间参考都相同.
amsr2Ref = georasterinfo(amsr2BtPath).RasterReference;
amsr2RowN = amsr2Ref.RasterSize(1); amsr2ColN = amsr2Ref.RasterSize(2);

% 读取MODIS LC图像的空间参考. 沙漠图像与土地覆盖图像的空间参考相同.
modisLcRef = georasterinfo(modisLcPath).RasterReference;
lcRowN = modisLcRef.RasterSize(1); lcColN = modisLcRef.RasterSize(2);

% 读取MODIS Snow图像的空间参考.
modisSnowRef = georasterinfo(modisSnowPath).RasterReference;
snowRowN = modisSnowRef.RasterSize(1); snowColN = modisSnowRef.RasterSize(2);

% 获取LC和Snow影像左上角第一个滑动窗口的边界行列号, 行列数.
% 滑动窗口的大小为AMSR2影像的像元尺寸. 滑动窗口的行列数是LC影像和AMSR2影像的分辨率倍数. 
% 边界行列号: [topRow, bottomRow, leftCol, rightCol],  行列数: [blockRowN, blockColN].
[snowBlockBdy, snowBlockSize] = getStartBlockRowCol(modisSnowRef, amsr2Ref);
[lcBlockBdy, lcBlockSize] = getStartBlockRowCol(modisLcRef, amsr2Ref);
lcBlockN = prod(lcBlockSize);

%% 合并地表覆盖、沙漠、分区和积雪分区影像.
% 获取沙漠影像和编码. [110, ..., 200]
desertLayer = readgeoraster(desertPath);
desertNodata = georasterinfo(desertPath).MissingDataIndicator;  % Nodata是255.
desertCodeList = unique(desertLayer);
desertCodeList = uint16(desertCodeList(desertCodeList ~= desertNodata));
desertCodeListN = length(desertCodeList);

% 获取地形分区影像.
zonesLayer = readgeoraster(zonesPath);
zonesNodata = georasterinfo(zonesPath).MissingDataIndicator;  % Nodata是128.

% 获取LC编码. [0: 其他(裸地: 100, 植被: 200), 水体: 300, 冰川: 400, 建筑: 500, 积雪: 600, 混合: 1000]
lcLayer = uint16(readgeoraster(modisLcPath)) * 100;
lcNodata = georasterinfo(modisLcPath).MissingDataIndicator * 100;
lcCodeList = unique(lcLayer);
lcCodeList = lcCodeList(lcCodeList ~= lcNodata);
extensiveLcCodeList = lcCodeList(1:2); % 裸地, 植被.
scatteredLcCodeList = lcCodeList(3:5);  % 水体, 冰川, 建筑.
lcCodeList = [0; scatteredLcCodeList; 600; 1000];
lcCodeListN = length(lcCodeList);

% 分年度合并LC, 并升尺度到AMSR2数据的分辨率.
for i = 1 : length(yearList)
    yearStr = num2str(yearList(i));

    % 保存环境变量综合分类体系数据的路径.
    outLcPctArrayMatPath = fullfile(ccsevPath, ['CCSEV_', yearStr, '.mat']);
    if exist(outLcPctArrayMatPath, 'file')
        continue
    end
    
    fprintf('创建环境变量综合分类体系(CCSEV)年份: %s\n', yearStr)
    % 获取地表覆盖类型影像. 移除原数据中的裸地和植被类型, 然后将剩余的地表覆盖类型与沙漠叠置. 当水体, 冰雪,
    %   建筑与沙漠有重叠时, 保持他们不变.
    modisLcName = sprintf('MCD12Q1.A%s001.006_gcs_reclsfy_resmpl.tif', yearStr);
    modisLcPath = fullfile(dataPath, 'MCD12Q1_2_CN_TIF', modisLcName);
    lcLayer = uint16(readgeoraster(modisLcPath)) * 100;
    lcLayer(ismember(lcLayer, extensiveLcCodeList)) = 0; % 将裸地, 植被编码替换为其他.
    desertIndex = (desertLayer ~= desertNodata) & ~ismember(lcLayer, scatteredLcCodeList);
    lcLayer(desertIndex) = desertLayer(desertIndex);

    % 叠置地表覆盖影像和分区影像. 除了积雪之外的其他地表覆盖类型与分区影像在此处结合. 因为积雪覆盖的面积每天
    %   都在变化, 所以积雪覆盖在下一步处理. 混合像元定义为像元内没有任何地表覆盖类型的面积超过60%的像元. 合
    %   并的地表覆盖编码定义为: 分区代码 + 地表覆盖代码. 例如: 310表示在分区10的水域(300).
    [othersPctLayer, desertPctLayer, waterPctLayer, glacierPctLayer, buildingPctLayer] = ...
        deal(zeros(amsr2RowN, amsr2ColN, 'single') * nan);
    desertCodeLayer = zeros(amsr2RowN, amsr2ColN, 'uint8');
    zonesLcLayer = uint16(zonesLayer);
    for ii = 1 : amsr2RowN
        for jj = 1 : amsr2ColN
            % 跳过Zonal图层的非中国区域像元.
            if zonesLcLayer(ii, jj) == zonesNodata
                continue
            end
            % 定位每一个LC滑动窗口.
            blockTopRow = lcBlockBdy(1) + lcBlockSize(1) * (ii - 1);
            blockBottomRow = lcBlockBdy(2) + lcBlockSize(1) * (ii - 1);
            blockLeftCol = lcBlockBdy(3) + lcBlockSize(2) * (jj - 1);
            blockRightCol = lcBlockBdy(4) + lcBlockSize(2) * (jj - 1);
            if blockRightCol > lcColN || blockBottomRow > lcRowN
                continue
            end
            % 窗口计算.
            lcBlock = reshape(lcLayer(blockTopRow:blockBottomRow,blockLeftCol:blockRightCol), [],1);
            othersPctLayer(ii, jj) = sum(lcBlock == lcCodeList(1)) / lcBlockN;
            desertPctLayer(ii, jj) = sum(ismember(lcBlock, desertCodeList)) / lcBlockN;
            waterPctLayer(ii, jj) = sum(lcBlock == lcCodeList(2)) / lcBlockN;
            glacierPctLayer(ii, jj) = sum(lcBlock == lcCodeList(3)) / lcBlockN;
            buildingPctLayer(ii, jj) = sum(lcBlock == lcCodeList(4)) / lcBlockN;
            lcPctList = [othersPctLayer(ii, jj), desertPctLayer(ii, jj), waterPctLayer(ii, jj), ...
                glacierPctLayer(ii, jj), buildingPctLayer(ii, jj)]';

            % 如果单个窗口(1个AMSR2像元)内的某土地覆盖类型面积超过窗口面积的60%，则将其视为纯像元, 否则视为
            %   混合像元.
            maxLcPct = max(lcPctList);
            maxIndex = find(lcPctList == maxLcPct);
            if length(maxIndex) == 1 % 最大值数量只有1个时, 单个土地覆盖类型面积比例有超过0.6的可能.
                if maxLcPct >= lcThreshold % 纯像元.
                    % maxIndex == 1 为其他地表覆盖类型(植被, 裸地), 不做处理.
                    if maxIndex == 2 % 沙漠.
                        zonesLcLayer(ii, jj) = mode(lcBlock); % 若众数个数大于1, mode只取较小的众数.
                    elseif maxIndex >= 3 % 水域, 冰川, 建筑.
                        zonesLcLayer(ii, jj) = zonesLcLayer(ii, jj) + mode(lcBlock);
                    end
                else  % 混合像元.
                    zonesLcLayer(ii, jj) = zonesLcLayer(ii, jj) + lcCodeList(6);
                end
            else % 最大值数量超过1个时, 肯定是混合像元.
                zonesLcLayer(ii, jj) = zonesLcLayer(ii, jj) + lcCodeList(6);
            end

            % 获取与AMSR2影像分辨率相同的沙漠编码数据层, 当窗口内存在多个沙漠分区, 取面积最大的沙漠分区为该
            %   像元的沙漠分区.
            if desertPctLayer(ii, jj) > 0
                desertPixelCount = zeros(desertCodeListN, 1);
                for n = 1 : desertCodeListN
                    desertPixelCount(n) = sum(lcBlock == desertCodeList(n));
                end
                maxDesertPixelCount = max(desertPixelCount);
                majorDesertIndex = find(desertPixelCount == maxDesertPixelCount, 1);
                desertCodeLayer(ii, jj) = desertCodeList(majorDesertIndex);
            end
        end
    end

    % 读取每天的MODIS积雪图像, 将其混合到分区地表覆盖图像中，并修正土地覆盖百分比.
    modisSnowYearPath = fullfile(dataPath, 'MYD10C1_2_CN_TIF', ['MYD10C1_', yearStr, 'XXX']);
    modisSnowDailyFileList = dir(fullfile(modisSnowYearPath, 'MYD10C1*.tif'));
    modisSnowDailyFileList = {modisSnowDailyFileList.name}';
    modisSnowDailyFileListN = length(modisSnowDailyFileList);
    lcDateList = cell(modisSnowDailyFileListN, 1);
    
    [otherPctArray, desertPctArray, glacierPctArray, waterPctArray, snowPctArray, ...
        buildingPctArray] = deal(zeros(amsr2RowN, amsr2ColN, modisSnowDailyFileListN, 'single'));
    zonesLcArray = zeros(amsr2RowN, amsr2ColN, modisSnowDailyFileListN, 'uint16');
    snowZonesLayer = uint16(zonesLayer) + lcCodeList(5); % snowZonesLayer表示积雪所在的分区.
    for j = 1 : modisSnowDailyFileListN
        modisSnowDailyFile = modisSnowDailyFileList{j};
        fprintf('修正积雪覆盖比例: %s\n', modisSnowDailyFile);

        % 获取有MODIS积雪覆盖数据的日期.
        yearDateStr = split(modisSnowDailyFile, '.');
        lcDateList{j} = yday2ymd(yearDateStr{2}(2:end));
        
        % 获取积雪比例数据层.
        modisSnowLayer = readgeoraster(fullfile(modisSnowYearPath, modisSnowDailyFile));
        snowPctLayer = zeros(amsr2RowN, amsr2ColN, 'single');
        for ii = 1 : amsr2RowN
            for jj = 1 : amsr2ColN
                % 跳过zonesLcLayer的非中国区域像元.
                if zonesLcLayer(ii, jj) == zonesNodata
                    continue
                end
                % 定位每一个snow滑动窗口.
                blockTopRow = snowBlockBdy(1) + snowBlockSize(1) * (ii - 1);
                blockBottomRow = snowBlockBdy(2) + snowBlockSize(1) * (ii - 1);
                blockLeftCol = snowBlockBdy(3) + snowBlockSize(2) * (jj - 1);
                blockRightCol = snowBlockBdy(4) + snowBlockSize(2) * (jj - 1);
                if blockRightCol > snowColN || blockBottomRow > snowRowN
                    continue
                end
                % 积雪滑动窗口计算.
                snowBlock = modisSnowLayer(blockTopRow:blockBottomRow, blockLeftCol:blockRightCol);
                snowPixelIndex = (snowBlock >= 0) & (snowBlock <= 100);
                snowPixelN = sum(snowPixelIndex, 'all');
                snowPctLayer(ii, jj) = sum(snowBlock(snowPixelIndex)) / (100 * snowPixelN); 

                % 如果原沙漠像元有一半以上面积由积雪覆盖, 则将此像元视为非沙漠像元.
                % 这一步考虑在反演AMSR2 LST的代码中加入.
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
        zonesLcLayer(pureSnowIndexLayer) = snowZonesLayer(pureSnowIndexLayer);
        zonesLcArray(:, :, j) = zonesLcLayer;

        % 修正地表覆盖比例数据层.
        otherPctArray(:, :, j) = othersPctLayer .* (1 - snowPctLayer);
        desertPctArray(:, :, j) = desertPctLayer .* (1 - snowPctLayer);
        waterPctArray(:, :, j) = waterPctLayer .* (1 - snowPctLayer);
        glacierPctArray(:, :, j) = glacierPctLayer .* (1 - snowPctLayer);
        buildingPctArray(:, :, j) = buildingPctLayer .* (1 - snowPctLayer);
    end

    % LC和分区影像代码的重分类.
    %   [阶梯1_东北东部(1): 1-3, 阶梯1_华北(2): 4-5, 阶梯1_华南(3): 6-10]
    %   [阶梯2_西南(4): 11-14, 阶梯2_西北东部(5): 15-25, 阶梯2_东北西部(6): 26-29]
    %   [阶梯2_西北西部(7): 30-46， 阶梯3_青藏高原(8): 47-62]
    fprintf('中国8大地形分区: %s\n', yearStr);
    zonesLcCodeList = unique(zonesLcArray);
    zonesLcCodeList = zonesLcCodeList(zonesLcCodeList ~= zonesNodata);
    fixedZonesLcArray = zonesLcArray;
    for j = 1 : length(regionNodes) - 1
        for k = 2 : lcCodeListN - 1
            lcIndexVector = (zonesLcCodeList >= lcCodeList(k) + regionNodes(j)) & ...
                (zonesLcCodeList < lcCodeList(k) + regionNodes(j+1));
            zonesLcCodeList(lcIndexVector) = lcCodeList(k) + j;
            lcIndexArray = (zonesLcArray >= lcCodeList(k) + regionNodes(j)) & ...
                (zonesLcArray < lcCodeList(k) + regionNodes(j+1));
            fixedZonesLcArray(lcIndexArray) = lcCodeList(k) + j;
        end
    end
    zonesLcCodeList = unique(zonesLcCodeList);

    % desertCodeLayer 受积雪变化影响, 每天也在变化, 需修正后才能用. 这一步考虑在反演AMSR2 LST时加入.
    fprintf('保存%s年的CCESV.\n', yearStr)
    save(outLcPctArrayMatPath, 'zonesLcLayer', 'otherPctArray', 'desertPctArray', ...
        'waterPctArray', 'glacierPctArray', 'buildingPctArray', 'snowPctArray', ...
        'fixedZonesLcArray', 'desertCodeLayer', 'zonesLcCodeList', 'lcCodeList', 'lcDateList');
end
