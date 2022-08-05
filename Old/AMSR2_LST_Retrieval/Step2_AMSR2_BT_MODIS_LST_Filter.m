%%
%{
2019/01/10 说明
  本实验将moddispixelpercentage阈值设置为0.6, 导致东北等地的AMSR2和MODIS像素匹配度较高的地方很少, 无法构建
回归模型. 采用类似土地覆被区域的回归模型的策略被证实是不恰当的, 因为它们会导致巨大的误差.
  计划使用一个新的根据景观的异质性来设置阈值策略. 在土地覆被类型均匀的缓坡或平坦地区(如东北地区), 阈值可以设
置的越小越好, 而在陡峭的地形或高度异质性的地区, 阈值可以设置的越大.
  由于AMSR2 BT图像没有应的QC层来控制质量, 因此在检索的AMSRE LST图像中会出现大量的异常像素. 设法想出处理这个
问题的计策. 利用AMSR2 LST块像素的STD或一个像素的时间序列.
%}

%% 标识与预设参数.
% 指定的轨道的标识. 1标识升轨(A, 白天), 2标识降轨(D, 晚上).
flg1 = 1;

orbit = {'A', 'D'};
orbit = orbit{flg1};
dayNight = {'Day', 'Night'};
dayNight = dayNight{flg1};

% 年份列表.
yearList = 2012 : 2019;
yearListN = length(yearList);

% 通道列表.
channelList = {'06', '07', '10', '18', '23', '36', '89'};
channelListN = length(channelList);

% 标记MODIS LST像元质量的QC值, 包括: [0, 2, 3, 8, 17, 25, 65, 73, 81, 89, 129, 137, 145, 153].
% QC值代表的含义见文件: 'F:\PaperFusionLST\Doc\MYD11A1_QC.csv'
qc = [0, 17];

% 筛选有效AMSR2像元的阈值, 定义为一个AMSR2像元范围内的MODIS像元的百分比.
modisPixelPercentage = 0.6;

% AMSR2 BT影像中的NoData.
amsr2Nodata = [65534, 65535];

%% 路径.
rootPath = 'E:\AMSR_MODIS_Fusion\';
codeFunctionPath = fullfile(rootPath, 'Code', 'Functions');
% codeFunctionPath = fullfile('D:\AMSR_MODIS_Fusion\Code\Functions');
dataPath = fullfile(rootPath, 'Data');
addpath(codeFunctionPath);

% AMSR2 BT, MODIS LST, SRTM DTM文件夹名.
modisLstPrjCnFolder = 'MYD11A1_2_PrjCN_TIF';
amsr2CnTifFolder = 'AMSR2_2_CN_TIF';
amsr2CnMatFolder = 'AMSR2_2_CN_Matlab';

% 输出数据路径, 同日期经过空间交集掩膜的AMSR2 BT和MODIS LST数据.
amsr2MaskTifFolder = 'AMSR2_3_MaskCn_TIF';
amsr2MaskTifPath = fullfile(dataPath, amsr2MaskTifFolder);
if ~exist(amsr2MaskTifPath, 'dir')
    mkdir(amsr2MaskTifPath)
end

amsr2MaskMatFolder = 'AMSR2_3_MaskCn_Matlab';
amsr2MaskMatPath = fullfile(dataPath, amsr2MaskMatFolder);
if ~exist(amsr2MaskMatPath, 'dir')
    mkdir(amsr2MaskMatPath)
end

modisLstMaskTifFolder = 'MYD11A1_3_MaskCn_TIF';
modisLstMaskTifPath = fullfile(dataPath, modisLstMaskTifFolder);
if ~exist(modisLstMaskTifPath, 'dir')
    mkdir(modisLstMaskTifPath)
end

modisLstMaskMatFolder = 'MYD11A1_3_MaskCn_Matlab';
modisLstMaskMatPath = fullfile(dataPath, modisLstMaskMatFolder);
if ~exist(modisLstMaskMatPath, 'dir')
    mkdir(modisLstMaskMatPath)
end

%% 读取数据参考和属性信息.
% 获取样例AMSR2 BT影像的空间参考信息.
[~, amsr2Ref] = readgeoraster(fullfile(dataPath, amsr2CnTifFolder, 'L3.TB6GHz_10\2012\07\', ...
    'GW1AM2_20120700_01M_EQMA_L3SGT06HA2220220_BtH.tif'));
amsr2RowN = amsr2Ref.RasterSize(1);
amsr2ColN = amsr2Ref.RasterSize(2);
amsr2CellsizeX = amsr2Ref.CellExtentInLongitude;
amsr2CellsizeY = amsr2Ref.CellExtentInLatitude;
amsr2XMin = amsr2Ref.LongitudeLimits(1);
amsr2XMax = amsr2Ref.LongitudeLimits(2);
amsr2YMin = amsr2Ref.LatitudeLimits(1);
amsr2YMax = amsr2Ref.LatitudeLimits(2);
amsr2PixelLeftSideXVector = amsr2XMin : amsr2CellsizeX : amsr2XMax - amsr2CellsizeX;
amsr2PixelTopSideYVector = amsr2YMax : - amsr2CellsizeY : amsr2YMin + amsr2CellsizeY;

% 获取样例MODIS LST影像的空间参考信息.
[~, modisLstRef] = readgeoraster(fullfile(dataPath, modisLstPrjCnFolder, 'MYD11A1_2002XXX_TIF', ...
    'MYD11A1.A2002185.LST_Day.tif'));
modisLstCellsizeX = modisLstRef.CellExtentInLongitude;
modisLstCellsizeY = modisLstRef.CellExtentInLatitude;
modisLstXMin = modisLstRef.LongitudeLimits(1);
modisLstXMax = modisLstRef.LongitudeLimits(2);
modisLstYMin = modisLstRef.LatitudeLimits(1);
modisLstYMax = modisLstRef.LatitudeLimits(2);
modisLstPixelLeftSideXVector = modisLstXMin : modisLstCellsizeX : modisLstXMax - modisLstCellsizeX;
modisLstPixelTopSideYVector = modisLstYMax : - modisLstCellsizeY : modisLstYMin + modisLstCellsizeY;

% 读取MODIS像元面积和SRTM坡度数据.
[srtmSlpLayer, srtmRef] = readgeoraster(fullfile(dataPath, 'SRTM', 'SRTM_0d01_CN_Slp.tif'));
modisAreaLayer = readgeoraster(fullfile(dataPath, modisLstPrjCnFolder, 'MYD11A1.PixelArea.tif'));

% 从MODIS, SRTM影像中获取与中国区AMSR2影像边界对齐的起始和结束的行号和列号.
[modisLStartRow, modisLEndRow, modisLStartCol, modisLEndCol] = getBdyRowCol(modisLstRef, amsr2Ref);
[srtmStartRow, srtmEndRow, srtmStartCol, srtmEndCol] = getBdyRowCol(srtmRef, amsr2Ref);

% 获取中国区域的范围栅格, MODIS像元面积影像, SRTM坡度数据和.
cnMaskLayer = readgeoraster(fullfile(dataPath, 'Zones', 'GeographicalZones_merge.tif')) ~= 128;

modisAreaLayer = modisAreaLayer(modisLStartRow:modisLEndRow, modisLStartCol:modisLEndCol);
modisAreaLayer = double(modisAreaLayer) * 0.01;

srtmSlpLayer = srtmSlpLayer(srtmStartRow:srtmEndRow, srtmStartCol:srtmEndCol);
cosGamaLayer = cosd(srtmSlpLayer);

% 获取用于升尺度MODIS像元的移动滑块的行列数.
blockColN = round(amsr2CellsizeX / modisLstCellsizeX);
blockRowN = round(amsr2CellsizeY / modisLstCellsizeY);
blockN = blockColN * blockRowN;

% 获取与中国区范围AMSR2影像左上角边界对齐的滑动块的左右边界列号与上下边界行号.
itemsPixelLeftSideXVector = modisLstPixelLeftSideXVector(modisLStartCol : modisLEndCol);
itemsPixelTopSideYVector = modisLstPixelTopSideYVector(modisLStartRow : modisLEndRow);
itemsColN = length(itemsPixelLeftSideXVector);
itemsRowN = length(itemsPixelTopSideYVector);

for i = 1 : amsr2RowN
    leftSideXDistance = abs(itemsPixelLeftSideXVector - amsr2PixelLeftSideXVector(i));
    minLeftSideDistance = min(leftSideXDistance);
    if minLeftSideDistance <= modisLstCellsizeX
        itemsStartLeftSideCol = find(leftSideXDistance == minLeftSideDistance);
        break;
    end
end
itemsStartRightSideCol = itemsStartLeftSideCol + blockColN - 1;

for i = 1 : amsr2ColN
    topSideYDistance = abs(itemsPixelTopSideYVector - amsr2PixelTopSideYVector(i));
    minTopSideYDistance = min(topSideYDistance);
    if minTopSideYDistance <= modisLstCellsizeY
        itemsStartTopSideRow = find(topSideYDistance == minTopSideYDistance);
        break;
    end
end
itemsStartBottomSideCol = itemsStartTopSideRow + blockRowN - 1;

%% 读取各年份的数据, 筛选相同日期数据, 升尺度MODIS LST, 并输出空间交集数据.
% 分年份读取MODIS LST和AMSR2 BT的文件.
for i = 1 : yearListN
    yearStr = num2str(yearList(i));
    disp(['处理年份:', yearStr]);

    % 输出的掩膜后年度MODIS LST和AMSR2 BT的矩阵的mat文件路径.
    modisLstMaskMatFileName = ['MYD11A1_MaskCn_', yearStr, '_', dayNight, '.mat'];
    modisLstMaskMatFilePath = fullfile(modisLstMaskMatPath, modisLstMaskMatFileName);
    amsr2MaskMatFileName = ['AMSR2_BT_MaskCn_', yearStr, '_', dayNight, '.mat'];
    amsr2MaskMatFilePath = fullfile(amsr2MaskMatPath, amsr2MaskMatFileName);
    modisMatFileExist = exist(modisLstMaskMatFilePath, 'file');
    amsr2MatFileExist = exist(amsr2MaskMatFilePath, 'file');

    % 读取MODIS LST, QC以及Emis 31, 32的文件列表.
    modisLstYearFolder = ['MYD11A1_', yearStr, 'XXX_TIF'];
    modisLstPrjYearPath = fullfile(dataPath, modisLstPrjCnFolder, modisLstYearFolder);

    modisLstYearList = dir(fullfile(modisLstPrjYearPath, ['*LST_', dayNight, '.tif']));
    modisLstYearList = {modisLstYearList.name}';
    modisLstYearListN = length(modisLstYearList);

    modisQcYearList = dir(fullfile(modisLstPrjYearPath, ['*QC_', dayNight, '.tif']));
    modisQcYearList = {modisQcYearList.name}';
    modisQcYearListN = length(modisQcYearList);

    modisEmis31YearList = dir(fullfile(modisLstPrjYearPath, '*Emis_31.tif'));
    modisEmis31YearList = {modisEmis31YearList.name}';
    modisEmis31YearListN = length(modisEmis31YearList);

    modisEmis32YearList = dir(fullfile(modisLstPrjYearPath, '*Emis_32.tif'));
    modisEmis32YearList = {modisEmis32YearList.name}';
    modisEmis32YearListN = length(modisEmis32YearList);

    % 判断MODIS LST, QC以及Emiss文件是否一一对应.
    if ~ isequal(modisLstYearListN, modisQcYearListN, modisEmis31YearListN, modisEmis32YearListN)
        error('%s年的MODIS LST, QC, Emis31, Emis32数据不匹配, 请检查!', yearStr)
    end
    [modisLstDateList, modisQcDateList, modisEmis31DateList, modisEmis32DateList] = ...
        deal(cell(modisLstYearListN, 1));
    for j = 1 : modisLstYearListN
        modisLstDate = split(modisLstYearList{j}, '.');
        modisLstDateList{j} = yday2ymd(modisLstDate{2}(2:end));
        modisQcDate = split(modisQcYearList{j}, '.');
        modisQcDateList{j} = yday2ymd(modisQcDate{2}(2:end));
        modisEmis31Date = split(modisEmis31YearList{j}, '.');
        modisEmis31DateList{j} = yday2ymd(modisEmis31Date{2}(2:end));
        modisEmis32Date = split(modisEmis32YearList{j}, '.');
        modisEmis32DateList{j} = yday2ymd(modisEmis32Date{2}(2:end));
    end
    if  ~ isequal(modisLstDateList, modisQcDateList, modisEmis31DateList, modisEmis32DateList)
        error('%s年的MODIS LST与QC数据不匹配, 请检查!', yearStr)
    end

    % 从mat文件中读取AMSR2 BT数据各通道文件列表, 数据本身及空间参考.
    % 变量 amsr2PathMatrix 为存储AMSR2 BT TIF文件路径字符串的cell矩阵. 矩阵的行表示不同的日期, 列表示不同
    %   的极化通道. 第一行为表头, 存储2极化和7通道的名称: ['10H', '10V', '18H', '18V', '23H', '23V',
    %   '36H', '36V', '6H', '6V', '7H', '7V', '89H', '89V']. 其他行依次存储这些极化通道对应的不同日期的
    %   AMSR2 BT TIF文件的路径字符串.
    amsr2CnMatName = sprintf('AMSR2_BT_%s_%s.mat', yearStr, dayNight);
    amsr2CnMatPath = fullfile(dataPath, amsr2CnMatFolder, amsr2CnMatName);
    if ~exist(amsr2CnMatPath, "file")
        error('%s %s AMSR2 BT Mat文件不存在, 请检查!', yearStr, dayNight);
    end
    load(amsr2CnMatPath, 'amsr2PathMatrix', 'amsr2Ref', 'amsr2H10YearArray', ...
        'amsr2V10YearArray', 'amsr2H18YearArray', 'amsr2V18YearArray', 'amsr2H23YearArray', ...
        'amsr2V23YearArray', 'amsr2H36YearArray', 'amsr2V36YearArray', 'amsr2H06YearArray', ...
        'amsr2V06YearArray', 'amsr2H07YearArray', 'amsr2V07YearArray', 'amsr2H89YearArray', ...
        'amsr2V89YearArray')
    amsr2PathMatrix = replace(amsr2PathMatrix, 'F:\AMSR_MODIS_Fusion\Data\', dataPath);

    % 输入AMSR2 BT数据日期列表.
    cpListN = length(amsr2PathMatrix(1, :));
    amsr2PathMatrix = amsr2PathMatrix(2:end, :);
    [~, amsr2Name, ~] = fileparts(amsr2PathMatrix);
    amsr2DateMatrix = split(amsr2Name, '_');
    amsr2DateList = amsr2DateMatrix(:, 1, 2);

    % 获取MODIS LST和AMSR2 BT均有数据日期的文件路径.
    [~, modisDateIndex, amsr2DateIndex] = intersect(modisLstDateList, amsr2DateList);
    dateFilterN = length(modisDateIndex);
    modisLstPathFilterList = cellfun(@fullfile, repmat({modisLstPrjYearPath}, dateFilterN, 1), ...
        modisLstYearList(modisDateIndex), 'UniformOutput', false);
    modisQcPathFilterList = cellfun(@fullfile, repmat({modisLstPrjYearPath}, dateFilterN, 1), ...
        modisQcYearList(modisDateIndex), 'UniformOutput', false);
    modisEmis31PathFilterList = cellfun(@fullfile, repmat({modisLstPrjYearPath}, dateFilterN, 1),...
        modisEmis31YearList(modisDateIndex), 'UniformOutput', false);
    modisEmis32PathFilterList = cellfun(@fullfile, repmat({modisLstPrjYearPath}, dateFilterN, 1),...
        modisEmis32YearList(modisDateIndex), 'UniformOutput', false);
    amsr2PathFilterMatrix = amsr2PathMatrix(amsr2DateIndex, :);
    dateFilterList = amsr2DateList(amsr2DateIndex);

    % 创建存储掩膜后MODIS LST和AMSR2 BT数据的年路径.
    modisLstMaskTifYearPath = fullfile(dataPath, modisLstMaskTifFolder, modisLstYearFolder);
    if ~exist(modisLstMaskTifYearPath, 'dir')
        mkdir(modisLstMaskTifYearPath)
    end

    amsr2MaskYearFolder = ['AMSR2_', yearStr, 'XXXX_TIF'];
    amsr2MaskTifYearPath = fullfile(dataPath, amsr2MaskTifFolder, amsr2MaskYearFolder);
    if ~exist(amsr2MaskTifYearPath, 'dir')
        mkdir(amsr2MaskTifYearPath)
    end

    % 创建存储掩膜后年度MODIS LST和AMSR2 BT的矩阵.
    if ~modisMatFileExist || ~amsr2MatFileExist
        modisLstMaskYearArray = zeros(amsr2RowN, amsr2ColN, dateFilterN, 'single') * nan;
        [amsr2HMaskYearArray, amsr2VMaskYearArray] = ...
            deal(zeros(amsr2RowN, amsr2ColN, channelListN * dateFilterN, 'uint16'));
        startN = 1;
    end

    % 读取MODIS LST和AMSR2 BT数据, 并掩膜输出.
    for j = 1 : dateFilterN
        % 掩膜后的MODIS LST文件输出路径.
        modisLstMaskFileName = split(modisLstPathFilterList{j}, '\');
        modisLstMaskFileName = modisLstMaskFileName{end};
        modisLstMaskFilePath = fullfile(modisLstMaskTifYearPath, modisLstMaskFileName);

        % 掩膜后的各通道AMSR2 BT文件输出路径.
        amsr2MaskYearDateFolder = ['AMSR2_', dateFilterList{j}];
        amsr2MaskYearDatePath = fullfile(amsr2MaskTifYearPath, amsr2MaskYearDateFolder);
        if ~exist(amsr2MaskYearDatePath, 'dir')
            mkdir(amsr2MaskYearDatePath);
        end
        amsr2MaskFileNameList = split(amsr2PathFilterMatrix(j, :), '\');
        amsr2MaskFileNameList = amsr2MaskFileNameList(:, :, end);
        amsr2MaskFilePathList = fullfile(amsr2MaskYearDatePath, amsr2MaskFileNameList);

        % 输出的TIF数据是否存在的标记.
        amsr2MaskFileExistList = false(cpListN, 1);
        for k = 1 : cpListN
            amsr2MaskFileExistList(k) = exist(amsr2MaskFilePathList{k}, 'file');
        end
        amsr2MaskFileN = sum(logical(amsr2MaskFileExistList));

        % 在输出的掩膜后TIF文件存在的前提下, 如果Mat文件也存在, 则直接进入下次循环, 否则将当前日期数据存储
        %   在年度矩阵中. 如果掩膜后的TIF文件不存在, 则Mat文件肯定不存在, 此时进行掩膜处理, 并保存TIF文件和
        %   Mat文件.
        if exist(modisLstMaskFilePath, 'file') && (amsr2MaskFileN == cpListN)
            if modisMatFileExist && amsr2MatFileExist
                continue
            else
                modisLstMaskYearArray(:, :, j) = readgeoraster(modisLstMaskFilePath);
                [amsr2HArray, amsr2VArray] = deal(zeros(amsr2RowN, amsr2ColN, channelListN));
                for k = 1 : channelListN
                    amsr2HArray(:, :, k) = readgeoraster(amsr2MaskFilePathList{2*k-1});
                    amsr2VArray(:, :, k) = readgeoraster(amsr2MaskFilePathList{2*k});
                end
                amsr2HMaskYearArray(:, :, startN:startN+channelListN-1) = amsr2HArray;
                amsr2VMaskYearArray(:, :, startN:startN+channelListN-1) = amsr2VArray;
                startN = startN + channelListN;
                continue
            end
        end

        % 读取MODIS LST, QC, Emis31, Emis32数据层.
        modisLstLayer = readgeoraster(modisLstPathFilterList{j});
        modisQcLayer = readgeoraster(modisQcPathFilterList{j});
        modisEmis31Layer = readgeoraster(modisEmis31PathFilterList{j});
        modisEmis32Layer = readgeoraster(modisEmis32PathFilterList{j});

        % 将这些数据裁剪到AMSR2影像定义的中国范围.
        modisLstLayer = modisLstLayer(modisLStartRow:modisLEndRow, modisLStartCol:modisLEndCol);
        modisQcLayer = modisQcLayer(modisLStartRow:modisLEndRow, modisLStartCol:modisLEndCol);
        modisEmis31Layer = modisEmis31Layer(modisLStartRow:modisLEndRow, modisLStartCol:modisLEndCol);
        modisEmis32Layer = modisEmis32Layer(modisLStartRow:modisLEndRow, modisLStartCol:modisLEndCol);

        % 将MODIS LST和发射率都转为实际值.
        modisLstLayer = double(modisLstLayer) * 0.02;
        modisLstLayer(~ismember(modisQcLayer, qc) | (modisLstLayer == 65535*0.02)) = nan;
        modisEmis31Layer = double(modisEmis31Layer) * 0.002 + 0.49;
        modisEmis31Layer(modisEmis31Layer == 255 * 0.002 + 0.49) = nan;
        modisEmis32Layer = double(modisEmis32Layer) * 0.002 + 0.49;
        modisEmis32Layer(modisEmis32Layer == 255 * 0.002 + 0.49) = nan;

        % MODIS LST影像升尺度.
        modisBbeLayer = 0.273 + 1.778 * modisEmis31Layer - 1.807 * modisEmis31Layer .* ...
            modisEmis32Layer - 1.037 * modisEmis32Layer + 1.774 * modisEmis32Layer .^ 2;
        numeratorLayer = modisAreaLayer .* modisBbeLayer .* modisLstLayer .^ 4 .* ...
            secd(srtmSlpLayer) ./ cosGamaLayer;
        numeratorLayer(numeratorLayer < 0) = nan;
        denominatorLayer = modisBbeLayer .* modisAreaLayer .* secd(srtmSlpLayer);
        [numeratorSumLayer, denominatorSumLayer] = deal(zeros(amsr2RowN, amsr2ColN) * nan);
        for ii = 1 : amsr2RowN
            for jj = 1 : amsr2ColN
                % 定位滑动块的位置.
                itemsBlockTopRow = itemsStartTopSideRow + blockRowN * (ii - 1);
                itemsBlockBottomRow = itemsStartBottomSideCol + blockRowN * (ii - 1);
                itemsBlockLeftCol = itemsStartLeftSideCol + blockColN * (jj - 1);
                itemsBlockRightCol = itemsStartRightSideCol + blockColN * (jj - 1);
                if itemsBlockRightCol > itemsColN || itemsBlockBottomRow > itemsRowN
                    continue;
                end

                % 块计算.
                numeratorBlock = numeratorLayer(itemsBlockTopRow : itemsBlockBottomRow, ...
                    itemsBlockLeftCol: itemsBlockRightCol);
                denominatorBlock = denominatorLayer(itemsBlockTopRow : itemsBlockBottomRow, ...
                    itemsBlockLeftCol: itemsBlockRightCol);
                blockAvailableIndex = ~isnan(numeratorBlock) & ~isnan(denominatorBlock);
                availablePixelRatio = sum(blockAvailableIndex(:)) / blockN;
                if availablePixelRatio < modisPixelPercentage
                    continue;
                end
                numeratorBlock(~blockAvailableIndex) = nan;
                denominatorBlock(~blockAvailableIndex) = nan;
                numeratorSumLayer(ii, jj) = sum(numeratorBlock(:), 'omitnan');
                denominatorSumLayer(ii, jj) = sum(denominatorBlock(:), 'omitnan');
            end
        end
        N = ones(amsr2RowN, amsr2ColN) * 4;
        upScaledModisLstLayer = nthroot(numeratorSumLayer ./ denominatorSumLayer, N);

        % AMSR2 BT质量控制, 排除PR值大于1和轨道间隙像元.
        [amsr2HArray, amsr2VArray, amsr2PrArray] = deal(zeros(amsr2RowN, amsr2ColN, channelListN));
        for k = 1 : channelListN % 顺序: {'06', '07', '10', '18', '23', '36', '89'}.
            amsr2HArray(:, :, k) = eval(sprintf('amsr2H%sYearArray(:, :, j)', channelList{k}));
            amsr2VArray(:, :, k) = eval(sprintf('amsr2V%sYearArray(:, :, j)', channelList{k}));
            amsr2PrArray(:, :, k) = amsr2HArray(:, :, k) ./ amsr2VArray(:, :, k);
        end
        amsr2PrIndexArray = amsr2PrArray > 1;
        amsr2PrIndexLayer = logical(sum(amsr2PrIndexArray, 3));
        amsr2NodataIndexArray = ismember(amsr2HArray,amsr2Nodata)|ismember(amsr2VArray,amsr2Nodata);
        amsr2NodataIndexLayer = logical(sum(amsr2NodataIndexArray, 3));
        amsr2NanIndexLayer = amsr2PrIndexLayer | amsr2NodataIndexLayer;
        amsr2NanIndexArray = repmat(amsr2NanIndexLayer, 1, 1, channelListN);
        amsr2HArray(amsr2NanIndexArray) = nan;
        amsr2VArray(amsr2NanIndexArray) = nan;

        % 获取中国范围内AMSR2 BT数据和升尺度后MODIS LST数据中的共有可用像元.
        modisNanIndexLayer = isnan(upScaledModisLstLayer);
        nanIndexLayer = modisNanIndexLayer | amsr2NanIndexLayer | ~cnMaskLayer;
        nanIndexArray = repmat(nanIndexLayer, 1, 1, channelListN);
        upScaledModisLstLayer(nanIndexLayer) = nan;
        amsr2HArray(nanIndexArray) = nan;
        amsr2VArray(nanIndexArray) = nan;

        % 输出质量控制且升尺度后的MODIS LST和AMSR2 BT数据.
        fprintf('输出%s %s的MODIS LST, AMSR2 BT掩膜数据.\n', dateFilterList{j}, dayNight)
        geotiffwrite(modisLstMaskFilePath, single(upScaledModisLstLayer), amsr2Ref, ...
            'TiffTags', struct('Compression','LZW'));
        for k = 1 : channelListN
            geotiffwrite(amsr2MaskFilePathList{2*k-1}, uint16(amsr2HArray(:, :, k)), amsr2Ref, ...
                'TiffTags', struct('Compression','LZW'));
            geotiffwrite(amsr2MaskFilePathList{2*k}, uint16(amsr2VArray(:, :, k)), amsr2Ref, ...
                'TiffTags', struct('Compression','LZW'));
        end

        % 年度掩膜MODIS LST和AMSR2 BT矩阵.
        if ~modisMatFileExist || ~amsr2MatFileExist
            modisLstMaskYearArray(:, :, j) = upScaledModisLstLayer;
            amsr2HMaskYearArray(:, :, startN:startN+channelListN-1) = amsr2HArray;
            amsr2VMaskYearArray(:, :, startN:startN+channelListN-1) = amsr2VArray;
            startN = startN + channelListN;
        end
    end

    % 存储掩膜后年度MODIS LST和AMSR2 BT的矩阵的mat文件.
    if ~modisMatFileExist || ~amsr2MatFileExist
        fprintf('保存%s年的MODIS LST, AMSR2 BT Mat文件.\n', yearStr)
        save(modisLstMaskMatFilePath, 'amsr2Ref', 'modisLstMaskYearArray', 'dateFilterList');
        save(amsr2MaskMatFilePath, 'amsr2Ref', 'amsr2HMaskYearArray', 'amsr2VMaskYearArray', ...
            'dateFilterList');
    end
end