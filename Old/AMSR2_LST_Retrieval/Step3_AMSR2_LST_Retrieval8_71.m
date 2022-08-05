%% AMSR2地表温度反演

%% 功能标记与预设参数
% 指定白天和晚上的标记. 1表示白天(升轨), 2表示晚上(降轨).
flg1 = 1;

% AMSR2数据的通道, 极化, 轨道.
channelList = ["10", "18", "23", "36", "06", "07", "89"];
channelListN = length(channelList);
polarize = ["H", "V"];
polarizeN = length(polarize);

% [pMatrix, cMatrix] = meshgrid(polarize, channelList);
% pcList = cellstr(reshape(pMatrix + cMatrix, [], 1));

% 回归模型中的变量个数(14个通道 + 2个二次项 + 1常数项).
variablesN = channelListN * polarizeN + 3;

% orbit = {'A', 'D'};  % A表示升轨, D表示降轨.
daynight = {'Day', 'Night'};
daynight = daynight{flg1};

% 有数据的年份月份列表(时间区间：2012/07/02-2019/12/31).
% startDate = [2012, 07, 02]; endDate = [2019, 12, 31];
% yearList = startDate(1) : endDate(1);
yearList = 2012 : 2019;
yearListN = length(yearList);

% 各季节[春, 夏, 秋, 冬]包含的月份, 和各月份的名称.
seasonMonthsList = {[3, 4, 5], [6, 7, 8], [9, 10, 11], [12, 1, 2]};
seasonNameList = {'Spring', 'Summer', 'Autumn', 'Winter'};
seasonMonthsListN = length(seasonMonthsList);
monthNameList = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov','Dec'};
monthNameListN = length(monthNameList);

% LC和分区影像代码的重分类.
%   [阶梯1_东北的东部(1): 1-4, 阶梯1_华北(2): 5-8, 阶梯1_华南(3): 9-13]
%   [阶梯2_西南(4): 14-17, 阶梯2_西北的东部(5): 18-27, 阶梯2_东北的西部(6): 28-30]
%   [阶梯2_西北的西部(7)： 31-53， 阶梯3_青藏高原(8): 54-71]
regionNodes = [1, 5, 9, 14, 18, 28, 31, 54, 71+1];
regionsN = length(regionNodes) - 1;

% AMSR2 BT影像中的NoData.
amsr2Nodata = [65534, 65535];

% 地表覆盖类型升尺度时纯像元阈值.
lcThreshold = 0.6;

%% 路径.
% 根目录.
rootPath = 'E:\AMSR_MODIS_Fusion';
dataPath = fullfile(rootPath, 'Data');
figurePath = fullfile(rootPath, 'Figures');
addpath(fullfile(rootPath, 'Code\Functions'));
if ~exist(figurePath, 'dir')
    mkdir(figurePath)
end

% 输入数据路径.
amsr2CnMatPath = fullfile(dataPath, 'AMSR2_2_CN_Matlab');
amsr2MaskMatPath = fullfile(dataPath, 'AMSR2_3_MaskCn_Matlab');
modisLstMaskMatPath = fullfile(dataPath, 'MYD11A1_3_MaskCn_Matlab');
ccsevMatPath = fullfile(dataPath, 'CCSEV_Matlab');

% 输出统计图路径.
lstScatterPath = fullfile(figurePath, 'LstScatter');
if ~exist(lstScatterPath, 'dir')
    mkdir(lstScatterPath)
end

% 输出反演的AMSR2 LST路径.
regressionMatPath = fullfile(dataPath, 'Regression_Matlab');
if ~exist(regressionMatPath, 'dir')
    mkdir(regressionMatPath)
end
amsr2LstPath = fullfile(dataPath, 'AMSR2_4_LstCn_TIF');
if ~exist(amsr2LstPath, 'dir')
    mkdir(amsr2LstPath)
end

%% 回归和输出.
% 分年度反演AMSR2 LST.
for i = 1 : yearListN
    yearStr = num2str(yearList(i));
    yearStr = '2013';

    % 存储纯像元分区回归系数与反演的AMSR2 LST影像数据的Mat文件路径.
    regressionPureMatName = sprintf('Regression_Pure_%s_%s.mat', yearStr, daynight);
    regressionPureMatPath = fullfile(regressionMatPath, regressionPureMatName);
    if  ~exist(regressionPureMatPath, 'file')

        % 创建输出AMSR2 BT(LST)和MODIS LST散点图的文件夹.
        lstScatterYearPath = fullfile(lstScatterPath, sprintf('%s %s', yearStr, daynight));
        if ~exist(lstScatterYearPath, 'dir')
            mkdir(lstScatterYearPath);
        end

        % 从Mat文件中读取中国区原始的AMSR2 BT数据.
        amsr2CnMatFileName = sprintf('AMSR2_BT_%s_%s.mat', yearStr, daynight);
        amsr2CnMatFilePath = fullfile(amsr2CnMatPath, amsr2CnMatFileName);
        load(amsr2CnMatFilePath, 'amsr2H*YearArray', 'amsr2V*YearArray');

        % 从Mat文件中读取Mask后的AMSR2 BT和MODIS LST数据.
        amsr2MaskMatFileName = sprintf('AMSR2_BT_MaskCn_%s_%s.mat', yearStr, daynight);
        amsr2MaskMatFilePath = fullfile(amsr2MaskMatPath, amsr2MaskMatFileName);
        load(amsr2MaskMatFilePath, 'amsr2HMaskYearArray', 'amsr2VMaskYearArray', 'dateFilterList');

        modisLstMaskFileName = sprintf('MYD11A1_MaskCn_%s_%s.mat', yearStr, daynight);
        modisLstMaskFilePath = fullfile(modisLstMaskMatPath, modisLstMaskFileName);
        load(modisLstMaskFilePath, 'amsr2Ref', 'modisLstMaskYearArray');

        % 从Mat文件中读取CCSEV.
        ccsevFileName = sprintf('CCSEV_%s.mat', yearStr);
        ccsevFilePath = fullfile(ccsevMatPath, ccsevFileName);
%         load(ccsevFilePath, 'zonesLcIdList', 'lcCodeList', 'dateList', 'fixedZonesLcArray', ...
%             'otherPctArray', 'waterPctArray', 'glacierPctArray', 'buildingPctArray', ...
%             'snowPctArray', 'desertPctArray', 'desertCodeLayer');
        load(ccsevFilePath, 'zonesLcIdList', 'lcCodeList', 'dateList', 'fixedZonesLcArray', ...
            'desertCodeLayer');

        % 获取共有日期的分区和数据矩阵.
        [validDateList, dateIndex1, dateIndex2] = intersect(dateFilterList, dateList);
        validYearMonthList = datetime(validDateList , 'InputFormat', 'yyyyMMdd').Month;
        validDateListN = length(validDateList);
        fixedZonesLcArray = fixedZonesLcArray(:, :, dateIndex2);
        modisLstMaskYearArray = double(modisLstMaskYearArray(:, :, dateIndex1)); %!!!本来就是double, 不用转!!!

        % 获取AMSR2 BT影像的行列数.
        [amsr2RowN, amsr2ColN] = size(amsr2HMaskYearArray, 1, 2);

        % 获取沙漠分区编码.
        desertCodeList = unique(desertCodeLayer);
        desertCodeList(isnan(desertCodeList)) = [];
        desertCodeListN = length(desertCodeList);

        % 将不同通道的AMSR2 BT数据拆分到单独的矩阵, 并存入元胞数组中.
        % 通道的排序: [10H, 10V, 18H, 18V, 23H, 23V, 36H, 36V, 06H, 06V, 07H, 07V, 89H, 89V]
        [amsr2HMaskYearCell, amsr2VMaskYearCell] = deal(cell(1, channelListN));
        for k = 1 : channelListN
            channelIndexVector = k : channelListN : channelListN * validDateListN;
            amsr2HChannelArray = amsr2HMaskYearArray(:, :, channelIndexVector);
            amsr2VChannelArray = amsr2VMaskYearArray(:, :, channelIndexVector);
%             amsr2HMaskYearCell{k} = double(amsr2HChannelArray(:, :, dateIndex1));
%             amsr2VMaskYearCell{k} = double(amsr2VChannelArray(:, :, dateIndex1));
            amsr2HMaskYearCell{k} = single(amsr2HChannelArray(:, :, dateIndex1));
            amsr2VMaskYearCell{k} = single(amsr2VChannelArray(:, :, dateIndex1));
        end
        clear amsr2HChannelArray amsr2VChannelArray;
        clear amsr2HMaskYearArray amsr2VMaskYearArray;

        % 纯像元和混合像元分别反演AMSR2 LST.
        zonesLcIdListN = length(zonesLcIdList);
        mixedLcN = sum(zonesLcIdList >= lcCodeList(end));
        pureLcIdN = zonesLcIdListN - mixedLcN;

        % 存储评价AMSR2 LST反演精度指标的矩阵.
        [rmseYearVector, nYearVector, r2YearVector] = deal(zeros(pureLcIdN, 1) * nan);
        [rmseSeasonArray, nSeasonArray, r2SeasonArray] = deal(zeros(pureLcIdN, 4) * nan);
        [rmseMonthArray, nMonthArray, r2MonthArray] = deal(zeros(pureLcIdN, 12) * nan);

        % 存储反演AMSR2 LST回归系数的矩阵.
        coefficientYearArray = zeros(pureLcIdN, variablesN);
        coefficientSeasonArray = zeros(pureLcIdN, variablesN, 4);
        coefficientMonthArray = zeros(pureLcIdN, variablesN, 12);

        % 存储分别使用年尺度, 季节尺度, 和月尺度反演后的AMSR2 LST影像数据的矩阵, 包括掩膜区和全中国区.
        [amsr2LstMaskYearArray1, amsr2LstCnYearArray1, amsr2LstMaskYearArray2, ...
            amsr2LstCnYearArray2, amsr2LstMaskYearArray3, amsr2LstCnYearArray3] = ...
            deal(zeros(amsr2RowN, amsr2ColN, validDateListN, 'single'));

        % =========================================================================================
        % 获取反演AMSR2 LST的系数, 精度, 以及掩膜区AMSR2 LST.
        for j = 1 : pureLcIdN
            zonesLcId = zonesLcIdList(j);
            zoneName = sprintf('Zone%d', zonesLcId);
            fprintf('分区%d, %s年.\n', zonesLcId, yearStr);

            % 获取当前分区的位置索引.
            zonesLcIndexArray = (fixedZonesLcArray == zonesLcId);

            % 从一整年的数组中提取当前分区的原始AMSR2 BT影像, 以及掩膜后的AMSR2 BT和MODIS LST影像.
            [amsr2HCnZoneYearCell, amsr2VCnZoneYearCell] = deal(cell(1, channelListN));
            [amsr2HMaskZoneYearCell,amsr2VMaskZoneYearCell] = deal(cell(1, channelListN));
            for k = 1 : channelListN
                amsr2HCnYearArray = single(eval(sprintf('amsr2H%sYearArray', channelList{k})));
                amsr2VCnYearArray = single(eval(sprintf('amsr2V%sYearArray', channelList{k})));
                amsr2HCnYearArray(ismember(amsr2HCnYearArray, amsr2Nodata)) = nan;
                amsr2VCnYearArray(ismember(amsr2VCnYearArray, amsr2Nodata)) = nan;
                amsr2HCnYearArray = amsr2HCnYearArray(:, :, dateIndex1);
                amsr2VCnYearArray = amsr2VCnYearArray(:, :, dateIndex1);
                amsr2HCnZoneYearCell{k} = setnan(amsr2HCnYearArray, ~zonesLcIndexArray) / 100;
                amsr2VCnZoneYearCell{k} = setnan(amsr2VCnYearArray, ~zonesLcIndexArray) / 100;
                amsr2HMaskZoneYearCell{k} = setnan(amsr2HMaskYearCell{k}, ~zonesLcIndexArray)/100;
                amsr2VMaskZoneYearCell{k} = setnan(amsr2VMaskYearCell{k}, ~zonesLcIndexArray)/100;
            end
            amsr2QdCnZoneYearCell = {...
                (amsr2VCnZoneYearCell{4} - amsr2VCnZoneYearCell{3}) .^ 2, ...
                (amsr2VCnZoneYearCell{4} - amsr2VCnZoneYearCell{2}) .^ 2};
            modisLstMaskZoneYearArray = setnan(modisLstMaskYearArray, ~zonesLcIndexArray);
            clear amsr2HCnYearArray amsr2VCnYearArray

            % -------------------------------------------------------------------------------------
            % 按年尺度回归AMSR2 BT和MODIS LST, 并获取反演的AMSR2 LST影像.
            % 按照回归函数的格式要求, 将参与回归的影像矩阵重新排列.
            valueIndex = find(~isnan(modisLstMaskZoneYearArray));
            modisLstMaskZoneYearVector = modisLstMaskZoneYearArray(valueIndex);
            amsr2HMaskZoneYearVector = cell(1, channelListN);
            amsr2VMaskZoneYearVector = cell(1, channelListN);
            for k = 1 : channelListN
                amsr2HMaskZoneYearVector{k} = amsr2HMaskZoneYearCell{k}(valueIndex);
                amsr2VMaskZoneYearVector{k} = amsr2VMaskZoneYearCell{k}(valueIndex);
            end
            amsr2QdMaskZoneYearVector = [...
                (amsr2VMaskZoneYearVector{4} - amsr2VMaskZoneYearVector{3}) .^ 2, ...
                (amsr2VMaskZoneYearVector{4} - amsr2VMaskZoneYearVector{2}) .^ 2];

            amsr2BtMaskZoneYearRecords = double(cell2mat([amsr2HMaskZoneYearVector, ...
                amsr2VMaskZoneYearVector, amsr2QdMaskZoneYearVector]));
            clear amsr2HMaskZoneYearVector amsr2VMaskZoneYearVector amsr2QdMaskZoneYearVector;
            clear amsr2Qd1MaskZoneYearVector amsr2Qd2MaskZoneYearVector;

            % 样本数大于2时才能执行逐步回归, 否则该分区不能反演AMSR2 LST.
            % 系数矩阵 [10H, 18H, 23H, 36H, 06H, 07H, 89H, 10V, 18V, 23V, 36V, 06V, 07V, 89V, 常数].
            modisLstMaskZoneYearVectorN = length(modisLstMaskZoneYearVector);
            if modisLstMaskZoneYearVectorN >= 2
                timestamp = sprintf('%s %s', yearStr, daynight);
                fprintf('分区%d 年尺度 回归: %s.\n', zonesLcId, timestamp)
                mdl = stepwiselm(amsr2BtMaskZoneYearRecords, modisLstMaskZoneYearVector, ...
                    Lower='constant', Upper='linear', Criterion='aic');
                lstRMSE = mdl.RMSE;
                lstR2 = mdl.Rsquared.Ordinary;
                amsr2LstMaskZoneYearVector = mdl.Fitted;
                variableIndex = find(mdl.VariableInfo.InModel == 1);
                pYear = mdl.Coefficients.Estimate;
                pYear = pYear([2:end, 1]);
                coefficientYearArray(j, [variableIndex; end]) = pYear;

                amsr2LstMaskZoneYearVector2 = single(zeros(amsr2RowN*amsr2ColN*validDateListN, 1));
                amsr2LstMaskZoneYearVector2(valueIndex) = amsr2LstMaskZoneYearVector;
                amsr2LstMaskZoneYearArray = reshape(amsr2LstMaskZoneYearVector2, ...
                    amsr2RowN, amsr2ColN, validDateListN);
                clear amsr2LstMaskZoneYearVector2

                % 根据回归系数计算分区内的AMSR2 LST(初始化矩阵为常数项截距).
                amsr2LstCnZoneYearArray = repmat(pYear(end), amsr2RowN, amsr2ColN, validDateListN);
                for k = 1 : channelListN
                    amsr2LstCnZoneYearArray = amsr2LstCnZoneYearArray + ...
                        amsr2HCnZoneYearCell{k} * coefficientYearArray(j, k) + ...
                        amsr2VCnZoneYearCell{k} * coefficientYearArray(j, k + channelListN);
                end
                amsr2LstCnZoneYearArray = amsr2LstCnZoneYearArray + ...
                    amsr2QdCnZoneYearCell{:, 1} * coefficientYearArray(j, k*2+1) + ...
                    amsr2QdCnZoneYearCell{:, 2} * coefficientYearArray(j, k*2+2);
                % nan值改为0, 便于算完所有分区AMSR2 LST后将结果相加.
                amsr2LstCnZoneYearArray(isnan(amsr2LstCnZoneYearArray)) = 0;

                % 输出AMSR2 BT(LST)和MODIS LST的年度散点图.
                zoneScatterYearName = sprintf('LstScatter_%s_%s.tif', zoneName, timestamp);
                zoneScatterYearPath = fullfile(lstScatterYearPath, zoneScatterYearName);
                if  ~exist(zoneScatterYearPath, 'file')
                    f = lstScatter(amsr2LstMaskZoneYearVector, modisLstMaskZoneYearVector, ...
                        timestamp, zoneName, [lstR2, lstRMSE]);
                    exportgraphics(f, zoneScatterYearPath);
                    close all;
                end
            else
                lstRMSE = nan; lstR2 = nan;
                amsr2LstMaskZoneYearVector = zeros(modisLstMaskZoneYearVectorN, 1) * nan;
                amsr2LstMaskZoneYearArray = zeros(amsr2RowN, amsr2ColN, validDateListN);
                amsr2LstCnZoneYearArray = zeros(amsr2RowN, amsr2ColN, validDateListN);
            end
            rmseYearVector(j) = lstRMSE;
            nYearVector(j) = modisLstMaskZoneYearVectorN;
            r2YearVector(j) = lstR2;

            % 将反演的当年AMSR2 LST保存到年度矩阵中.
            amsr2LstMaskYearArray1 = amsr2LstMaskYearArray1 + amsr2LstMaskZoneYearArray;
            amsr2LstCnYearArray1 = amsr2LstCnYearArray1 + amsr2LstCnZoneYearArray;

            % -------------------------------------------------------------------------------------
            % 按季节尺度回归AMSR2 BT和MODIS LST, 并获取反演的AMSR2 LST影像.
            for n = 1 : seasonMonthsListN
                seasonIndex = ismember(validYearMonthList, seasonMonthsList{n});
                seasonIndexN = sum(seasonIndex);
                seasonName = seasonNameList{n};

                % 筛选当前季节的数据矩阵.
                modisLstMaskZoneSeasonArray = modisLstMaskZoneYearArray(:, :, seasonIndex);
                [amsr2HMaskZoneSeasonCell, amsr2VMaskZoneSeasonCell] = deal(cell(1, channelListN));
                [amsr2HCnZoneSeasonCell, amsr2VCnZoneSeasonCell] = deal(cell(1, channelListN));
                for k = 1 : channelListN
                    amsr2HMaskZoneSeasonCell{k} = amsr2HMaskZoneYearCell{k}(:,:,seasonIndex);
                    amsr2VMaskZoneSeasonCell{k} = amsr2VMaskZoneYearCell{k}(:,:,seasonIndex);
                    amsr2HCnZoneSeasonCell{k} = amsr2HCnZoneYearCell{k}(:,:,seasonIndex);
                    amsr2VCnZoneSeasonCell{k} = amsr2VCnZoneYearCell{k}(:,:,seasonIndex);
                end
                amsr2QdCnZoneSeasonCell = {...
                    (amsr2VCnZoneSeasonCell{4} - amsr2VCnZoneSeasonCell{3}) .^ 2, ...
                    (amsr2VCnZoneSeasonCell{4} - amsr2VCnZoneSeasonCell{2}) .^ 2};

                % 按照回归函数的格式要求, 将参与回归的影像矩阵重新排列.
                valueIndex = find(~isnan(modisLstMaskZoneSeasonArray));
                modisLstMaskZoneSeasonVector = modisLstMaskZoneSeasonArray(valueIndex);
                [amsr2HMaskZoneSeasonVector, amsr2VMaskZoneSeasonVector] = ...
                    deal(cell(1, channelListN));
                for k = 1 : channelListN
                    amsr2HMaskZoneSeasonVector{k} = amsr2HMaskZoneSeasonCell{k}(valueIndex);
                    amsr2VMaskZoneSeasonVector{k} = amsr2VMaskZoneSeasonCell{k}(valueIndex);
                end
                amsr2QdMaskZoneSeasonVector = [ ...
                    (amsr2VMaskZoneSeasonVector{4} - amsr2VMaskZoneSeasonVector{3}) .^ 2, ...
                    (amsr2VMaskZoneSeasonVector{4} - amsr2VMaskZoneSeasonVector{2}) .^ 2];

                amsr2BtMaskZoneSeasonRecords = double(cell2mat([amsr2HMaskZoneSeasonVector, ...
                    amsr2VMaskZoneSeasonVector, amsr2QdMaskZoneSeasonVector]));
                clear amsr2HMaskZoneSeasonCell amsr2VMaskZoneSeasonCell amsr2QdMaskZoneSeasonVector;
                clear amsr2HMaskZoneSeasonVector amsr2VMaskZoneSeasonVector;

                % 样本数大于2时才能执行逐步回归, 否则该分区不能反演AMSR2 LST.
                % 系数 [10H, 18H, 23H, 36H, 06H, 07H, 89H, 10V, 18V, 23V, 36V, 06V, 07V, 89V, 常数].
                modisLstMaskZoneSeasonVectorN = length(modisLstMaskZoneSeasonVector);
                if modisLstMaskZoneSeasonVectorN >= 2
                    timestamp = sprintf('%s %s %s', yearStr, daynight, seasonName);
                    fprintf('分区%d 季节尺度 回归: %s.\n', zonesLcId, timestamp)
                    mdl = stepwiselm(amsr2BtMaskZoneSeasonRecords, modisLstMaskZoneSeasonVector,...
                        Lower='constant', Upper='linear', Criterion='aic');
                    lstRMSE = mdl.RMSE;
                    lstR2 = mdl.Rsquared.Ordinary;
                    amsr2LstMaskZoneSeasonVector = mdl.Fitted;
                    variableIndex = find(mdl.VariableInfo.InModel == 1);
                    pYear = mdl.Coefficients.Estimate;
                    pYear = pYear([2:end, 1]);
                    coefficientSeasonArray(j, [variableIndex; end], n) = pYear;

                    amsr2LstMaskZoneSeasonVector2 = zeros(amsr2RowN * amsr2ColN * seasonIndexN, 1);
                    amsr2LstMaskZoneSeasonVector2(valueIndex) = amsr2LstMaskZoneSeasonVector;
                    amsr2LstMaskZoneSeasonArray = reshape(amsr2LstMaskZoneSeasonVector2, ...
                        amsr2RowN, amsr2ColN, seasonIndexN);
                    clear amsr2LstMaskZoneSeasonVector2

                    % 根据回归系数计算分区内的AMSR2 LST(初始化矩阵为常数项截距).
                    amsr2LstCnZoneSeasonArray = repmat(pYear(end),amsr2RowN,amsr2ColN,seasonIndexN);
                    for k = 1 : channelListN
                        amsr2LstCnZoneSeasonArray = amsr2LstCnZoneSeasonArray + ...
                            amsr2HCnZoneSeasonCell{k} * coefficientSeasonArray(j,k,n) + ...
                            amsr2VCnZoneSeasonCell{k} * coefficientSeasonArray(j,k+channelListN,n);
                    end
                    amsr2LstCnZoneSeasonArray = amsr2LstCnZoneSeasonArray + ...
                        amsr2QdCnZoneSeasonCell{:, 1} * coefficientSeasonArray(j, k*2+1, n) + ...
                        amsr2QdCnZoneSeasonCell{:, 2} * coefficientSeasonArray(j, k*2+2, n);
                    % nan值改为0, 便于算完所有分区AMSR2 LST后将结果相加.
                    amsr2LstCnZoneSeasonArray(isnan(amsr2LstCnZoneSeasonArray)) = 0;

                    % 输出AMSR2 BT(LST)和MODIS LST的季节散点图.
                    zoneScatterSeasonName = sprintf('LstScatter_%s_%s.tif', zoneName, timestamp);
                    zoneScatterSeasonPath = fullfile(lstScatterYearPath, zoneScatterSeasonName);
                    if ~exist(zoneScatterSeasonPath, 'file')
                        f = lstScatter(amsr2LstMaskZoneSeasonVector,modisLstMaskZoneSeasonVector,...
                            timestamp, zoneName, [lstR2, lstRMSE]);
                        exportgraphics(f, zoneScatterSeasonPath);
                        close all;
                    end
                else
                    lstRMSE = nan; lstR2 = nan;
                    amsr2LstMaskZoneSeasonVector = zeros(modisLstMaskZoneSeasonVectorN, 1) * nan;
                    amsr2LstMaskZoneSeasonArray = zeros(amsr2RowN, amsr2ColN, seasonIndexN);
                    amsr2LstCnZoneSeasonArray = zeros(amsr2RowN, amsr2ColN, seasonIndexN);
                end
                rmseSeasonArray(j, n) = lstRMSE;
                nSeasonArray(j, n) = modisLstMaskZoneSeasonVectorN;
                r2SeasonArray(j, n) = lstR2;

                % 将反演的当前季节AMSR2 LST保存到年度矩阵中.
                zonesLcSeasonIndexArray = zonesLcIndexArray(:, :, seasonIndex);
                zonesLcSeasonIndexArray2 = false(amsr2RowN, amsr2ColN, validDateListN);
                zonesLcSeasonIndexArray2(:, :, seasonIndex) = zonesLcSeasonIndexArray;
                amsr2LstMaskYearArray2(zonesLcSeasonIndexArray2) = ...
                    amsr2LstMaskZoneSeasonArray(zonesLcSeasonIndexArray);
                amsr2LstCnYearArray2(zonesLcSeasonIndexArray2) = ...
                    amsr2LstCnZoneSeasonArray(zonesLcSeasonIndexArray);
                clear zonesLcSeasonIndexArray zonesLcSeasonIndexArray2
            end

            % -------------------------------------------------------------------------------------
            % 按月尺度回归AMSR2 BT和MODIS LST, 并获取反演的AMSR2 LST影像.
            for n = 1 : monthNameListN
                monthIndex = (validYearMonthList == n);
                monthIndexN = sum(monthIndex);
                monthName = monthNameList{n};

                % 筛选当前月份的数据矩阵.
                modisLstMaskZoneMonthArray = modisLstMaskZoneYearArray(:, :, monthIndex);
                [amsr2HMaskZoneMonthCell, amsr2VMaskZoneMonthCell] = deal(cell(1, channelListN));
                [amsr2HCnZoneMonthCell, amsr2VCnZoneMonthCell] = deal(cell(1, channelListN));
                for k = 1 : channelListN
                    amsr2HMaskZoneMonthCell{k} = amsr2HMaskZoneYearCell{k}(:,:,monthIndex);
                    amsr2VMaskZoneMonthCell{k} = amsr2VMaskZoneYearCell{k}(:,:,monthIndex);
                    amsr2HCnZoneMonthCell{k} = amsr2HCnZoneYearCell{k}(:,:,monthIndex);
                    amsr2VCnZoneMonthCell{k} = amsr2VCnZoneYearCell{k}(:,:,monthIndex);
                end
                amsr2QdCnZoneMonthCell = {...
                    (amsr2VCnZoneMonthCell{4} - amsr2VCnZoneMonthCell{3}) .^ 2, ...
                    (amsr2VCnZoneMonthCell{4} - amsr2VCnZoneMonthCell{2}) .^ 2};

                % 按照回归函数的格式要求, 将参与回归的影像矩阵重新排列.
                valueIndex = find(~isnan(modisLstMaskZoneMonthArray));
                modisLstMaskZoneMonthVector = modisLstMaskZoneMonthArray(valueIndex);
                [amsr2HMaskZoneMonthVector, amsr2VMaskZoneMonthVector] = deal(cell(1,channelListN));
                for k = 1 : channelListN
                    amsr2HMaskZoneMonthVector{k} = amsr2HMaskZoneMonthCell{k}(valueIndex);
                    amsr2VMaskZoneMonthVector{k} = amsr2VMaskZoneMonthCell{k}(valueIndex);
                end
                amsr2QdMaskZoneMonthVector = [ ...
                    (amsr2VMaskZoneMonthVector{4} - amsr2VMaskZoneMonthVector{3}) .^ 2, ...
                    (amsr2VMaskZoneMonthVector{4} - amsr2VMaskZoneMonthVector{2}) .^ 2];

                amsr2BtMaskZoneMonthRecords = double(cell2mat([amsr2HMaskZoneMonthVector, ...
                    amsr2VMaskZoneMonthVector, amsr2QdMaskZoneMonthVector]));
                clear amsr2HMaskZoneMonthCell amsr2VMaskZoneMonthCell amsr2QdMaskZoneMonthVector;
                clear amsr2HMaskZoneMonthVector amsr2VMaskZoneMonthVector;

                % 样本数大于2时才能执行逐步回归, 否则该分区不能反演AMSR2 LST.
                % 系数 [10H, 18H, 23H, 36H, 06H, 07H, 89H, 10V, 18V, 23V, 36V, 06V, 07V, 89V, 常数].
                modisLstMaskZoneMonthVectorN = length(modisLstMaskZoneMonthVector);
                if modisLstMaskZoneMonthVectorN >= 2
                    timestamp = sprintf('%s %s %s', yearStr, daynight, monthName);
                    fprintf('分区%d 月尺度 回归: %s.\n', zonesLcId, timestamp)
                    mdl = stepwiselm(amsr2BtMaskZoneMonthRecords, modisLstMaskZoneMonthVector, ...
                        Lower='constant', Upper='linear', Criterion='aic');
                    lstRMSE = mdl.RMSE;
                    lstR2 = mdl.Rsquared.Ordinary;
                    amsr2LstMaskZoneMonthVector = mdl.Fitted;
                    variableIndex = find(mdl.VariableInfo.InModel == 1);
                    pYear = mdl.Coefficients.Estimate;
                    pYear = pYear([2:end, 1]);
                    coefficientMonthArray(j, [variableIndex; end], n) = pYear;

                    amsr2LstMaskZoneMonthVector2 = zeros(amsr2RowN * amsr2ColN * monthIndexN, 1);
                    amsr2LstMaskZoneMonthVector2(valueIndex) = amsr2LstMaskZoneMonthVector;
                    amsr2LstMaskZoneMonthArray = reshape(amsr2LstMaskZoneMonthVector2, ...
                        amsr2RowN, amsr2ColN, monthIndexN);
                    clear amsr2LstMaskZoneMonthVector2

                    % 根据回归系数计算分区内的AMSR2 LST(初始化矩阵为常数项截距).
                    amsr2LstCnZoneMonthArray = repmat(pYear(end),amsr2RowN,amsr2ColN,monthIndexN);
                    for k = 1 : channelListN
                        amsr2LstCnZoneMonthArray = amsr2LstCnZoneMonthArray + ...
                            amsr2HCnZoneMonthCell{k} * coefficientMonthArray(j, k, n) + ...
                            amsr2VCnZoneMonthCell{k} * coefficientMonthArray(j ,k+channelListN, n);
                    end
                    amsr2LstCnZoneMonthArray = amsr2LstCnZoneMonthArray + ...
                        amsr2QdCnZoneMonthCell{:, 1} * coefficientMonthArray(j, k*2+1, n) + ...
                        amsr2QdCnZoneMonthCell{:, 2} * coefficientMonthArray(j, k*2+2, n);
                    % nan值改为0, 便于算完所有分区AMSR2 LST后将结果相加.
                    amsr2LstCnZoneMonthArray(isnan(amsr2LstCnZoneMonthArray)) = 0;

                    % 输出AMSR2 BT(LST)和MODIS LST的月份散点图.
                    zoneScatterMonthName = sprintf('LstScatter_%s_%s.tif', zoneName, timestamp);
                    zoneScatterMonthPath = fullfile(lstScatterYearPath, zoneScatterMonthName);
                    if ~exist(zoneScatterMonthPath, 'file')
                        f = lstScatter(amsr2LstMaskZoneMonthVector, modisLstMaskZoneMonthVector,...
                            timestamp, zoneName, [lstR2, lstRMSE]);
                        exportgraphics(f, zoneScatterMonthPath);
                        close all;
                    end
                else
                    lstRMSE = nan; lstR2 = nan;
                    amsr2LstMaskZoneMonthVector = zeros(modisLstMaskZoneMonthVectorN, 1) * nan;
                    amsr2LstMaskZoneMonthArray = zeros(amsr2RowN, amsr2ColN, monthIndexN);
                    amsr2LstCnZoneMonthArray = zeros(amsr2RowN, amsr2ColN, monthIndexN);
                end
                rmseMonthArray(j, n) = lstRMSE;
                nMonthArray(j, n) = modisLstMaskZoneMonthVectorN;
                r2MonthArray(j, n) = lstR2;

                % 将反演的当前月份AMSR2 LST保存到年度矩阵中.
                zonesLcMonthIndexArray = zonesLcIndexArray(:, :, monthIndex);
                zonesLcMonthIndexArray2 = false(amsr2RowN, amsr2ColN, validDateListN);
                zonesLcMonthIndexArray2(:, :, monthIndex) = zonesLcMonthIndexArray;
                amsr2LstMaskYearArray3(zonesLcMonthIndexArray2) = ...
                    amsr2LstMaskZoneMonthArray(zonesLcMonthIndexArray);
                amsr2LstCnYearArray3(zonesLcMonthIndexArray2) = ...
                    amsr2LstCnZoneMonthArray(zonesLcMonthIndexArray);
                clear zonesLcMonthIndexArray zonesLcMonthIndexArray2
            end

            % 输出分区AMSR2 BT和MODIS LST影像(原AMSRE_LST_Retrieval_ZQ_16.m程序里有此步骤, 考虑去除).

        end
        save(regressionPureMatPath, 'rmseYearVector', 'nYearVector', 'r2YearVector', ...
            'rmseSeasonArray', 'nSeasonArray', 'r2SeasonArray', ...
            'rmseMonthArray', 'nMonthArray', 'r2MonthArray', ...
            'coefficientYearArray', 'coefficientSeasonArray', 'coefficientMonthArray', ...
            'amsr2LstMaskYearArray1', 'amsr2LstCnYearArray1', 'amsr2LstMaskYearArray2', ...
            'amsr2LstCnYearArray2', 'amsr2LstMaskYearArray3', 'amsr2LstCnYearArray3');
    end

%     system('shutdown -s -t 60')
    
    if 1 == 1
        continue
    end
    
    % ==============================================================================================
    % 模型优化.
    % 模型优化基于月尺度模型. 当样本数太少导致月尺度模型无法回归或精度较差时, 使用所在季节或年份的模型替换月
    % 尺度模型. !!!!! 此优化策略在某些分区会导致反演的AMSR2 LST存在较大误差, 需要改进. !!!!!

    load(regressionPureMatPath, 'rmseYearVector', 'nYearVector', 'r2YearVector', ...
        'rmse*Array', 'r2*Array', 'n*Array', 'coefficient*Array');

    % 分区范围内的平均像元数. 由于每日的积雪变化, 普通分区(1-71)的面积每日都在变化, 所以使用分区的平均像元个
    %   数来控制模型.
    zonesLcPixelCountVector = zeros(zonesLcIdListN, 1);
    for j = 1 : zonesLcIdListN
        zoneLcIdConut = sum(fixedZonesLcArray == zonesLcIdList(j), 'all');
        zonesLcPixelCountVector(j) = zoneLcIdConut / validDateListN;
    end

    % 将像元样本数少于1倍分区平均像元数的月份模型系数替换为所在季节的模型系数. 如果季节模型的像元样本数小于3
    %   倍的分区平均像元数, 则使用年尺度模型系数替换. 如果一年内的像素数仍然很少, 则借用其他相似土地类型分区
    %   的系数.
    fixedCoefficientMonthArray = coefficientMonthArray;
    fixedRmseMonthArray = rmseMonthArray;
    fixedR2MonthArray = r2MonthArray;
    fixedNMonthArray = nMonthArray;
    for j = 1 : pureLcIdN
        for k = 1 : monthNameListN
            monthSampleRatio = nMonthArray(j, k) / zonesLcPixelCountVector(j);
            if monthSampleRatio < 1
                for n = 1 : seasonMonthsListN
                    if ismember(k, seasonMonthsList{n}) % 用于确定当前月所属的季节.
                        break
                    end
                end
                seasonSampleRatio = nSeasonArray(j, n) / zonesLcPixelCountVector(j);
                if seasonSampleRatio < 3
                    fixedCoefficientMonthArray(j, :, k) = coefficientYearArray(j, :);
                    fixedRmseMonthArray(j, k) = rmseYearVector(j);
                    fixedR2MonthArray(j, k) = r2YearVector(j);
                    fixedNMonthArray(j, k) = nYearVector(j);
                else
                    fixedCoefficientMonthArray(j, :, k) = coefficientSeasonArray(j, :, n);
                    fixedRmseMonthArray(j, k) = rmseSeasonArray(j, n);
                    fixedR2MonthArray(j, k) = r2SeasonArray(j, n);
                    fixedNMonthArray(j, k) = nSeasonArray(j, n);
                end
            end
        end
    end

    % 处理一整年都没有足够像元样本的分区.
    % 建筑: 5001-5005,5007-5008(zonesLcIdList: 91-97), 水体: 3001-3008(zonesLcIdList: 81-88).
    % 2012年: 96, 97: nan. 


    % ==============================================================================================
    % 使用优化后的模型反演纯像元分区(1-71; 1000, 1100, ..., 2000; 3000, 4000, ..., 6000)的AMSR2 LST.
    amsr2LstCnPureLcYearArray = zeros(amsr2RowN, amsr2ColN, validDateListN);
    for j = 1 : pureLcIdN
        zonesLcId = zonesLcIdList(j);
        zoneName = sprintf('Zone %d', zonesLcId);
        fprintf('反演Zone %d的AMSR2地表温度.\n', zonesLcId);

        % 从年度矩阵中获取分区AMSR2 BT和MODIS LST影像.
        zonesLcIndexArray = (fixedZonesLcArray == zonesLcIdList(j));
        [amsr2HCnZoneYearCell, amsr2VCnZoneYearCell] = deal(cell(channelListN, 1));
        for k = 1 : channelListN
            tempArray = eval(sprintf('amsr2H%sYearArray', channelList{k}));
            amsr2HCnZoneYearCell{k} = setnan(tempArray, ~zonesLcIndexArray);
            tempArray = eval(sprintf('amsr2V%sYearArray', channelList{k}));
            amsr2VCnZoneYearCell{k} = setnan(tempArray, ~zonesLcIndexArray);
        end

        % 反演纯像元分区每个月的AMSR2 LST.
        for k = 1 : monthNameListN
            monthIndex = (validYearMonthList == k);
            monthIndexN = sum(monthIndex);

            % 提取当前分区各通道的AMSR2 BT.
            [amsr2HCnZoneMonthCell, amsr2VCnZoneMonthCell] = deal(cell(channelListN, 1));
            for n = 1 : channelListN
                amsr2HCnZoneMonthCell{n} = amsr2HCnZoneYearCell{n}(:, :, monthIndex);
                amsr2VCnZoneMonthCell{n} = amsr2VCnZoneYearCell{n}(:, :, monthIndex);
            end

            % 反演AMSR2 LST.
            pFixed = fixedCoefficientMonthArray(j, :, k);
            amsr2LstCnZoneMonthArray = ones(amsr2RowN, amsr2ColN, monthIndexN) * pFixed(end);
            for n = 1 : channelListN
                amsr2LstCnZoneMonthArray = amsr2LstCnZoneMonthArray + ...
                    amsr2HCnZoneMonthCell{n} * pFixed(n) + ...
                    amsr2VCnZoneMonthCell{n} * pFixed(n + channelListN);
            end
            amsr2LstCnZoneMonthArray(isnan(amsr2LstCnZoneMonthArray)) = 0;

            % 将反演的当前纯像元分区和月份的AMSR2 LST保存到年度矩阵中.
            zonesLcMonthIndexArray = zonesLcIndexArray(:, :, monthIndex);
            zonesLcMonthIndexArray2 = false(amsr2RowN, amsr2ColN, validDateListN);
            zonesLcMonthIndexArray2(:, :, monthIndex) = zonesLcMonthIndexArray;
            amsr2LstCnPureLcYearArray(zonesLcMonthIndexArray2) = ...
                amsr2LstCnZoneMonthArray(zonesLcMonthIndexArray);
        end
    end
    % ----------------------------------------------------------------------------------------------
    % 反演混合像元的AMSR2 LST.
    % 在第一次循环中, 组合包含水体, 冰川, 建筑, 积雪和其他地表覆盖类型的不完全LST. 在第二次循环中, 补上剩下
    %   的沙漠LST部分. 沙漠的比例模型与其他地表覆盖类型不同, 因此他们不能在同一个循环中处理.

    % 循环1: 遍历每个普通分区(1-71), 并寻找其中的混合像元, 然后通过对水域, 冰川, 建筑, 积雪和其他地表覆盖类
    %   型的加权平均计算不完全组分温度.
    zonesIdList = zonesLcIdList(zonesLcIdList < 1000);
    zonesIdListN = length(zonesIdList);
    amsr2LstCnMixPixelYearArray = zeros(amsr2RowN, amsr2ColN, validDateListN);
    for j = 1 : zonesIdListN
        % 寻找年度分区矩阵内每个分区的混合像元. 例: 1 + 10000 表示分区1中的混合像元.
        zoneId = zonesIdList(j);
        mixLcId = zoneId + lcCodeList(end);
        mixPixelIndexArray = (fixedZonesLcArray == mixLcId);
        fprintf('反演混合像元%d的AMSR2 LST', mixLcId);
        % 跳过没有混合像元的分区.
        if sum(mixPixelIndexArray, 'all') == 0
            continue
        end
        % 获取当前分区全年所有通道混合像元的AMSR2 BT影像.
        [amsr2HCnMixPixelYearArray, amsr2VCnMixPixelYearArray] = deal(cell(1, channelListN));
        for k = 1 : channelListN
            tempArray = eval(sprintf('amsr2H%sYearArray', channelList{k}));
            amsr2HCnMixPixelYearArray{k} = setnan(tempArray, ~mixPixelIndexArray);
            tempArray = eval(sprintf('amsr2V%sYearArray', channelList{k}));
            amsr2VCnMixPixelYearArray{k} = setnan(tempArray, ~mixPixelIndexArray);
        end
        % 获取全年当前分区中每个地表覆盖类型(水域, 冰川, 建筑, 积雪, 其他)的面积比例.
        otherPctMixPixelYearArray = setnan(otherPctArray, ~mixPixelIndexArray);
        waterPctMixPixelYearArray = setnan(waterPctArray, ~mixPixelIndexArray);
        glacierPctMixPixelYearArray = setnan(glacierPctArray, ~mixPixelIndexArray);
        buildingPctMixPixelYearArray = setnan(buildingPctArray, ~mixPixelIndexArray);
        snowPctMixPixelYearArray = setnan(snowPctArray, ~mixPixelIndexArray);

        % 使用优化后的月尺度模型计算混合像元的AMSR2 LST.
        for k = 1 : monthNameListN
            monthIndex = (validYearMonthList == k);
            monthIndexN = sum(monthIndex);
            mixPixelMonthIndexArray = mixPixelIndexArray(:, :, monthIndex);
            % 跳过没有混合像元的月份.
            if sum(mixPixelMonthIndexArray, 'all') == 0
                continue
            end
            % 获取当前月份的AMSR2 BT影像.
            [amsr2HCnMixPixelMonthArray, amsr2VCnMixPixelMonthArray] = deal(cell(1, channelListN));
            for m = 1 : channelListN
                amsr2HCnMixPixelMonthArray{m} = amsr2HCnMixPixelYearArray{m}(:, :, monthIndex);
                amsr2VCnMixPixelMonthArray{m} = amsr2VCnMixPixelYearArray{m}(:, :, monthIndex);
            end
            % 获取当前月份和分区中每个地表覆盖类型(水域, 冰川, 建筑, 积雪, 其他)的面积比例.
            otherPctMixPixelMonthArray = otherPctMixPixelYearArray(:, :, monthIndex);
            waterPctMixPixelMonthArray = waterPctMixPixelYearArray(:, :, monthIndex);
            glacierPctMixPixelMonthArray = glacierPctMixPixelYearArray(:, :, monthIndex);
            buildingPctMixPixelMonthArray = buildingPctMixPixelYearArray(:, :, monthIndex);
            snowPctMixPixelMonthArray = snowPctMixPixelYearArray(:, :, monthIndex);

            % 寻找特定地表覆盖类型区AMSR2 LST的反演参数.
            % regionNodes: [1, 5, 9, 14, 18, 28, 31, 57, 71+1]
            % regionIDs: [1, 2, 3, 4, 5, 6, 7, 8]
            % 将zoneId转换为regionId, 并寻找地表覆盖在zonesLcIdList中的索引.
            for m = 1 : regionsN
                if (regionNodes(m) <= zoneId) && (zoneId < regionNodes(m+1))
                    regionId = m;
                    break
                end
            end

            % 如果该分区包含特定的地表覆盖, 获取其回归系数, 否则将其系数设为0.
            waterLcIdIndex = (zonesLcIdList == (lcCodeList(2) + regionId));
            glacierLcIdIndex = (zonesLcIdList == (lcCodeList(3) + regionId));
            buildingLcIdIndex = (zonesLcIdList == (lcCodeList(4) + regionId));
            snowLcIdIndex = (zonesLcIdList == (lcCodeList(5) + regionId));

            pOtherMonth = fixedCoefficientMonthArray(j, :, k);
            [pWaterMonth, pGlacierMonth, pSnowMonth, pBuildingMonth] = deal(zeros(1, variablesN));
            if sum(waterLcIdIndex) == 1
                pWaterMonth = fixedCoefficientMonthArray(waterLcIdIndex, :, i);
            end
            if sum(glacierLcIdIndex) == 1
                pGlacierMonth = fixedCoefficientMonthArray(glacierLcIdIndex, :, i);
            end
            if sum(snowLcIdIndex) == 1
                pSnowMonth = fixedCoefficientMonthArray(snowLcIdIndex, :, i);
            end
            if sum(buildingLcIdIndex) == 1
                pBuildingMonth = fixedCoefficientMonthArray(buildingLcIdIndex, :, i);
            end
            % 存在没有某个地表覆盖类型像元的区域. 需要将其他分区的系数赋给该分区. 举例:
            % pBuildingMonth = fixedCoefficientMonthArray(zonesLcIdList == 10, :, i);

            % 回归.
            coefficientsAmsr2Array = zeros(amsr2RowN, amsr2ColN, monthIndexN, variablesN);
            for m = 1 : variablesN
                coefficientsAmsr2Array(:, :, :, m) = pOtherMonth(m).*otherPctMixPixelMonthArray+...
                    pWaterMonth(m) .* waterPctMixPixelMonthArray + ...
                    pGlacierMonth(m) .* glacierPctMixPixelMonthArray + ...
                    pBuildingMonth(m) .* buildingPctMixPixelMonthArray + ...
                    pSnowMonth(m) .* snowPctMixPixelMonthArray;
            end
            amsr2PartialLstCnMixedLcMonthArray = coefficientsAmsr2Array(:, :, :, end);
            for m = 1 : channelListN
                amsr2PartialLstCnMixedLcMonthArray = amsr2PartialLstCnMixedLcMonthArray + ...
                    coefficientsAmsr2Array(:, :, :, m) .* amsr2HCnMixPixelMonthArray{m} + ...
                    coefficientsAmsr2Array(:, :, :, m+7) .* amsr2VCnMixPixelMonthArray{m};
            end
            amsr2PartialLstCnMixedLcMonthArray(isnan(amsr2PartialLstCnMixedLcMonthArray)) = 0;

            mixPixelMonthIndexArray2 = false(amsr2RowN, amsr2ColN, validDateListN);
            mixPixelMonthIndexArray2(:, :, monthIndex) = mixPixelMonthIndexArray;
            amsr2LstCnMixPixelYearArray(mixPixelMonthIndexArray2) = ...
                amsr2PartialLstCnMixedLcMonthArray(mixPixelMonthIndexArray);
        end
    end

    % 循环2. 将剩余的沙漠部分LST添加到循环1获取的不完全LST中. 遍历每个沙漠分区(1000, 1100, ..., 2000).
    desertCodeArray = repmat(desertCodeLayer, [1 1 validDateListN]);
    mixedDesertIndexArray = (0 < desertPctArray) & (desertPctArray <= lcThreshold);
    mixedLcsIndexArray = (fixedZonesLcArray > lcCodeList(end));
    for j = 1 : desertCodeListN
        % 确定每个沙漠区的混合像元.
        desertId = desertCodeList(n);
        fprintf('反演沙漠区%d中混合像元的AMSR2 LST.', desertId);

        desertCodeIndexArray = (desertCodeArray == desertId);
        desertMixedPixelIndexArray = mixedDesertIndexArray & mixedLcsIndexArray & ...
            desertCodeIndexArray;

        mixedDesertPctYearArray = setnan(desertPctArray, ~desertMixedPixelIndexArray);
        [amsr2HCnMixedDesertYearCell, amsr2VCnMixedDesertYearCell] = deal(cell(channelListN, 1));
        for k = 1 : channelListN
            tempArray = eval(sprintf('amsr2H%sYearArray', channelList{k}));
            amsr2HCnMixedDesertYearCell{k} = setnan(tempArray, ~desertMixedPixelIndexArray);
            tempArray = eval(sprintf('amsr2V%sYearArray', channelList{k}));
            amsr2VCnMixedDesertYearCell{k} = setnan(tempArray, ~desertMixedPixelIndexArray);
        end

        for k = 1 : monthNameListN
            monthIndex = (validYearMonthList == k);
            monthIndexN = sum(monthIndex);
            desertMixedPixelMonthIndexArray = desertMixedPixelIndexArray(:, :, monthIndex);
            % 跳过没有混合像元沙漠分区的月份.
            if sum(desertMixedPixelMonthIndexArray, 'all') == 0
                continue
            end

            mixedDesertPctMonthArray = mixedDesertPctYearArray(:, :, monthIndex);
            [amsr2HCnMixedDesertMonthCell, amsr2VCnMixedDesertMonthCell] = deal(cell(channelListN, 1));
            for m = 1 : channelListN
                amsr2HCnMixedDesertMonthCell{m} = amsr2HCnMixedDesertYearCell{m}(:, :, monthIndex);
                amsr2VCnMixedDesertMonthCell{m} = amsr2VCnMixedDesertYearCell{m}(:, :, monthIndex);
            end

            % 确定沙漠区反演AMSR2 LST的回归系数.
            pDesertMonth = zeros(1, variablesN);
            desertLcIdIndex = find(zonesLcIdList == desertId);
            if ~isempty(desertLcIdIndex)
                pDesertMonth = fixedCoefficientMonthArray(desertLcIdIndex, :, i);
            end

            % 回归.
            amsr2LstCnPart2MonthArray = repmat(pDesertMonth(end), [amsr2RowN amsr2ColN monthIndexN]);
            for m = 1 : channelListN
                amsr2LstCnPart2MonthArray = amsr2LstCnPart2MonthArray + ...
                    pDesertMonth(m) .* amsr2HCnMixedDesertMonthCell{m} + ...
                    pDesertMonth(m+channelListN) .* amsr2VCnMixedDesertMonthCell{m};
            end
            amsr2LstCnPart2MonthArray = amsr2LstCnPart2MonthArray .* mixedDesertPctMonthArray;
            amsr2LstCnPart2MonthArray(isnan(amsr2LstCnPart2MonthArray)) = 0;

            desertMixPixelMonthIndexArray2 = false(amsr2RowN, amsr2ColN, validDateListN);
            desertMixPixelMonthIndexArray2(:, :, monthIndex) = desertMixedPixelMonthIndexArray;
            amsr2LstCnMixPixelYearArray(desertMixPixelMonthIndexArray2) = ...
                amsr2LstCnMixPixelYearArray(desertMixPixelMonthIndexArray2) + ...
                amsr2LstCnPart2MonthArray(desertMixedPixelMonthIndexArray);
        end
    end

    % 合并纯像元和混合像元.
    amsr2LstCnYearArray1 = amsr2LstCnPureLcYearArray + amsr2LstCnMixPixelYearArray;

    % 去掉AMSR2 LST中的异常值.

    % ==============================================================================================
    % 修复AMSR2 LST影像中轨道间隙的数据空缺.
    % 插值AMSR2 LST轨道间隙的日期间隔.
    halfInterval = 4;
    interval = halfInterval * 2 + 1;

    % 给AMSR2 LST矩阵添加首尾图层, 便于时空插值.
    amsr2LstCnHeadArray = amsr2LstCnYearArray1(:, :, end-halfInterval+1 : end);
    amsr2LstCnTailArray = amsr2LstCnYearArray1(:, :, 1 : halfInterval);
    amsr2LstCnYearExtArray = cat(3, amsr2LstCnHeadArray, amsr2LstCnYearArray1, amsr2LstCnTailArray);

    % 使用像元前后几天的时间序列的插值来填补每日AMSR2 LST影像的轨道间隙空缺.
    for j = 1 : validDateListN
        dateStr = validDateList{j};
        fprintf('填补%s的AMSR2 LST轨道间隙.\n', dateStr);

        % 准备用于插值的数据.
        fixedZonesLcLayer = fixedZonesLcArray(:, :, 1);
        amsr2LstCnDailyLayer = amsr2LstCnYearArray1(:, :, j);
        amsr2LstCnIntervalArray = amsr2LstCnYearExtArray(:, :, j : j+halfInterval*2);
        amsr2LstCnDailyGapIndexList = find(isnan(amsr2LstCnDailyLayer) & ~isnan(fixedZonesLcLayer));
        amsr2LstCnDailyGapIndexListN = length(amsr2LstCnDailyGapIndexList);

        % 插值.
        for k = 1 : amsr2LstCnDailyGapIndexListN
            gapIndex = amsr2LstCnDailyGapIndexList(k);
            yVector = zeros(interval, 1) * nan;
            for m = 1 : interval
                amsr2LstCnDailyExtLayer = amsr2LstCnIntervalArray(:, :, m);
                yVector(m) = amsr2LstCnDailyExtLayer(gapIndex);
            end
            xVector = 1 : interval;
            xVector(isnan(yVector)) = [];
            yVector(isnan(yVector)) = [];
            if isempty(xVector)
                amsr2LstGapValue = nan;
            elseif length(xVector) == 1
                amsr2LstGapValue = yVector;
            else
                amsr2LstGapValue = interp1(xVector, yVector, (interval+1)/2);
            end
            amsr2LstCnDailyLayer(gapIndex) = amsr2LstGapValue;
        end
        % 用每日插值后的数据替换原来有轨道间隙的数据.
        amsr2LstCnYearArray1(:, :, j) = amsr2LstCnDailyLayer;
    end

    % ==============================================================================================
    % 输出修复轨道间隙空缺后的AMSR2 LST.
    for j = 1 : validDateListN
        amsr2LstCnDailyLayer = amsr2LstCnYearArray1(:, :, j);
        geotiffwrite(amsr2LstCnPath, amsr2LstCnDailyLayer, amsr2Ref);
    end

end

























%% 自定义函数.
% AMSR2和MODIS地表温度的散点图.
function f = lstScatter(amsr2LstVector, modisLstVector, zoneName, timestamp, r2Rmse)
validIndex =  ~isnan(amsr2LstVector) & ~isnan(modisLstVector);
amsr2LstVector = amsr2LstVector(validIndex);
modisLstVector = modisLstVector(validIndex);

lstBias = mean(amsr2LstVector - modisLstVector);
lstMAE = mean(abs(amsr2LstVector - modisLstVector));

f = figure;
plot(amsr2LstVector, modisLstVector, '.k', [220, 360], [220, 360], 'r');
xlabel('AMSR2 LST'); ylabel('MODIS LST');
title(sprintf('AMSR2 LST vs MODIS LST in %s %s', zoneName, timestamp));

txt1 = ['N: ', num2str(sum(~isnan(modisLstVector)))];
txt2 = ['R^2: ', num2str(r2Rmse(1), '%.3f')];
txt3 = ['Bias: ', num2str(lstBias, '%.3f')];
txt4 = ['MAE: ', num2str(lstMAE, '%.3f')];
txt5 = ['RMSE: ', num2str(r2Rmse(2), '%.3f')];
text(0.7, 0.29, txt1, 'Units', 'normalized', 'FontSize', 12);
text(0.7, 0.23, txt2, 'Units', 'normalized', 'FontSize', 12);
text(0.7, 0.17, txt3, 'Units', 'normalized', 'FontSize', 12);
text(0.7, 0.11, txt4, 'Units', 'normalized', 'FontSize', 12);
text(0.7, 0.05, txt5, 'Units', 'normalized', 'FontSize', 12);
end
