%% AMSR2地表温度反演.

%% 功能标记与预设参数.
% 指定白天和晚上的标记. 1表示白天(升轨), 2表示晚上(降轨).
flg1 = 1;

% AMSR2数据的通道, 极化.
channelList = ["10", "18", "23", "36", "06", "07", "89"];
channelListN = length(channelList);
polarize = ["H", "V"];
polarizeN = length(polarize);

% 回归模型中的变量个数(14个通道 + 2个二次项 + 1常数项).
variablesN = channelListN * polarizeN + 3;

% 数据年份列表(时间区间2012/07/02-2019/12/31).
yearList = 2012 : 2019;
yearListN = length(yearList);

% 各月份的名称.
monthNameList = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov','Dec'};
monthNameListN = length(monthNameList);

% 昼夜标记.
daynight = {'Day', 'Night'};
daynight = daynight{flg1};

% LC和分区影像代码的重分类. 主要针对水域, 冰川, 建筑这三个LC区中LST的反演.
%   [阶梯1_东北东部(1): 1-3, 阶梯1_华北(2): 4-5, 阶梯1_华南(3): 6-10]
%   [阶梯2_西南(4): 11-14, 阶梯2_西北东部(5): 15-25, 阶梯2_东北西部(6): 26-29]
%   [阶梯2_西北西部(7): 30-46， 阶梯3_青藏高原(8): 47-62]
% 重分类时, 使用了右开区间, 导致最后一个分区编码71排除在分类结果中, 62+1, 使其包含在重分类结果中.
regionNodes = [1, 4, 6, 11, 15, 26, 30, 47, 62+1]; 
regionsN = length(regionNodes) - 1;

% AMSR2 BT影像中的NoData.
amsr2Nodata = [65534, 65535];

%% 路径.
% 根目录.
rootPath = 'E:\AMSR_MODIS_Fusion';
dataPath = fullfile(rootPath, 'Data');
funcPath = fullfile(rootPath, 'Code\Functions');
figPath = fullfile(rootPath, 'Figures');
addpath(funcPath);
if ~exist(figPath, 'dir')
    mkdir(figPath)
end

% 输入数据路径.
amsr2CnMatPath = fullfile(dataPath, 'AMSR2_2_CN_Matlab');
amsr2MaskMatPath = fullfile(dataPath, 'AMSR2_3_MaskCn_Matlab');
modisLstMaskMatPath = fullfile(dataPath, 'MYD11A1_3_MaskCn_Matlab');
ccsevMatPath = fullfile(dataPath, 'CCSEV_Matlab');
regressionMatPath = fullfile(dataPath, 'Regression_Matlab'); 

% 输出反演的AMSR2 LST路径.
amsr2LstMatPath = fullfile(dataPath, 'AMSR2_4_LST_Matlab');
if ~exist(amsr2LstMatPath, 'dir')
    mkdir(amsr2LstMatPath)
end
amsr2LstPath = fullfile(dataPath, 'AMSR2_4_LstCn_TIF');
if ~exist(amsr2LstPath, 'dir')
    mkdir(amsr2LstPath)
end
amsr2LstGapFilledPath = fullfile(dataPath, 'AMSR2_4_LstGapfilled_TIF');
if ~exist(amsr2LstGapFilledPath, 'dir')
    mkdir(amsr2LstGapFilledPath)
end

%% 回归和输出.
% 分年度反演AMSR2 LST.
for i = 1 : yearListN
    yearStr = num2str(yearList(i));
%     yearStr = '2013';

    % 输出文件的路径.
    amsr2LstYearFolder = sprintf('AMSR2_LST_%sXXXX_TIF', yearStr);
    amsr2LstCnYearPath = fullfile(amsr2LstPath, amsr2LstYearFolder);
    if ~exist(amsr2LstCnYearPath, 'dir')
        mkdir(amsr2LstCnYearPath);
    end
    amsr2LstGapFilledYearPath = fullfile(amsr2LstGapFilledPath, amsr2LstYearFolder);
    if ~exist(amsr2LstGapFilledYearPath, 'dir')
        mkdir(amsr2LstGapFilledYearPath);
    end

    % 从Mat文件中读取中国区原始的AMSR2 BT数据.
    amsr2CnMatFileName = sprintf('AMSR2_BT_%s_%s.mat', yearStr, daynight);
    amsr2CnMatFilePath = fullfile(amsr2CnMatPath, amsr2CnMatFileName);
    load(amsr2CnMatFilePath, 'amsr2H*YearArray', 'amsr2V*YearArray');

    % 从Mat文件中读取Mask后的AMSR2 BT数据.
    amsr2MaskMatFileName = sprintf('AMSR2_BT_MaskCn_%s_%s.mat', yearStr, daynight);
    amsr2MaskMatFilePath = fullfile(amsr2MaskMatPath, amsr2MaskMatFileName);
    load(amsr2MaskMatFilePath, 'lstDateFilterList', 'amsr2Ref');

    % 从Mat文件中读取CCSEV.
    ccsevFilePath = fullfile(ccsevMatPath, sprintf('CCSEV_%s.mat', yearStr));
    load(ccsevFilePath, 'zonesLcCodeList', 'lcCodeList', 'lcDateList', 'fixedZonesLcArray');

    % 从Mat文件中读取纯像元分区的回归系数.
%     regressionPureMatName = sprintf('Regression_Pure_%s_%s.mat', yearStr, daynight);
    regressionPureMatName = sprintf('Regression_Pure_2013_%s.mat', daynight);
    regressionPureMatPath = fullfile(regressionMatPath, regressionPureMatName);
    load(regressionPureMatPath, 'fixedCoefficientMonthArray');

    % 获取共有日期的分区和数据矩阵.
    [validDateList, lstDateIndex, lcDateIndex] = intersect(lstDateFilterList, lcDateList);
    validYearMonthList = datetime(validDateList , 'InputFormat', 'yyyyMMdd').Month;
    validDateListN = length(validDateList);
    fixedZonesLcArray = fixedZonesLcArray(:, :, lcDateIndex);

    % 获取ASMR2 LST数据的行列数.
    [amsr2RowN, amsr2ColN] = size(fixedZonesLcArray, [1 2]);

    % ----------------------------------------------------------------------------------------------
    % 使用优化后的模型反演纯像元分区(1-62; 100, 110, ..., 200; 300, 400, ..., 600)的AMSR2 LST.
    amsr2LstPureYearMatName = sprintf('AMSR2_LstPure_%s_%s.mat', daynight, yearStr);
    amsr2LstPureYearMatPath = fullfile(amsr2LstMatPath, amsr2LstPureYearMatName);
    if ~exist(amsr2LstPureYearMatPath, 'file')
        amsr2LstCnPureLcYearArray = zeros(amsr2RowN, amsr2ColN, validDateListN, 'single');
        pureZonesLcCodeList = zonesLcCodeList(zonesLcCodeList < lcCodeList(end));  % < 1000
        for j = 1 : length(pureZonesLcCodeList)
            zonesLcCode = pureZonesLcCodeList(j);
            fprintf('反演分区%d %s年%s的AMSR2地表温度.\n', zonesLcCode, yearStr, daynight);

            % 从年度矩阵中获取分区AMSR2 BT和MODIS LST影像.
            zonesLcIndexArray = (fixedZonesLcArray == zonesLcCode);
            [amsr2HCnZoneYearCell, amsr2VCnZoneYearCell] = deal(cell(channelListN, 1));
            for k = 1 : channelListN
                amsr2HCnYearArray = single(eval(sprintf('amsr2H%sYearArray', channelList{k})));
                amsr2VCnYearArray = single(eval(sprintf('amsr2V%sYearArray', channelList{k})));
                amsr2HCnYearArray(ismember(amsr2HCnYearArray, amsr2Nodata)) = nan;
                amsr2VCnYearArray(ismember(amsr2VCnYearArray, amsr2Nodata)) = nan;
                amsr2HCnYearArray = amsr2HCnYearArray(:, :, lstDateIndex) / 100;
                amsr2VCnYearArray = amsr2VCnYearArray(:, :, lstDateIndex) / 100;
                amsr2HCnZoneYearCell{k} = setnan(amsr2HCnYearArray, ~zonesLcIndexArray);
                amsr2VCnZoneYearCell{k} = setnan(amsr2VCnYearArray, ~zonesLcIndexArray);
            end
            clear amsr2HCnYearArray amsr2VCnYearArray

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
                amsr2QdCnZoneMonthCell = {...
                    (amsr2VCnZoneMonthCell{4} - amsr2VCnZoneMonthCell{3}) .^ 2, ...
                    (amsr2VCnZoneMonthCell{4} - amsr2VCnZoneMonthCell{2}) .^ 2};

                % 反演AMSR2 LST.
                fixedP = single(fixedCoefficientMonthArray(j, :, k));
                amsr2LstCnZoneMonthArray = repmat(fixedP(end), amsr2RowN, amsr2ColN, monthIndexN);
                for n = 1 : channelListN
                    amsr2LstCnZoneMonthArray = amsr2LstCnZoneMonthArray + ...
                        amsr2HCnZoneMonthCell{n} * fixedP(n) + ...
                        amsr2VCnZoneMonthCell{n} * fixedP(n+channelListN);
                end
                amsr2LstCnZoneMonthArray = amsr2LstCnZoneMonthArray + ...
                    amsr2QdCnZoneMonthCell{:, 1} * fixedP(n*2+1) + ...
                    amsr2QdCnZoneMonthCell{:, 2} * fixedP(n*2+2);
                amsr2LstCnZoneMonthArray(isnan(amsr2LstCnZoneMonthArray)) = 0;

                % 将反演的当前纯像元分区和月份的AMSR2 LST保存到年度矩阵中.
                zonesLcMonthIndexArray = zonesLcIndexArray(:, :, monthIndex);
                zonesLcMonthIndexArray2 = false(amsr2RowN, amsr2ColN, validDateListN);
                zonesLcMonthIndexArray2(:, :, monthIndex) = zonesLcMonthIndexArray;
                amsr2LstCnPureLcYearArray(zonesLcMonthIndexArray2) = ...
                    amsr2LstCnZoneMonthArray(zonesLcMonthIndexArray);
                clear zonesLcMonthIndexArray zonesLcMonthIndexArray2 amsr2LstCnZoneMonthArray
            end
        end
        save(amsr2LstPureYearMatPath, 'amsr2LstCnPureLcYearArray', 'validDateList', 'amsr2Ref');
    end

    % ----------------------------------------------------------------------------------------------
    % 反演混合像元的AMSR2 LST.
    % 在第一次循环中, 组合包含水体, 冰川, 建筑, 积雪和其他地表覆盖类型的不完全LST. 在第二次循环中, 补上剩下
    %   的沙漠LST部分. 沙漠的比例模型与其他地表覆盖类型不同, 因此他们不能在同一个循环中处理.

    % 循环1: 遍历每个普通分区(1-62), 并寻找其中的混合像元, 然后通过对水域, 冰川, 建筑, 积雪和其他地表覆盖类
    %   型的加权平均计算不完全组分温度.
    amsr2LstCnMix1YearMatName = sprintf('AMSR2_LstMix1_%s_%s.mat', daynight, yearStr);
    amsr2LstCnMix1YearMatPath = fullfile(amsr2LstMatPath, amsr2LstCnMix1YearMatName);
    if ~exist(amsr2LstCnMix1YearMatPath, 'file')
        % 读取LC的面积比例矩阵.
        load(ccsevFilePath, 'otherPctArray', 'buildingPctArray', 'glacierPctArray', ...
            'waterPctArray', 'snowPctArray');

        amsr2LstCnMix1YearArray = zeros(amsr2RowN, amsr2ColN, validDateListN, 'single');
        mixZonesCodeList = zonesLcCodeList(zonesLcCodeList > lcCodeList(end));  % > 1000
        for j = 1 : length(mixZonesCodeList)
            % 寻找年度分区矩阵内每个分区的混合像元. 例: 1 + 1000 表示分区1中的混合像元.
            mixZonesCode = mixZonesCodeList(j);
            zoneCode = mixZonesCode - lcCodeList(end);
            fprintf('反演%s年分区%d %s的混合像元AMSR2地表温度.\n', yearStr, mixZonesCode, daynight);

            % 跳过没有混合像元的分区.
            mixedPixelIndexArray = (fixedZonesLcArray == mixZonesCode);
            if sum(mixedPixelIndexArray, 'all') == 0
                continue
            end

            % 获取当前分区全年所有通道混合像元的AMSR2 BT影像.
            [amsr2HCnMixPixelYearArray, amsr2VCnMixPixelYearArray] = deal(cell(1, channelListN));
            for k = 1 : channelListN
                amsr2HCnYearArray = single(eval(sprintf('amsr2H%sYearArray', channelList{k})));
                amsr2VCnYearArray = single(eval(sprintf('amsr2V%sYearArray', channelList{k})));
                amsr2HCnYearArray(ismember(amsr2HCnYearArray, amsr2Nodata)) = nan;
                amsr2VCnYearArray(ismember(amsr2VCnYearArray, amsr2Nodata)) = nan;
                amsr2HCnYearArray = amsr2HCnYearArray(:, :, lstDateIndex) / 100;
                amsr2VCnYearArray = amsr2VCnYearArray(:, :, lstDateIndex) / 100;
                amsr2HCnMixPixelYearArray{k} = setnan(amsr2HCnYearArray, ~mixedPixelIndexArray);
                amsr2VCnMixPixelYearArray{k} = setnan(amsr2VCnYearArray, ~mixedPixelIndexArray);
            end
            clear amsr2HCnYearArray amsr2VCnYearArray

            % 获取全年当前分区中每个地表覆盖类型(水域, 冰川, 建筑, 积雪, 其他)的面积比例.
            otherPctMixPixelYearArray = setnan(otherPctArray, ~mixedPixelIndexArray);
            waterPctMixPixelYearArray = setnan(waterPctArray, ~mixedPixelIndexArray);
            glacierPctMixPixelYearArray = setnan(glacierPctArray, ~mixedPixelIndexArray);
            buildingPctMixPixelYearArray = setnan(buildingPctArray, ~mixedPixelIndexArray);
            snowPctMixPixelYearArray = setnan(snowPctArray, ~mixedPixelIndexArray);

            % 使用优化后的月尺度模型计算混合像元的AMSR2 LST.
            for k = 1 : monthNameListN
                monthIndex = (validYearMonthList == k);
                monthIndexN = sum(monthIndex);

                % 跳过没有混合像元的月份.
                mixedPixelMonthIndexArray = mixedPixelIndexArray(:, :, monthIndex);
                if sum(mixedPixelMonthIndexArray, 'all') == 0
                    continue
                end

                % 获取当前月份的AMSR2 BT影像.
                [amsr2HCnMixPixelMonthCell,amsr2VCnMixPixelMonthCell] = deal(cell(1,channelListN));
                for m = 1 : channelListN
                    amsr2HCnMixPixelMonthCell{m} = amsr2HCnMixPixelYearArray{m}(:, :, monthIndex);
                    amsr2VCnMixPixelMonthCell{m} = amsr2VCnMixPixelYearArray{m}(:, :, monthIndex);
                end
                amsr2QdCnMixPixelMonthCell = {...
                    (amsr2VCnMixPixelMonthCell{4} - amsr2VCnMixPixelMonthCell{3}) .^ 2, ...
                    (amsr2VCnMixPixelMonthCell{4} - amsr2VCnMixPixelMonthCell{2}) .^ 2};

                % 获取当前月份和分区中每个地表覆盖类型(水域, 冰川, 建筑, 积雪, 其他)的面积比例.
                otherPctMixPixelMonthArray = otherPctMixPixelYearArray(:, :, monthIndex);
                waterPctMixPixelMonthArray = waterPctMixPixelYearArray(:, :, monthIndex);
                glacierPctMixPixelMonthArray = glacierPctMixPixelYearArray(:, :, monthIndex);
                buildingPctMixPixelMonthArray = buildingPctMixPixelYearArray(:, :, monthIndex);
                snowPctMixPixelMonthArray = snowPctMixPixelYearArray(:, :, monthIndex);

                % 寻找特定地表覆盖类型区AMSR2 LST的反演参数.
                % regionNodes: [1, 4, 6, 11, 15, 26, 30, 47, 62+1]
                % regionCodes: [1, 2, 3, 4, 5, 6, 7, 8]
                % 将zoneCode转换为regionCode, 并寻找地表覆盖在zonesLcCodeList中的索引.
                for m = 1 : regionsN
                    if (regionNodes(m) <= zoneCode) && (zoneCode < regionNodes(m+1))
                        regionCode = m;
                        break
                    end
                end

                % 如果该分区包含特定的地表覆盖, 获取其回归系数, 否则将其系数设为0.
                waterLcCodeIndex = (zonesLcCodeList == (lcCodeList(2) + regionCode));
                glacierLcCodeIndex = (zonesLcCodeList == (lcCodeList(3) + regionCode));
                buildingLcCodeIndex = (zonesLcCodeList == (lcCodeList(4) + regionCode));
                snowLcCodeIndex = (zonesLcCodeList == (lcCodeList(5) + regionCode));

                pOtherMonth = fixedCoefficientMonthArray(j, :, k);
                [pWaterMonth,pGlacierMonth,pSnowMonth,pBuildingMonth] = deal(zeros(1, variablesN));
                if sum(waterLcCodeIndex) == 1
                    pWaterMonth = fixedCoefficientMonthArray(waterLcCodeIndex, :, k);
                end
                if sum(glacierLcCodeIndex) == 1
                    pGlacierMonth = fixedCoefficientMonthArray(glacierLcCodeIndex, :, k);
                end
                if sum(snowLcCodeIndex) == 1
                    pSnowMonth = fixedCoefficientMonthArray(snowLcCodeIndex, :, k);
                end
                if sum(buildingLcCodeIndex) == 1
                    pBuildingMonth = fixedCoefficientMonthArray(buildingLcCodeIndex, :, k);
                end

                % 存在没有某个地表覆盖类型像元的区域. 需要将其他分区的系数赋给该分区. 举例:
                % pBuildingMonth = fixedCoefficientMonthArray(zonesLcCodeList == 10, :, k);

                % 回归.
                coefficientsAmsr2Array = zeros(amsr2RowN,amsr2ColN,monthIndexN,variablesN,'single');
                for m = 1 : variablesN
                    coefficientsAmsr2Array(:,:,:,m) = pOtherMonth(m).*otherPctMixPixelMonthArray+...
                        pWaterMonth(m) .* waterPctMixPixelMonthArray + ...
                        pGlacierMonth(m) .* glacierPctMixPixelMonthArray + ...
                        pBuildingMonth(m) .* buildingPctMixPixelMonthArray + ...
                        pSnowMonth(m) .* snowPctMixPixelMonthArray;
                end
                amsr2PartialLstCnMixedLcMonthArray = coefficientsAmsr2Array(:, :, :, end);
                for m = 1 : channelListN
                    amsr2PartialLstCnMixedLcMonthArray = amsr2PartialLstCnMixedLcMonthArray + ...
                        amsr2HCnMixPixelMonthCell{m}.*coefficientsAmsr2Array(:,:,:,m) + ...
                        amsr2VCnMixPixelMonthCell{m}.*coefficientsAmsr2Array(:,:,:,m+channelListN);
                end
                amsr2PartialLstCnMixedLcMonthArray = amsr2PartialLstCnMixedLcMonthArray + ...
                    amsr2QdCnMixPixelMonthCell{:, 1} .* coefficientsAmsr2Array(:, :, :, m*2+1) + ...
                    amsr2QdCnMixPixelMonthCell{:, 2} .* coefficientsAmsr2Array(:, :, :, m*2+2);
                amsr2PartialLstCnMixedLcMonthArray(isnan(amsr2PartialLstCnMixedLcMonthArray)) = 0;

                mixPixelMonthIndexArray2 = false(amsr2RowN, amsr2ColN, validDateListN);
                mixPixelMonthIndexArray2(:, :, monthIndex) = mixedPixelMonthIndexArray;
                amsr2LstCnMix1YearArray(mixPixelMonthIndexArray2) = ...
                    amsr2PartialLstCnMixedLcMonthArray(mixedPixelMonthIndexArray);
                clear mixPixelMonthIndexArray2 mixedPixelMonthIndexArray;
                clear amsr2PartialLstCnMixedLcMonthArray;
            end
        end
        save(amsr2LstCnMix1YearMatPath, 'amsr2LstCnMix1YearArray', 'validDateList', 'amsr2Ref')        
    end

    % 循环2. 将剩余的沙漠部分LST添加到循环1获取的不完全LST中. 遍历每个沙漠分区(100, 110, ..., 200).
    % 获取沙漠分区编码.
    amsr2LstCnMix2YearMatName = sprintf('AMSR2_LstMix2_%s_%s.mat', daynight, yearStr);
    amsr2LstCnMix2YearMatPath = fullfile(amsr2LstMatPath, amsr2LstCnMix2YearMatName);
    if  ~exist(amsr2LstCnMix2YearMatPath, 'file')
        % 读取沙漠分区编码数据.
        load(ccsevFilePath, 'desertCodeLayer', 'desertPctArray');
        load(amsr2LstCnMix1YearMatPath, 'amsr2LstCnMix1YearArray')

        desertCodeList = unique(desertCodeLayer);
        desertCodeList(desertCodeList == 0) = [];
        desertCodeArray = repmat(desertCodeLayer, [1 1 validDateListN]);
        desertPctArray = desertPctArray(:, :, lcDateIndex);
        mixedLcsIndexArray = (fixedZonesLcArray > lcCodeList(end));
        amsr2LstCnMix2YearArray = amsr2LstCnMix1YearArray;
        for j = 1 : length(desertCodeList)
            desertCode = desertCodeList(j);
            fprintf('反演%s年沙漠区%d %s混合像元的AMSR2地表温度.\n', yearStr, desertCode, daynight);

            % 确定每个沙漠区的混合像元.
            desertCodeMixIndexArray = (desertCodeArray == desertCode) & mixedLcsIndexArray;

            % 获取每个沙漠区的混合像元AMSR2 BT数据.
            desertPctMixYearArray = setnan(desertPctArray, ~desertCodeMixIndexArray);
            [amsr2HCnMixDesertYearCell, amsr2VCnMixDesertYearCell] = deal(cell(channelListN, 1));
            for k = 1 : channelListN
                amsr2HCnYearArray = single(eval(sprintf('amsr2H%sYearArray', channelList{k})));
                amsr2VCnYearArray = single(eval(sprintf('amsr2V%sYearArray', channelList{k})));
                amsr2HCnYearArray(ismember(amsr2HCnYearArray, amsr2Nodata)) = nan;
                amsr2VCnYearArray(ismember(amsr2VCnYearArray, amsr2Nodata)) = nan;
                amsr2HCnYearArray = amsr2HCnYearArray(:, :, lstDateIndex) / 100;
                amsr2VCnYearArray = amsr2VCnYearArray(:, :, lstDateIndex) / 100;
                amsr2HCnMixDesertYearCell{k} = setnan(amsr2HCnYearArray, ~desertCodeMixIndexArray);
                amsr2VCnMixDesertYearCell{k} = setnan(amsr2VCnYearArray, ~desertCodeMixIndexArray);
            end
            clear amsr2HCnYearArray amsr2VCnYearArray

            % 反演每个沙漠区混合像元的沙漠地表温度组分.
            for k = 1 : monthNameListN
                monthIndex = (validYearMonthList == k);
                monthIndexN = sum(monthIndex);
                desertCodeMixMonthIndexArray = desertCodeMixIndexArray(:, :, monthIndex);

                % 跳过没有混合像元沙漠分区的月份.
                if sum(desertCodeMixMonthIndexArray, 'all') == 0
                    continue
                end

                % 获取每个月份当前沙漠区的混合像元AMSR2 BT数据.
                desertPctMixMonthArray = desertPctMixYearArray(:, :, monthIndex);
                [amsr2HCnMixDesertMonthCell,amsr2VCnMixDesertMonthCell] = deal(cell(channelListN,1));
                for m = 1 : channelListN
                    amsr2HCnMixDesertMonthCell{m} = amsr2HCnMixDesertYearCell{m}(:, :, monthIndex);
                    amsr2VCnMixDesertMonthCell{m} = amsr2VCnMixDesertYearCell{m}(:, :, monthIndex);
                end
                amsr2QdCnMixDesertMonthCell = {...
                    (amsr2VCnMixDesertMonthCell{4} - amsr2VCnMixDesertMonthCell{3}) .^ 2, ...
                    (amsr2VCnMixDesertMonthCell{4} - amsr2VCnMixDesertMonthCell{2}) .^ 2};

                % 确定沙漠区反演AMSR2 LST的回归系数.
                pMonth = zeros(1, variablesN);
                desertCodeIndex = find(zonesLcCodeList == desertCode);
                if ~isempty(desertCodeIndex)
                    pMonth = fixedCoefficientMonthArray(desertCodeIndex, :, k);
                end

                % 回归.
                amsr2LstCnMix2MonthArray = repmat(pMonth(end), [amsr2RowN amsr2ColN monthIndexN]);
                for m = 1 : channelListN
                    amsr2LstCnMix2MonthArray = amsr2LstCnMix2MonthArray + ...
                        amsr2HCnMixDesertMonthCell{m} .* pMonth(m) + ...
                        amsr2VCnMixDesertMonthCell{m} .* pMonth(m+channelListN);
                end
                amsr2LstCnMix2MonthArray = desertPctMixMonthArray .* (amsr2LstCnMix2MonthArray + ...
                    amsr2QdCnMixDesertMonthCell{:, 1} .* pMonth(m*2+1) + ...
                    amsr2QdCnMixDesertMonthCell{:, 2} .* pMonth(m*2+2));
                amsr2LstCnMix2MonthArray(isnan(amsr2LstCnMix2MonthArray)) = 0;

                desertCodeMixMonthIndexArray2 = false(amsr2RowN, amsr2ColN, validDateListN);
                desertCodeMixMonthIndexArray2(:, :, monthIndex) = desertCodeMixMonthIndexArray;
                amsr2LstCnMix2YearArray(desertCodeMixMonthIndexArray2) = ...
                    amsr2LstCnMix2YearArray(desertCodeMixMonthIndexArray2) + ...
                    amsr2LstCnMix2MonthArray(desertCodeMixMonthIndexArray);
                clear desertCodeMixMonthIndexArray2 desertCodeMixMonthIndexArray
                clear amsr2LstCnMix2MonthArray
            end
        end
        save(amsr2LstCnMix2YearMatPath, 'amsr2LstCnMix2YearArray', 'validDateList', 'amsr2Ref');
    end

    % ----------------------------------------------------------------------------------------------
    % 合并AMSR2 LST纯像元和混合像元, 去掉异常值.
    amsr2LstCnYearMatName = sprintf('AMSR2_LstCn_%s_%s.mat', daynight, yearStr);
    amsr2LstCnYearMatPath = fullfile(amsr2LstMatPath, amsr2LstCnYearMatName);
    if ~exist(amsr2LstCnYearMatPath, 'file')
        load(amsr2LstPureYearMatPath, 'amsr2LstCnPureLcYearArray');
        load(amsr2LstCnMix2YearMatPath, 'amsr2LstCnMix2YearArray');
        amsr2LstCnYearArray = amsr2LstCnPureLcYearArray + amsr2LstCnMix2YearArray;

        % 去掉AMSR2 LST中的异常值.

        save(amsr2LstCnYearMatPath, 'amsr2LstCnYearArray', 'validDateList', 'amsr2Ref');
    end

    % 输出反演的AMSR2 LST为TIF格式.
    fprintf('输出反演的%s年%s的AMSR2地表温度.\n', yearStr, daynight);
    load(amsr2LstCnYearMatPath, 'amsr2LstCnYearArray')
    for j = 1 : validDateListN
        amsr2LstCnName = sprintf('AMSR2_LST_%s_%s.tif', daynight, validDateList{j});
        amsr2LstCnPath = fullfile(amsr2LstCnYearPath, amsr2LstCnName);
        if exist(amsr2LstCnPath, 'file')
            continue
        end
        geotiffwrite(amsr2LstCnPath, amsr2LstCnYearArray(:, :, j), amsr2Ref, ...
            TiffTags=struct('Compression','LZW'));
    end

    continue

    % ----------------------------------------------------------------------------------------------
    % 修复AMSR2 LST影像中轨道间隙的数据空缺, 并输出每日的TIF文件.
    amsr2LstGapFilledMatName = sprintf('AMSR2_LstGapfilled_%s_%s.mat', daynight, yearStr);
    amsr2LstGapFilledMatPath = fullfile(amsr2LstMatPath, amsr2LstGapFilledMatName);
    if ~exist(amsr2LstGapFilledMatPath, 'file')
        load(amsr2LstCnYearMatPath, 'amsr2LstCnYearArray');
        % 用于插值的时间间隔区间.
        halfInterval = 4;
        interval = halfInterval * 2 + 1;

        % 给AMSR2 LST矩阵添加首尾图层, 便于时空插值.
        amsr2LstCnHeadArray = amsr2LstCnYearArray(:, :, end-halfInterval+1 : end);
        amsr2LstCnTailArray = amsr2LstCnYearArray(:, :, 1 : halfInterval);
        amsr2LstCnYearExtArray = cat(3,amsr2LstCnHeadArray,amsr2LstCnYearArray,amsr2LstCnTailArray);

        % 使用像元前后几天的时间序列的插值来填补每日AMSR2 LST影像的轨道间隙空缺.
        amsr2LstCnGapfilledYearArray = zeros(size(amsr2LstCnYearArray), 'single');
        for j = 1 : validDateListN
            fprintf('填补%s的AMSR2 LST轨道间隙.\n', validDateList{j});

            % 准备用于插值的数据.
            fixedZonesLcLayer = fixedZonesLcArray(:, :, j);
            amsr2LstCnDailyLayer = amsr2LstCnYearArray(:, :, j);
            amsr2LstCnIntervalArray = amsr2LstCnYearExtArray(:, :, j : j+halfInterval*2);
            amsr2LstCnGapIndexList = find(amsr2LstCnDailyLayer == 0 & fixedZonesLcLayer ~= 128);
            amsr2LstCnGapIndexListN = length(amsr2LstCnGapIndexList);

            % 插值.
            for k = 1 : amsr2LstCnGapIndexListN
                gapIndex = amsr2LstCnGapIndexList(k);
                yVector = zeros(interval, 1) * nan;
                for m = 1 : interval
                    amsr2LstCnDailyExtLayer = amsr2LstCnIntervalArray(:, :, m);
                    yVector(m) = amsr2LstCnDailyExtLayer(gapIndex);
                end
                xVector = (1 : interval)';
                xVector(yVector == 0) = [];
                yVector(yVector == 0) = [];
                if isempty(xVector)
                    amsr2LstGapValue = nan;
                elseif length(xVector) == 1
                    amsr2LstGapValue = yVector;
                else
                    amsr2LstGapValue = interp1(xVector, yVector, halfInterval+1);
                end
                amsr2LstCnDailyLayer(gapIndex) = amsr2LstGapValue;
            end
            amsr2LstCnGapfilledYearArray(:, :, j) = amsr2LstCnDailyLayer;
        end
        save(amsr2LstGapFilledMatPath, 'amsr2LstCnGapfilledYearArray', 'validDateList', 'amsr2Ref');
    end

    % 输出修复轨道间隙空缺后的AMSR2 LST为TIF格式.
    fprintf('输出填补过轨道间隙的%s年%s的AMSR2地表温度.\n', yearStr, daynight);
    load(amsr2LstGapFilledMatPath, 'amsr2LstCnGapfilledYearArray')
    for j = 1 : validDateListN
        amsr2LstCnName = sprintf('AMSR2_LST_%s_%s.tif', daynight, validDateList{j});
        amsr2LstCnPath = fullfile(amsr2LstGapFilledYearPath, amsr2LstCnName);
        if exist(amsr2LstCnPath, 'file')
            continue
        end
        geotiffwrite(amsr2LstCnPath, amsr2LstCnGapfilledYearArray(:, :, j), amsr2Ref, ...
            TiffTags=struct('Compression','LZW'));
    end
end

% system('shutdown -s -t 60');


