%% AMSR2地表温度反演模型评价.

%% 功能标记与预设参数.
% 指定白天和晚上的标记. 1表示白天(升轨), 2表示晚上(降轨).
flg1 = 1;

% 数据年份列表(时间区间2012/07/02-2019/12/31).
yearList = 2012 : 2019;
yearListN = length(yearList);

% 各月份的名称.
monthList = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
monthListN = length(monthList);

% 昼夜标记.
daynight = {'Day', 'Night'};
daynight = daynight{flg1};

% 各统计变量名称.
staVarNameList = {'monthSampleRatioArray', 'fixedR2MonthArray', 'fixedRmseMonthArray'};
staVarNameListN = length(staVarNameList);

% 各统计图名称.
staNameList = {'Samples Ratio', 'R2', 'RMSE'};

% 直方图各统计变量的分级.
edgesList = {[0, 1, 3:3:15], 0:0.2:1, 0:6};

% 直方图各统计变量图表属性.
yLimList = {[0 40], [0 60], [0 60]};  % Y轴范围.
staTextY = [36, 56, 56];  % 精度评价指标Y轴位置.
monthTextY = [-3, -6, -5];  % 月份标注Y轴位置.
staNameX = [96, 84, 96];  % 统计指标名称X轴位置.

% 空间分布图各变量图表属性.
meanValueRange = {[0 15], [0 1], [0 5]};  % 均值取值范围.
stdValueRange = {[0 15], [0 0.2], [0 1]}; % 标准差取值范围.
meanValueTicks = {0:3:15, 0:0.2:1, 0:5};
stdValueTicks = {0:3:15, 0:0.05:0.2, 0:0.2:1};


%% 路径.
% 根目录.
rootPath = 'F:\AMSR_MODIS_Fusion';
dataPath = fullfile(rootPath, 'Data');
figPath = fullfile(rootPath, 'Figures');

% 模型数据路径.
regressionMatPath = fullfile(dataPath, 'Regression_Matlab');

% CCSEV数据路径.
ccsevMatPath = fullfile(dataPath, 'CCSEV_Matlab');

% 分区影像数据路径.
zonesPath = fullfile(dataPath, 'Zones', 'GeographicalZones_62_Merged.tif');
[zonesLayer, zonesRef] = readgeoraster(zonesPath);
zonesLayer = single(zonesLayer);

% 图表输出路径.
modelEstPath = fullfile(figPath, 'ModelEstimation');
if ~exist(modelEstPath, 'dir')
    mkdir(modelEstPath)
end

%% 模型评估.
for i = 8 : yearListN
    yearStr = num2str(yearList(i));
    yearStr = '2013';

    % 模型Mat数据路径.
    regressionPureMatName = sprintf('Regression_Pure_%s_%s.mat', yearStr, daynight);
    regressionPureMatPath = fullfile(regressionMatPath, regressionPureMatName);
    if ~exist(regressionPureMatPath, 'file')
        error('%s年%s的模型文件不存在, 请检查!', yearStr, daynight);
    end

    % 分精度评价指标制图.
    for j = 1 : staVarNameListN
        staVarName = staVarNameList{j};
        staName = staNameList{j};
        load(regressionPureMatPath, 'pureLcCodeList', staVarName);

        % 不统计积雪覆盖区(601, 602, ..., 608).
        pureLcCodeIndex = pureLcCodeList < 600;
        staVarArray = eval(staVarName);
        staVarArray = staVarArray(pureLcCodeIndex, :);

        % 直方图制图.
        timestamp = sprintf('%s %s', yearStr, daynight);
        modelHistFigPath = fullfile(modelEstPath, sprintf('Hist_%s_%s.png', staName, timestamp));
        if j == 1, staName = replace(staName, ' ', '\n'); end
        if ~exist(modelHistFigPath, 'file')
            f1 = figure; f1.Position = [10 100 2000 500];
            edges = edgesList{j}; edgesN = length(edges);
            xtickNum = 1 : edgesN-1; staTextX = 1; monthTextX = 3;
            [staMean, staStd] = deal(zeros(monthListN, 1));
            for k = 1 : monthListN
                staMonth = staVarArray(:, k);
                staMean(k) = mean(staMonth);
                staStd(k) = std(staMean);
                yConut = histcounts(staMonth, edges);
                b1 = bar(xtickNum, yConut); b1.BarWidth = 1;

                meanStr = sprintf('Mean: %.2f', staMean(k));
                stdStr = sprintf('STD: %.2f', staStd(k));
                text(staTextX, staTextY(j), meanStr, FontSize=11, FontWeight='bold');
                text(staTextX, staTextY(j)-2, stdStr, FontSize=11, FontWeight='bold');
                text(monthTextX, monthTextY(j), monthList{k}, FontSize=11, FontWeight='bold');

                xtickNum = xtickNum + (edgesN + 1);
                staTextX = staTextX + (edgesN + 1);
                monthTextX = monthTextX + (edgesN + 1);

                hold on
            end
            xlabelStr = repmat([string(edges), ''], 1, monthListN);
            ax = gca; ax.FontWeight = 'bold'; grid on;
            ax.YLim = yLimList{j}; ax.YLabel.String = sprintf('Zones Number（%s）', daynight);
            ax.XTick = 0.5:length(xlabelStr); ax.XTickLabel = xlabelStr; ax.XTickLabelRotation = 90;
            if j == 1, staName = replace(staName, ' ', '\n'); end
            text(staNameX(j), 0, sprintf(staName), FontWeight='bold');

            exportgraphics(f1, modelHistFigPath);
            close all
        end

        % 空间分布图.
        if j == 1, staName = replace(staName, '\n', ' '); end
        modelZonesMeanFigName = sprintf('Zones_%s_Mean_%s.png', staName, timestamp);
        modelZonesMeanFigPath = fullfile(modelEstPath, modelZonesMeanFigName);
        modelZonesStdFigName = sprintf('Zones_%s_Std_%s.png', staName, timestamp);
        modelZonesStdFigPath = fullfile(modelEstPath, modelZonesStdFigName);
        if ~exist(modelZonesMeanFigPath, 'file') || ~exist(modelZonesStdFigPath, 'file')
            ccsevMatFilePath = fullfile(ccsevMatPath, sprintf('CCSEV_%s', yearStr));
            load(ccsevMatFilePath, 'zonesLcLayer');
            pureLcCodeList = pureLcCodeList(pureLcCodeIndex);
            pureLcCodeListN = length(pureLcCodeList);

            staMeanZonesLayer = single(zonesLcLayer);
            staMeanZonesLayer(staMeanZonesLayer == 128) = nan;
            zonesLayer = mod(staMeanZonesLayer, 100);
            zonesLayer(staMeanZonesLayer <= 200) = staMeanZonesLayer(staMeanZonesLayer <= 200);
            staStdZonesLayer = staMeanZonesLayer;

            staVarMeanVector = mean(staVarArray, 2);
            staVarStdVector = std(staVarArray, 0, 2);
            for k = 1 : pureLcCodeListN
                staMeanZonesLayer(zonesLayer == pureLcCodeList(k)) = staVarMeanVector(k);
                staStdZonesLayer(zonesLayer == pureLcCodeList(k)) = staVarStdVector(k);
            end
%             geotiffwrite(modelZonesMeanFigPath, staMeanZonesLayer, zonesRef);
%             geotiffwrite(modelZonesStdFigPath, staStdZonesLayer, zonesRef);

            f2 = figure;
            imshow(staMeanZonesLayer, meanValueRange{j}, Colormap=turbo);
            colorbar(FontSize=12, FontWeight='bold', Ticks=meanValueTicks{j});
            titleStr = sprintf('Annual Mean of %s in %s', staNameList{j}, daynight);
            title(titleStr, FontSize=12, FontWeight='bold');
            exportgraphics(f2, modelZonesMeanFigPath);

            f3 = figure;
            imshow(staStdZonesLayer, stdValueRange{j}, Colormap=turbo);
            colorbar(FontSize=12, FontWeight='bold', Ticks=stdValueTicks{j});
            titleStr = sprintf('Annual STD of %s in %s', staNameList{j}, daynight);
            title(titleStr, FontSize=12, FontWeight='bold');
            exportgraphics(f3, modelZonesStdFigPath);

            close all
        end

    end

end

