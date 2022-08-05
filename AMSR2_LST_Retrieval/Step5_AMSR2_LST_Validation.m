%% AMSR2地表温度的验证.

%% 功能标记和预设参数.
% 指定白天和晚上的标记. 1表示白天(升轨), 2表示晚上(降轨).
flg1 = 1;

% 数据年份列表(时间区间2012/07/02-2019/12/31).
yearList = 2012 : 2019;
yearListN = length(yearList);

% 昼夜标记.
daynight = {'Day', 'Night'};
daynight = daynight{flg1};

%% 路径.
% 根目录.
rootPath = 'F:\AMSR_MODIS_Fusion';
dataPath = fullfile(rootPath, 'Data');
figPath = fullfile(rootPath, 'Figures');

% 输入数据路径.
modisLstMaskMatPath = fullfile(dataPath, 'MYD11A1_3_MaskCn_Matlab');
amsr2LstMatPath = fullfile(dataPath, 'AMSR2_4_LST_Matlab');
regressionMatPath = fullfile(dataPath, 'Regression_Matlab'); 

% 输出的统计指标数据路径.
amsr2ModisScatterPath = fullfile(figPath, 'AMSR2_MODIS_LST_Scatter');
if ~exist(amsr2ModisScatterPath, 'dir')
    mkdir(amsr2ModisScatterPath)
end

%% 统计作图.
for i = 1 : yearListN
    yearStr = num2str(yearList(i));
    yearStr = '2013';

    amsr2ModisScatterYearPath = fullfile(amsr2ModisScatterPath, yearStr);
    if ~exist(amsr2ModisScatterYearPath, 'dir')
        mkdir(amsr2ModisScatterYearPath)
    end

    % 从Mat文件中读取Mask后的MODIS LST数据.
    modisLstMaskFileName = sprintf('MYD11A1_MaskCn_%s_%s.mat', yearStr, daynight);
    modisLstMaskFilePath = fullfile(modisLstMaskMatPath, modisLstMaskFileName);
    load(modisLstMaskFilePath, 'modisLstMaskYearArray');

    % 从Mat文件中读取反演的AMSR2 LST数据.
    amsr2LstCnYearMatName = sprintf('AMSR2_LstCn_%s_%s.mat', daynight, yearStr);
    amsr2LstCnYearMatPath = fullfile(amsr2LstMatPath, amsr2LstCnYearMatName);
    load(amsr2LstCnYearMatPath, 'amsr2LstCnYearArray', 'validDateList', 'amsr2Ref');

    % !!! 逐步回归得到的AMSR2 LST值与MODIS LST的散点图 !!!
    if strcmp(yearStr, '2013')
        regressionPureMatName = sprintf('Regression_Pure_%s_%s.mat', yearStr, daynight);
        regressionPureMatPath = fullfile(regressionMatPath, regressionPureMatName);
        load(regressionPureMatPath, 'amsr2LstMaskYearArray2');
        amsr2LstCnYearArray = amsr2LstMaskYearArray2;
    end

    fprintf('输出%s年AMSR2和MODIS地表温度的散点图.\n', yearStr);
    % 输出AMSR2和MODIS LST之间的年度散点图.
    timestamp = sprintf('%s %s', yearStr, daynight);
    amsr2ModisLstScatterName = sprintf('AMSR2_MODIS_LST_Scatter_%s.tif', timestamp);
    amsr2ModisLstScatterPath = fullfile(amsr2ModisScatterPath, amsr2ModisLstScatterName);
    if ~exist(amsr2ModisLstScatterPath, 'file')
        validIndexVector = find((amsr2LstCnYearArray ~= 0) & ~isnan(modisLstMaskYearArray));
        amsr2LstCnVector = amsr2LstCnYearArray(validIndexVector);
        modisLstMaskVector = modisLstMaskYearArray(validIndexVector);

        f = lstScatter(amsr2LstCnVector, modisLstMaskVector, timestamp);
        exportgraphics(f, amsr2ModisLstScatterPath);
        close all
    end

    % 输出AMSR2和MODIS LST之间的每日散点图.
    for j = 1 : length(validDateList)
        timestamp = sprintf('%s %s', validDateList{j}, daynight);
        amsr2ModisLstScatterName = sprintf('AMSR2_MODIS_LST_Scatter_%s.tif', timestamp);
        amsr2ModisLstScatterPath = fullfile(amsr2ModisScatterYearPath, amsr2ModisLstScatterName);
        if ~exist(amsr2ModisLstScatterPath, 'file')
            amsr2LstCnDailyLayer = amsr2LstCnYearArray(:, :, j);
            modisLstMaskDailyLayer = modisLstMaskYearArray(:, :, j);

            validIndexVector = find((amsr2LstCnDailyLayer ~= 0) & ~isnan(modisLstMaskDailyLayer));
            amsr2LstCnVector = amsr2LstCnDailyLayer(validIndexVector);
            modisLstMaskVector = modisLstMaskDailyLayer(validIndexVector);

            f = lstScatter(amsr2LstCnVector, modisLstMaskVector, timestamp);
            exportgraphics(f, amsr2ModisLstScatterPath);
            close all
        end
    end

end

%% 自定义函数.
% AMSR2和MODIS地表温度的散点图.
function f = lstScatter(amsr2LstVector, modisLstVector, timestamp)
lstBias = mean(amsr2LstVector - modisLstVector);
lstMAE = mean(abs(amsr2LstVector - modisLstVector));
lstR = corrcoef(amsr2LstVector, modisLstVector);
lstR2 = lstR(1, 2) .^ 2;
lstRMSE = sqrt(sum((amsr2LstVector - modisLstVector).^2) / length(amsr2LstVector));

f = figure; f.Visible = false;
plot(amsr2LstVector, modisLstVector, '.k', [220, 360], [220, 360], 'r');
xlabel('AMSR2 LST (K)'); ylabel('MODIS LST (K)');

txt0 = timestamp;
txt1 = ['N: ', num2str(sum(~isnan(modisLstVector)))];
txt2 = ['R^2: ', num2str(lstR2, '%.3f')];
txt3 = ['Bias: ', num2str(lstBias, '%.3f')];
txt4 = ['MAE: ', num2str(lstMAE, '%.3f')];
txt5 = ['RMSE: ', num2str(lstRMSE, '%.3f')];
text(0.7, 0.35, txt0, Units='normalized', FontSize=12, Color='r');
text(0.7, 0.29, txt1, Units='normalized', FontSize=12);
text(0.7, 0.23, txt2, Units='normalized', FontSize=12);
text(0.7, 0.17, txt3, Units='normalized', FontSize=12);
text(0.7, 0.11, txt4, Units='normalized', FontSize=12);
text(0.7, 0.05, txt5, Units='normalized', FontSize=12);
end
