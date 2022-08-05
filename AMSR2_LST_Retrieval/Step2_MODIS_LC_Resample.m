%% 将MODIS LC升尺度到AMSR2 BT数据的0.1度分辨率.

%% 路径.
% 数据根目录.
rootPath = 'E:\AMSR_MODIS_Fusion';
dataPath = fullfile(rootPath, 'Data');
funcPath = fullfile(rootPath, 'Code\Functions');
addpath(funcPath);

% MODIS LC路径.
modisLcPath = fullfile(dataPath, 'MCD12Q1_2_CN_TIF');

% 中国范围栅格路径.
extentCnPath = fullfile(dataPath, 'Zones', 'ExtentCN_0d1.tif');

%% 读取数据参考与属性.
[extentCnLayer, extentCnRef] = readgeoraster(extentCnPath);
extentCnRowN = extentCnRef.RasterSize(1); extentCnColN = extentCnRef.RasterSize(2);

%% 升尺度MODIS LC.
modisLcGcsList = dir(fullfile(modisLcPath, 'MCD12Q1*gcs.tif'));
modisLcGcsList = {modisLcGcsList.name}';
modisLcGcsListN = length(modisLcGcsList);
for i = 1 : modisLcGcsListN
    modisLcName = modisLcGcsList{i};
    modisLcUpscaledName = replace(modisLcName, 'gcs.tif', 'gcs_upscaled.tif');
    modisLcUpscaledPath = fullfile(modisLcPath, modisLcUpscaledName); 
    if exist(modisLcUpscaledPath, 'file')
        continue
    end

    fprintf('升尺度MODIS LC: %s\n', modisLcName);
    modisLcFilePath = fullfile(modisLcPath, modisLcName);
    [modisLcLayer, modisLcRef] = readgeoraster(modisLcFilePath);
    lcRowN = modisLcRef.RasterSize(1); lcColN = modisLcRef.RasterSize(2);

    [lcBlockBdy, lcBlockSize] = getStartBlockRowCol(modisLcRef, extentCnRef);
    upscaledModisLcLayer = zeros(extentCnRowN, extentCnColN, 'uint8');
    for ii = 1 : extentCnRowN
        for jj = 1 : extentCnColN
            % 跳过ExtentCn图层的非中国区像元.
            if extentCnLayer(ii, jj) ~= 1
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
            lcBlock = modisLcLayer(lcBlockTopRow:lcBlockBottomRow, lcBlockLeftCol:lcBlockRightCol);
            upscaledModisLcLayer(ii, jj) = mode(lcBlock, 'all');
        end
    end
    geotiffwrite(modisLcUpscaledPath, upscaledModisLcLayer, extentCnRef, ...
        TiffTags=struct('Compression','LZW'));
end
