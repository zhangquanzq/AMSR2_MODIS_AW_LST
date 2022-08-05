function [startBlockBdy, blockSize] = getStartBlockRowCol(ras1Ref, ras2Ref)
%GETSTARTBLOCKROWCOL 获取影像左上角第一个滑动窗口边界的行列号.
% ras1Ref: 进行窗口计算的影像数据的空间参考.
% ras2Ref: 确定滑动窗口大小的影像数据的空间参考. 需要保证ras1的空间分辨率是ras2的整数倍.
% startBlockBdy: 左上角第一个滑动窗口在ras1上的位置, 定义为: [topRow, bottomRow, leftCol, rightCol].
% blockSize: 滑动窗口的行列数. 定义为: [blockRowN, blockColN]

ras1CellsizeX = ras1Ref.CellExtentInLongitude;
ras1CellsizeY = ras1Ref.CellExtentInLatitude;
ras1XMin = ras1Ref.LongitudeLimits(1);
ras1XMax = ras1Ref.LongitudeLimits(2);
ras1Min = ras1Ref.LatitudeLimits(1);
ras1Max = ras1Ref.LatitudeLimits(2);
ras1PixelLeftSideXVector = ras1XMin : ras1CellsizeX : ras1XMax - ras1CellsizeX;
ras1PixelTopSideYVector = ras1Max : - ras1CellsizeY : ras1Min + ras1CellsizeY;

ras2RowN = ras2Ref.RasterSize(1);
ras2ColN = ras2Ref.RasterSize(2);
ras2CellsizeX = ras2Ref.CellExtentInLongitude;
ras2CellsizeY = ras2Ref.CellExtentInLatitude;
ras2XMin = ras2Ref.LongitudeLimits(1);
ras2XMax = ras2Ref.LongitudeLimits(2);
ras2YMin = ras2Ref.LatitudeLimits(1);
ras2YMax = ras2Ref.LatitudeLimits(2);
ras2PixelLeftSideXVector = ras2XMin : ras2CellsizeX : ras2XMax - ras2CellsizeX;
ras2PixelTopSideYVector = ras2YMax : - ras2CellsizeY : ras2YMin + ras2CellsizeY;

blockColN = round(ras2CellsizeX / ras1CellsizeX);
blockRowN = round(ras2CellsizeY / ras1CellsizeY);

for i = 1 : ras2ColN
    leftSideXDistance = abs(ras1PixelLeftSideXVector - ras2PixelLeftSideXVector(i));
    minLeftSideDistance = min(leftSideXDistance);
    if minLeftSideDistance <= ras1CellsizeX
        startLeftCol = find(leftSideXDistance == minLeftSideDistance);
        break
    end
end
startRightCol = startLeftCol + blockColN - 1;

for i = 1 : ras2RowN
    topSideYDistance = abs(ras1PixelTopSideYVector - ras2PixelTopSideYVector(i));
    minTopSideYDistance = min(topSideYDistance);
    if minTopSideYDistance <= ras1CellsizeY
        startTopRow = find(topSideYDistance == minTopSideYDistance);
        break
    end
end
startBottomRow = startTopRow + blockRowN - 1;

startBlockBdy = [startTopRow, startBottomRow, startLeftCol, startRightCol];
blockSize = [blockRowN, blockColN];
end
