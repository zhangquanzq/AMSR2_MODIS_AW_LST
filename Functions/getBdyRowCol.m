function [startRow, endRow, startCol, endCol] = getBdyRowCol(clipRef, baseRef)
% 获取待裁剪影像中与参考影像边界对齐的起始和结束行列号. 输入参数为两个影像的空间参考对象.
clipRowN = clipRef.RasterSize(1);
clipColN = clipRef.RasterSize(2);
clipCellsizeX = clipRef.CellExtentInLatitude;
clipCellsizeY = clipRef.CellExtentInLatitude;
clipXMin = clipRef.LongitudeLimits(1);
clipXMax = clipRef.LongitudeLimits(2);
clipYMin = clipRef.LatitudeLimits(1);
clipYMax = clipRef.LatitudeLimits(2);

baseXMin = baseRef.LongitudeLimits(1);
baseXMax = baseRef.LongitudeLimits(2);
baseYMax = baseRef.LatitudeLimits(2);
baseYMin = baseRef.LatitudeLimits(1);

if clipYMax > baseYMax
    startRow = 1 + ceil((clipYMax - baseYMax) / clipCellsizeY);
else
    startRow = 1;
end

if clipYMin < baseYMin
    endRow = clipRowN - ceil((baseYMin - clipYMin) / clipCellsizeY);
else
    endRow = clipRowN;
end

if baseXMin > clipXMin
    startCol = 1 + ceil((baseXMin - clipXMin) / clipCellsizeX);
else
    startCol = 1;
end

if baseXMax < clipXMax
    endCol = clipColN - ceil((clipXMax - baseXMax) / clipCellsizeX);
else
    endCol = clipColN;
end

end