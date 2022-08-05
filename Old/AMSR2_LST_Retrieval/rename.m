dataPath = 'F:\AMSR_MODIS_Fusion\Data\CCSEV_Matlab\';
dataList = dir([dataPath, 'CCSEV_*']);
dataList = {dataList.name}';
dataListN = length(dataList);

for i = 1 : dataListN
    dataName = dataList{i};
    yearStr = dataName(7:10);
    matPath = fullfile(dataPath, dataName);
%     dayNInYear = yeardays((str2double(yearStr)));
%     newDataName = replace(dataName, 'CCSEV', 'CCSEV2');
    if exist(matPath, "file")
        varList = who('-file', matPath);
        disp(yearStr)
        disp(varList)
        if ~ismember('dateList', varList)
            firstDate = datetime([yearStr, '0101'], 'InputFormat', 'yyyyMMdd');
            lastDate = datetime([yearStr, '1231'], 'InputFormat', 'yyyyMMdd');
            dateList = cellstr(datestr(firstDate: lastDate, 'yyyymmdd'));
            load(matPath, 'otherPctArray');
            if length(dateList) == size(otherPctArray, 3)
                save(matPath, 'dateList', '-append');
            else
                disp(yearStr);
            end
        end        
    end
end

