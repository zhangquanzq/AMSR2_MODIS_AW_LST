# coding=utf-8
# 裁剪中国区域的MODIS积雪覆盖数据.
import arcpy
import os
from os import path
from glob import glob

# 预设参数.
cnExtent = '66 17 136 55'  # [xMin, yMin, xMax, yMax] in Lat and Lon.
yearList = range(2012, 2020)

# 路径.
rootPath = r'F:\AMSR_MODIS_Fusion'
dataPath = path.join(rootPath, 'Data')

modisSnowWorldHdfPath = path.join(dataPath, 'MYD10C1_1_World_HDF')
modisSnowCnTifPath = path.join(dataPath, 'MYD10C1_2_CN_TIF')
if not path.exists(modisSnowCnTifPath):
    os.mkdir(modisSnowCnTifPath)

# 按年份将HDF格式中的MODIS Snow数据保存为TIF格式, 并裁剪为中国范围.
for yearNum in yearList:
    # 创建存储输出数据的年度文件夹.
    modisSnowYearTifPath = path.join(modisSnowCnTifPath, 'MYD10C1_{0}XXX'.format(yearNum))
    if not path.exists(modisSnowYearTifPath):
        os.mkdir(modisSnowYearTifPath)

    # 处理年内所有日期的MODIS Snow数据.
    modisSnowYearHdfPath = path.join(modisSnowWorldHdfPath, 'MYD10C1_{0}XXX'.format(yearNum))
    modisSnowHdfPathList = glob(path.join(modisSnowYearHdfPath, 'MYD10C1*.hdf'))
    for modisSnowHdfPath in modisSnowHdfPathList:
        # 临时文件, 为从HDF中提取的全球范围TIF格式MODIS Snow文件.
        modisSnowTifTempPath = path.join(modisSnowCnTifPath, 'modisSnowTmp.tif')
        if arcpy.Exists(modisSnowTifTempPath):
            arcpy.Delete_management(modisSnowTifTempPath)

        # 从HDF文件中提取全球范围的TIF格式MODIS Snow数据, 并裁剪到中国范围, 删除全球范围文件.
        modisSnowTifName = path.basename(modisSnowHdfPath).replace('.hdf', '.tif')
        modisSnowTifPath = path.join(modisSnowYearTifPath, modisSnowTifName)
        if not arcpy.Exists(modisSnowTifPath):
            print(u'提取中国区MODIS积雪覆盖数据: {0}.'.format(modisSnowTifName))
            arcpy.ExtractSubDataset_management(modisSnowHdfPath, modisSnowTifTempPath, '0')
            arcpy.Clip_management(modisSnowTifTempPath, cnExtent, modisSnowTifPath)
            arcpy.Delete_management(modisSnowTifTempPath)
