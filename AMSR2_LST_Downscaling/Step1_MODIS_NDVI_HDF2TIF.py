# coding=utf-8
# 拼接研究区的MODIS地表覆盖数据.
import arcpy
import os
import numpy as np
from os import path
from glob import glob

# 预设参数.
modisType = 'MYD13A2'
region = 'CN'
yearList = range(2022, 2023)
layerNameList = ['NDVI', 'QA']
layerNumList = ['0', '2']

# ArcPy环境设置.
arcpy.env.overwriteOutput = True
arcpy.env.pyramid = 'NONE'

# 路径.
rootPath1 = r'K:\AMSR2_LST_Retrieval'
rootPath2 = r'K:\AMSR2_LST_Downscaling'
dataPath1 = path.join(rootPath1, 'Data')
dataPath2 = path.join(rootPath2, 'Data')
ndviTileRegionPath = path.join(dataPath2, '{0}_1_Tile{1}_HDF'.format(modisType, region))

ndviPrjRegionPath = path.join(dataPath2, '{0}_2_Prj{1}_TIF'.format(modisType, region))
if not path.exists(ndviPrjRegionPath):
    os.mkdir(ndviPrjRegionPath)

prjPath = path.join(dataPath2, 'MODIS_Sinusoidal.prj')
refRasPath = path.join(dataPath1, r'MYD11A1_2_PrjCN_TIF\MYD11A1_2012XXX_TIF',
                       'MYD11A1.A2012001.LST_Day.tif')

# 按年份提取HDF格式文件中的MODIS NDVI数据并保存为TIF格式, 然后拼接MODIS NDVI.
for yearNum in yearList:
    ndviYearTifPath = path.join(ndviPrjRegionPath,'{0}_{1}XXX_TIF'.format(modisType, yearNum))
    if not path.exists(ndviYearTifPath):
        os.mkdir(ndviYearTifPath)

    ndviHdfPathList = glob(path.join(ndviTileRegionPath, '{0}_{1}XXX_HDF'.
                                     format(modisType, yearNum), modisType + '*.hdf'))
    ndviDateList = [path.basename(ndviHdfPath)[9:16] for ndviHdfPath in ndviHdfPathList]
    ndviDateList = np.unique(ndviDateList)

    for ndviDate in ndviDateList:
        print(u'提取和镶嵌{0} {1}的{2}数据.'.format(ndviDate, region, modisType))
        # 获取当天的NDVI数据路径列表.
        ndviHdfPathList2 = [ndviHdfPath for ndviHdfPath in ndviHdfPathList if 'A' + ndviDate in
                            ndviHdfPath]
        # 创建临时文件夹.
        tempPath = path.join(rootPath2, 'temp_{0}_{1}'.format(modisType, ndviDate))
        if not path.exists(tempPath):
            os.mkdir(tempPath)

        # 提取和镶嵌各NDVI数据层.
        for i in range(len(layerNameList)):
            layer, layerNum = layerNameList[i], layerNumList[i]
            # 判断镶嵌后的MODIS NDVI文件是否存在.
            ndviPrjName = '{0}.A{1}.061_{2}.tif'.format(modisType, ndviDate, layer)
            ndviPrjPath = path.join(ndviYearTifPath, ndviPrjName)
            if arcpy.Exists(ndviPrjPath):
                continue

            # 从HDF文件中提取各NDVI数据层, 并保存为TIF格式.
            for ndviHdfPath in ndviHdfPathList2:
                ndviTifName = path.basename(ndviHdfPath).replace('.hdf', '_{0}.tif'.format(layer))
                ndviTifPath = path.join(tempPath, ndviTifName)
                if not arcpy.Exists(ndviTifPath):
                    arcpy.ExtractSubDataset_management(ndviHdfPath, ndviTifPath, layerNum)

            # 使用地理数据库的镶嵌数据集拼接MODIS LC数据, 并导出.
            gdbName = 'modisMosaic.gdb'
            mosaicName = layer + 'Mosaic'
            gdbPath = path.join(tempPath, gdbName)
            mosaicPath = path.join(gdbPath, mosaicName)
            if arcpy.Exists(gdbPath):
                arcpy.Delete_management(gdbPath)
            arcpy.CreateFileGDB_management(tempPath, gdbName)
            arcpy.CreateMosaicDataset_management(gdbPath, mosaicName, prjPath)
            arcpy.AddRastersToMosaicDataset_management(mosaicPath, 'Raster Dataset', tempPath,
                                                       filter='*{0}.tif'.format(layer))

            ndviTempPath = ndviPrjPath.replace('.tif', '2.tif')
            arcpy.env.snapRaster = refRasPath
            arcpy.ProjectRaster_management(mosaicPath, ndviTempPath, refRasPath, 'NEAREST', '0.01')
            arcpy.ClearEnvironment('snapRaster')
            arcpy.Clip_management(ndviTempPath, '#', ndviPrjPath, refRasPath)
            arcpy.Delete_management(ndviTempPath)

        # 删除临时文件夹.
        arcpy.Delete_management(tempPath)
