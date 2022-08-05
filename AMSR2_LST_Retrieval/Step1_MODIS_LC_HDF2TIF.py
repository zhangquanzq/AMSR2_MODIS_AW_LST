# coding=utf-8
# 拼接研究区的MODIS地表覆盖数据.
import arcpy
import os
from os import path
from glob import glob

# 预设参数.
modisType = 'MCD12Q1'
region = 'CN'
yearList = range(2012, 2020)

# ArcPy环境设置.
arcpy.env.overwriteOutput = True
arcpy.env.pyramid = 'NONE'

# 路径.
rootPath = r'F:\AMSR_MODIS_Fusion'
dataPath = path.join(rootPath, 'Data')
modisLcTileRegionPath = path.join(dataPath, '{0}_1_Tile{1}_HDF'.format(modisType, region))

modisLcMosaicRegionPath = path.join(dataPath, '{0}_2_Mosaic{1}_TIF'.format(modisType, region))
if not path.exists(modisLcMosaicRegionPath):
    os.mkdir(modisLcMosaicRegionPath)

prjPath = path.join(dataPath, 'MODIS_Sinusoidal.prj')

# 按年份提取HDF格式文件中的MODIS LC数据并保存为TIF格式, 然后拼接MODIS LC.
for yearNum in yearList:
    # 判断镶嵌后的MODIS LC数据是否存在.
    modisLcMosaicName = '{0}.A{1}001.006.tif'.format(modisType, yearNum)
    modisLcMosaicPath = path.join(modisLcMosaicRegionPath, modisLcMosaicName)
    if arcpy.Exists(modisLcMosaicPath):
        continue

    # 创建临时文件夹.
    tempPath = path.join(rootPath, 'temp_{0}_{1}'.format(modisType, yearNum))
    if not path.exists(tempPath):
        os.mkdir(tempPath)

    # 从HDF文件中提取LC数据, 并保存为TIF格式.
    print(u'提取{0}年{1}的{2}数据.'.format(yearNum, region, modisType))
    modisLcYearFolderName = '{0}_{1}XXX_HDF'.format(modisType, yearNum)
    modisLcYearFolderPath = path.join(modisLcTileRegionPath, modisLcYearFolderName)
    modisLcHdfPathList = glob(path.join(modisLcYearFolderPath, '*.hdf'))
    for modisLcHdfPath in modisLcHdfPathList:
        modisLcTifPath = path.join(tempPath, path.basename(modisLcHdfPath).replace('hdf', 'tif'))
        if not arcpy.Exists(modisLcTifPath):
            arcpy.ExtractSubDataset_management(modisLcHdfPath, modisLcTifPath, '0')

    # 使用地理数据库的镶嵌数据集拼接MODIS LC数据, 并导出.
    print(u'镶嵌{0}年{1}的{2}数据.'.format(yearNum, region, modisType))
    gdbName = 'modisMosaic.gdb'
    mosaicName = 'modisMosaic'
    gdbPath = path.join(tempPath, gdbName)
    mosaicPath = path.join(gdbPath, mosaicName)
    if arcpy.Exists(gdbPath):
        arcpy.Delete_management(gdbPath)
    arcpy.CreateFileGDB_management(tempPath, gdbName)
    arcpy.CreateMosaicDataset_management(gdbPath, mosaicName, prjPath)
    arcpy.AddRastersToMosaicDataset_management(mosaicPath, 'Raster Dataset', tempPath,
                                               filter='*.tif')
    arcpy.CopyRaster_management(mosaicPath, modisLcMosaicPath)

    # 删除临时文件夹.
    arcpy.Delete_management(tempPath)

# 投影, 重分类, 重采样.
arcpy.env.workspace = modisLcMosaicRegionPath
modisLcMosaicList = arcpy.ListRasters('MCD12Q1*006.tif')
for modisLcMosaic in modisLcMosaicList:
    print(u'投影，重分类，重采样 {0}'.format(modisLcMosaic))
    snapDir = path.join(rootPath, r'AMSR2_2_CN_TIF\L3.TB6GHz_10\2012\07')
    snapRasterPath = path.join(snapDir, 'GW1AM2_20120700_01M_EQMA_L3SGT06HA2220220_BtH.tif')
    arcpy.env.snapRaster = snapRasterPath
    modisLcGcs = modisLcMosaic.replace('.tif', '_gcs.tif')
    modisLcReclassify = modisLcGcs.replace('.tif', '_reclsfy.tif')
    modisLcResmpl = modisLcReclassify.replace('.tif', '_resmpl.tif')
    if not arcpy.Exists(modisLcGcs):
        arcpy.ProjectRaster_management(modisLcMosaic, modisLcGcs, arcpy.SpatialReference(4326),
                                       'NEAREST', '0.005')
    if not arcpy.Exists(modisLcReclassify):
        lcRemap = arcpy.sa.RemapValue([[1, 2], [2, 2], [3, 2], [4, 2], [5, 2], [6, 2], [7, 2],
                                       [8, 2], [9, 2], [10, 2], [11, 3], [12, 2], [13, 5], [14, 2],
                                       [15, 4], [16, 1], [17, 3]])
        arcpy.sa.Reclassify(modisLcGcs, 'Value', lcRemap, 'NODATA').save(modisLcReclassify)
    if not arcpy.Exists(modisLcResmpl):
        arcpy.env.extent = arcpy.Raster(snapRasterPath).extent
        arcpy.Resample_management(modisLcReclassify, modisLcResmpl, '0.01', 'NEAREST')
        arcpy.ClearEnvironment('extent')
