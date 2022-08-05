# coding=utf-8
# 将AMSR2 H5格式的各通道亮温数据输出为TIF格式.
# 此Python脚本与Matlab脚本'Step1_AMSR2_BT_H5toTIF_Mat.m'的功能重复, 但运行速率低, 建议仅在没有Matlab时使用.
import arcpy
import numpy as np
import os
from os import path
from glob import glob

# 预设参数.
# AMSR2数据的极化.
polarize = ['H', 'V']

# 输出AMSR2 TIF文件的坐标系统.
amsr2Ref = arcpy.SpatialReference(4326)  # WGS 1984

# ArcPy环境设置.
arcpy.env.overwriteOutput = True
arcpy.env.pyramid = 'NONE'

# 路径.
rootPath = r'E:\AMSR_MODIS_Fusion\Data'
amsr2H5Path = path.join(rootPath, 'AMSR2_1_H5')

amsr2WorldTifPath = path.join(rootPath, 'AMSR2_2_World_TIF')
if not path.exists(amsr2WorldTifPath):
    os.mkdir(amsr2WorldTifPath)

# 将AMSR2 H5格式转为TIF格式.
for amsr2BandFolderPath in glob(path.join(amsr2H5Path, 'L3*')):
    amsr2TifBandFolderPath = path.join(amsr2WorldTifPath, path.basename(amsr2BandFolderPath))
    if not path.exists(amsr2TifBandFolderPath):
        os.mkdir(amsr2TifBandFolderPath)

    # 读取AMSR2 H5年份文件夹列表.
    for amsr2YearFolderPath in glob(path.join(amsr2BandFolderPath, '*')):
        amsr2TifYearFolderPath = path.join(amsr2TifBandFolderPath, path.basename(amsr2YearFolderPath))
        if not path.exists(amsr2TifYearFolderPath):
            os.mkdir(amsr2TifYearFolderPath)

        # 读取AMSR2 H5月份文件夹列表.
        for amsr2MonthFolderPath in glob(path.join(amsr2YearFolderPath, '*')):
            amsr2TifMonthFolderPath = path.join(amsr2TifYearFolderPath, path.basename(amsr2MonthFolderPath))
            if not path.exists(amsr2TifMonthFolderPath):
                os.mkdir(amsr2TifMonthFolderPath)

            # 读取AMSR2 H5每日文件列表.
            for amsr2H5DailyPath in glob(path.join(amsr2MonthFolderPath, 'GW1AM2*01D*.h5')):
                amsr2H5DailyName = path.basename(amsr2H5DailyPath)

                # 将AMSR2 H5的H, V极化数据转为TIF格式.
                for i in range(len(polarize)):
                    amsr2BtTifDailyName = amsr2H5DailyName.replace('.h5', '_Bt{0}.tif'.format(polarize[i]))
                    amsr2BtTifDailyPath = path.join(amsr2TifMonthFolderPath, amsr2BtTifDailyName)
                    if not path.exists(amsr2BtTifDailyPath):
                        try:
                            arcpy.ExtractSubDataset_management(amsr2H5DailyPath, amsr2BtTifDailyPath, str(i))
                            amsr2BtArray = arcpy.RasterToNumPyArray(amsr2BtTifDailyPath)
                            amsr2BtArray = np.hsplit(amsr2BtArray, 2)
                            amsr2BtArray = np.hstack((amsr2BtArray[1], amsr2BtArray[0]))
                            amsr2BtHRaster = arcpy.NumPyArrayToRaster(amsr2BtArray, arcpy.Point(-180, -90), 0.1, 0.1)
                            amsr2BtHRaster.save(amsr2BtTifDailyPath)
                            arcpy.DefineProjection_management(amsr2BtTifDailyPath, amsr2Ref)
                            print(amsr2BtTifDailyName)
                        except:
                            print('有问题数据: {0}'.format(amsr2BtTifDailyName))
