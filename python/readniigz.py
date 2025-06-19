# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 22:38:43 2023

@author: dreamen-i9
"""
import os
import SimpleITK
import numpy as np
from skimage.transform import resize
import pandas as pd
import cv2
import math
from calculate import calculate_areas, calculate_centers

# dcm_path = 
# dcm = SimpleITK.GetArrayFromImage(SimpleITK.ReadImage(dcm_path))[ind, :, :]

pathdir = r"ROI"
wrongpatient = []
dic = dict()
result = dict()
result2 = []
for patientdirname in os.listdir(pathdir):
    filepath = os.path.join(pathdir, patientdirname)
    for file in os.listdir(filepath):
        if "ROI" in file:
            ROIpath = os.path.join(filepath,file)
            ROI = SimpleITK.GetArrayFromImage(SimpleITK.ReadImage(ROIpath))
            #print(ROI.shape)
        elif "T2" in file:
            T2path = os.path.join(filepath,file)
            T2dcm = SimpleITK.GetArrayFromImage(SimpleITK.ReadImage(T2path))
            #print(T2dcm.shape)
        else:
            wrongpatient.append(file)

    '''
    淋巴结处理
    '''
    areaList = [] #保存各层面积
    centerList = [] #保存集合中心
    LINBA = ROI.astype(np.uint8)
    for i in range(ROI.shape[0]):
        # 首先进行图片二值化
        for j in range(512):
            for k in range(512):
                if(ROI[i,j,k]==2): LINBA[i,j,k]=1
                else: LINBA[i,j,k]=0
        areaList += calculate_areas(LINBA[i,:,:])
        centerList += calculate_centers(LINBA[i,:,:], i)

    average_area = sum(areaList)/len(areaList) #平均面积
    total_area = sum(areaList) #总面积
    max_area = max(areaList)
    min_area = min(areaList)
    dic[patientdirname] = [max_area, min_area, average_area, total_area, centerList]

    areaList2 = []  # 保存各层面积
    centerList2 = []  # 保存集合中心
    CANCER = ROI.astype(np.uint8)
    for i in range(CANCER.shape[0]):
        # 首先进行图片二值化
        for j in range(512):
            for k in range(512):
                if(ROI[i,j,k]==1): LINBA[i,j,k]=1
                else: LINBA[i,j,k]=0
        areaList2 += calculate_areas(CANCER[i,:,:])
        centerList2 += calculate_centers(CANCER[i,:,:], i)

    index = 0
    maxA = 0
    f = 0
    for area in areaList2:
        if(area>maxA):
            maxA = area
            index = f
        f += 1

    position = centerList2[index]

    distanceList = []
    for center in centerList:
        distance = math.sqrt(
            math.sqrt((center[0]-position[0])**2+(center[1]-position[1])**2)*math.fabs(center[2]-position[2])
        )
        distanceList.append(distance)

    max_d = max(distanceList)
    mean_d = sum(distanceList)/len(distanceList)
    min_d = min(distanceList)
    sum_d = sum(distanceList)

    result2.append([patientdirname, max_area, min_area, average_area, total_area, max_d, min_d, mean_d, sum_d])



'''
    T2dcm = T2dcm.astype(np.uint8)
    for i in range(T2dcm.shape[0]):
        # 首先进行图片二值化
        for j in range(512):
            for k in range(512):
                if (T2dcm[i, j, k] == 1):
                    T2dcm[i, j, k] = 1
                else:
                    T2dcm[i, j, k] = 0
        ret, labels, stats, centroids = cv2.connectedComponentsWithStats(T2dcm[i, :, :], connectivity=8)
        # 求最大区域
        maxT2 = 0
        max_index = 0
        f = 0
        for stat in stats:
            if(stat[4]>maxT2):
                maxT2 = stat[4]
                max_index = f
            f += 1
        T2_position = [centroids[max_index, 0], centroids[max_index, 1]]

        ROI_positions = dic[patientdirname][4]

        def distance(a_list, b_list):
            return math.sqrt((a_list[0]-b_list[0])**2+(a_list[1]-b_list[1])**2)

        distance_list = []
        for ROI_position in ROI_positions:
            distance_list.append(distance(T2_position, ROI_position))

    sum_distance = sum(distance_list)
    max_distance = max(distance_list)
    min_distance = min(distance_list)
    mean_distance = np.mean(distance_list)

    result[patientdirname] = [max_area, average_area, total_area, min_area,
                              max_distance, mean_distance, sum_distance, min_distance]
    result2.append([patientdirname, dic[patientdirname][0], dic[patientdirname][1], dic[patientdirname][2], dic[patientdirname][3],
                    max_distance, mean_distance, sum_distance, min_distance])
'''
#%%
import csv
column=['patientname', 'Lmax', 'Lmin', 'Lmean', 'Lsum', 'Dmax', 'Dmin', 'Dmean', 'Dsum']
res = pd.DataFrame(columns=column, data=result2)
res.to_csv('test.csv')



#print(dic)



        
                                        
            
    




'''
seg_path = r'F:\3-guanxu-colorectal\data\MRI\MR1712280010.nii.gz'

seg = SimpleITK.GetArrayFromImage(SimpleITK.ReadImage(seg_path))


AA = "train80_snr5"

if 'train' and '80' and 'snr5' in AA:
    print("这三个字符串同时存在于AA中")
else:
    print("这三个字符串不同时存在于AA中")
#%%
 # 导入所需库
import numpy as np
import cv2
import pandas as pd

 # 读取CT图像并转换为数组
ct_image = cv2.imread('ct_image.png', 0)

 # 定义肿瘤和淋巴结的像素值
tumor_value = 1
lymph_node_value = 2

 # 计算肿瘤和淋巴结的勾画区域

 # 计算平均距离和距离的和

 # 计算淋巴结勾画区域的面积和平均面积

 # 将特征保存到CSV文件中
data = {'Feature': ['Average Distance', 'Sum of Distances', 'Lymph Node Area', 'Average Area'],
         'Value': [average_distance, sum_of_distances, lymph_node_area, average_area]}
df = pd.DataFrame(data)
df.to_csv('image_features.csv', index=False)


 # 导入所需库
import numpy as np
import pandas as pd
from PIL import Image

 # 读取CT图像并转换为数组
ct_image = np.array(Image.open('ct_image.png'))

 # 定义肿瘤和淋巴结的像素值
tumor_value = 1
lymph_node_value = 2

 # 计算肿瘤和淋巴结的勾画区域
tumor_area = np.sum(ct_image == tumor_value)
lymph_node_area = np.sum(ct_image == lymph_node_value)

 # 计算不同淋巴结和肿瘤平均距离和距离之和的代码是什么？
 # 这部分需要使用图像处理算法来计算，具体实现需要根据您选择的方法来编写。

 # 计算淋巴结勾画区域的总面积和平均面积的代码是什么
 # 同样，计算面积需要根据具体的方法来实现，可以使用图像处理库中的函数来求解。



 # 导入所需库
import numpy as np
import pandas as pd
from PIL import Image
import cv2

 # 读取CT图像并转换为数组
ct_image = np.array(Image.open('ct_image.png'))

 


import os
import pydicom

 # 文件夹路径
folder_path = 'your_folder_path'

 # 遍历文件夹中的每个患者文件夹
for patient_folder in os.listdir(folder_path):
    patient_folder_path = os.path.join(folder_path, patient_folder)
    for file in os.listdir(patient_folder_path):
        if file.endswith('.dcm'):
            file_path = os.path.join(patient_folder_path, file)
            ds = pydicom.dcmread(file_path)
            # 对每个DICOM文件进行处理
             
 # 导入所需库
import numpy as np
import pandas as pd
from PIL import Image
import cv2

 # 读取CT图像并转换为数组
ct_image = np.array(Image.open('ct_image.png'))


 # 计算不同淋巴结和肿瘤平均距离和距离之和的代码需要根据具体的图像处理算法来编写。

 # 计算淋巴结勾画区域的总面积和平均面积的代码需要根据具体的方法来实现，可以使用图像处理库中的函数来求解。




             
import numpy as np
from scipy.spatial.distance import cdist

 # 找到像素值为1和2的位置
tumor_indices = np.argwhere(ct_image == 1)
lymph_indices = np.argwhere(ct_image == 2)

 # 计算平均距离和距离之和
average_distance_tumor = np.mean(cdist(tumor_indices, tumor_indices))
sum_of_distances_tumor = np.sum(cdist(tumor_indices, tumor_indices))

average_distance_lymph = np.mean(cdist(lymph_indices, lymph_indices))
sum_of_distances_lymph = np.sum(cdist(lymph_indices, lymph_indices))








#example
import numpy as np
from scipy import ndimage

# 创建一个示例的三维矩阵（假设为matrix）
matrix = np.array([[[1, 2, 2], [1, 1, 2]], [[2, 2, 3], [3, 3, 1]], [[1, 1, 1], [2, 2, 2]]])

# 找到像素值为2的区域位置
mask = (matrix == 2)
labeled, num_features = ndimage.label(mask)

# 计算不同像素值为1的区域的面积
areas = []
for i in range(1, num_features + 1):
    area = np.sum(labeled == i)
    areas.append(area)

# 计算所有面积之和
total_area = sum(areas)

print("不同像素值为1的区域的面积：", areas)
print("所有面积之和：", total_area)
'''





