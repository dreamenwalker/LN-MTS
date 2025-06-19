import numpy as np
import cv2
from shapely.geometry import Polygon

def calculate_areas(img):

    # 二值化图像
    #_, thresh = cv2.threshold(img, 1, 255, cv2.THRESH_BINARY)
    thresh = img
    # 查找轮廓
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # 计算每个白色区域的面积
    areas = [cv2.contourArea(contour) for contour in contours]
    # 使用连通组件标记函数来标记不同的区域
    num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(img.astype(np.uint8), connectivity=8)
    pixel_counts = []
    for label in range(1, num_labels):  # 跳过背景标签 0
        area = stats[label, cv2.CC_STAT_AREA]
        pixel_counts.append(area)
    return pixel_counts

def calculate_centers(img, i):
    # 读取图像
    thresh = img

    # 二值化图像
    #_, thresh = cv2.threshold(img, 1, 255, cv2.THRESH_BINARY)

    # 查找轮廓
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    centers = []
    for contour in contours:
        # 计算轮廓的矩
        M = cv2.moments(contour)

        # 计算中心点
        if M["m00"] != 0:
            cx = int(M["m10"] / M["m00"])
            cy = int(M["m01"] / M["m00"])
            centers.append((cx, cy, i))

    return centers