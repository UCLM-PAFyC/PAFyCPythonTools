import os
import shutil
from os.path import exists
import datetime
import glob
from math import floor, ceil, sqrt, isnan, modf, trunc
import csv
import re
import cv2 as cv
import math
from pathlib import Path

input_dxf_file = "D:/Aicedrone/20230125_Rail/qAicedrone/paper/topographic_surveying/taquimetrico sorolla_original.dxf"
output_dxf_file = "D:/Aicedrone/20230125_Rail/qAicedrone/paper/topographic_surveying/taquimetrico_sorolla_gcps.dxf"
input_gcps = "D:/Aicedrone/20230125_Rail/qAicedrone/paper/topographic_surveying/topographic_surveying_check_points_offsets.csv"
gcps_max_dis = 1.5
tile_size = 10.
tiles_round = [-1, 0, 1]
points_by_tilex_tiley = {}
with open(input_gcps) as file:
    cont = 0
    for line in file:
        cont = cont + 1
        if cont == 1:
            continue
        values = line.split(',')
        x = float(values[1])
        y = float(values[2])
        z = float(values[3])
        z_target = float(values[4])
        tile_x = math.floor(x / tile_size)
        tile_y = math.floor(y / tile_size)
        if not tile_x in points_by_tilex_tiley:
            points_by_tilex_tiley[tile_x] = {}
            points_by_tilex_tiley[tile_x][tile_y] = []
        elif not tile_y in points_by_tilex_tiley[tile_x]:
            points_by_tilex_tiley[tile_x][tile_y] = []
        point = [x, y, z, z_target]
        points_by_tilex_tiley[tile_x][tile_y].append(point)
input_dxf = open(input_dxf_file, 'r')
output_dxf = open(output_dxf_file, "w")
cont_line = 0
while True:
    line = input_dxf.readline().strip()
    cont_line = cont_line + 1
    if line.replace(' ', '') == 'EOF':
        output_dxf.write(line)
        output_dxf.write('\n')
        break
    output_dxf.write(line)
    output_dxf.write('\n')
    value = line.replace(' ','')
    if value == '10':
        line = input_dxf.readline().strip()
        cont_line = cont_line + 1
        if line.replace(' ', '') == 'EOF':
            output_dxf.write(line)
            output_dxf.write('\n')
            break
        str_lines_to_write = []
        str_lines_to_write.append(line)
        x = float(line.replace(' ',''))
        line = input_dxf.readline().strip()
        str_lines_to_write.append(line)
        cont_line = cont_line + 1
        if line.replace(' ', '') == 'EOF':
            for str_line in str_lines_to_write:
                output_dxf.write(str_line)
                output_dxf.write('\n')
            break
        value = line.replace(' ','')
        if value == '20':
            line = input_dxf.readline().strip()
            str_lines_to_write.append(line)
            cont_line = cont_line + 1
            if line.replace(' ', '') == 'EOF':
                for str_line in str_lines_to_write:
                    output_dxf.write(str_line)
                    output_dxf.write('\n')
                break
            y = float(line.replace(' ',''))
            line = input_dxf.readline().strip()
            str_lines_to_write.append(line)
            cont_line = cont_line + 1
            if line.replace(' ', '') == 'EOF':
                for str_line in str_lines_to_write:
                    output_dxf.write(str_line)
                    output_dxf.write('\n')
                break
            value = line.replace(' ','')
            if value == '30':
                line = input_dxf.readline().strip()
                str_lines_to_write.append(line)
                cont_line = cont_line + 1
                if line.replace(' ', '') == 'EOF':
                    for str_line in str_lines_to_write:
                        output_dxf.write(str_line)
                        output_dxf.write('\n')
                    break
                z = float(line.replace(' ',''))
                new_z = None
                new_x = None
                new_y = None
                tile_x = math.floor(x / tile_size)
                tile_y = math.floor(y / tile_size)
                tiles_x = []
                tiles_y = []
                for tile_dis_x in tiles_round:
                    tile_x_near = tile_x + round(tile_size * tile_dis_x)
                    for tile_dis_y in tiles_round:
                        tile_y_near = tile_y + round(tile_size * tile_dis_y)
                        tiles_x.append(tile_x_near)
                        tiles_y.append(tile_y_near)
                for it in range(len(tiles_x)):
                    tile_x_i = tiles_x[it]
                    tile_y_i = tiles_y[it]
                    if tile_x_i in points_by_tilex_tiley:
                        if tile_y_i in points_by_tilex_tiley[tile_x_i]:
                            points = points_by_tilex_tiley[tile_x_i][tile_y_i]
                            min_dis_3d = 10000000.
                            for i in range(len(points)):
                                diff_x = x - points[i][0]
                                diff_y = y - points[i][1]
                                diff_z = z - points[i][3]
                                dis_2d = math.sqrt(diff_x ** 2. + diff_y ** 2.)
                                dis_3d = math.sqrt(diff_x ** 2. + diff_y ** 2. + diff_z ** 2.)
                                if dis_3d < gcps_max_dis:
                                    if dis_3d < min_dis_3d:
                                        min_dis_3d = dis_3d
                                        new_z = points[i][3]
                                        new_x = points[i][0]
                                        new_y = points[i][1]
                if new_z:
                    str_new_x = f'{new_x:.3f}'
                    str_new_y = f'{new_y:.3f}'
                    str_new_z = f'{new_z:.3f}'
                    str_lines_to_write[0] = str_new_x
                    str_lines_to_write[2] = str_new_y
                    str_lines_to_write[4] = str_new_z
        for str_line in str_lines_to_write:
            output_dxf.write(str_line)
            output_dxf.write('\n')
    # if cont_line == 100000:
    #     break
input_dxf.close()
output_dxf.close()

