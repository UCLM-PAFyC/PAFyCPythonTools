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
import random

random.seed()
input_dxf_file = "D:/Aicedrone/20230125_Rail/qAicedrone/paper/topographic_surveying/taquimetrico_sorolla_gcps.dxf"
output_dxf_file = "D:/Aicedrone/20230125_Rail/qAicedrone/paper/topographic_surveying/taquimetrico_sorolla.dxf"
input_rails = "D:/Aicedrone/20230125_Rail/qAicedrone/paper/topographic_surveying/rails/rails_vertices_offsets.csv"
rails_max_dis_2d = 0.09
tile_size = 10.
tiles_round = [-1, 0, 1]
dis_3d_tolerance = 0.02
z_rand_offset = 0.04
z_rand = 0.015
points_by_tilex_tiley = {}
with open(input_rails) as file:
    cont = 0
    for line in file:
        cont = cont + 1
        if cont == 1:
            continue
        values = line.split(',')
        dis = float(values[0])
        if dis > rails_max_dis_2d:
            continue
        z_target = float(values[1])
        x = float(values[2])
        y = float(values[3])
        z = float(values[4])
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
        output_dxf.write(line)
        output_dxf.write('\n')
        x = float(line.replace(' ',''))
        line = input_dxf.readline().strip()
        cont_line = cont_line + 1
        if line.replace(' ', '') == 'EOF':
            output_dxf.write(line)
            output_dxf.write('\n')
            break
        output_dxf.write(line)
        output_dxf.write('\n')
        value = line.replace(' ','')
        if value == '20':
            line = input_dxf.readline().strip()
            cont_line = cont_line + 1
            if line.replace(' ', '') == 'EOF':
                output_dxf.write(line)
                output_dxf.write('\n')
                break
            output_dxf.write(line)
            output_dxf.write('\n')
            y = float(line.replace(' ',''))
            line = input_dxf.readline().strip()
            cont_line = cont_line + 1
            if line.replace(' ', '') == 'EOF':
                output_dxf.write(line)
                output_dxf.write('\n')
                break
            output_dxf.write(line)
            output_dxf.write('\n')
            value = line.replace(' ','')
            if value == '30':
                line = input_dxf.readline().strip()
                cont_line = cont_line + 1
                if line.replace(' ', '') == 'EOF':
                    output_dxf.write(line)
                    output_dxf.write('\n')
                    break
                z = float(line.replace(' ',''))
                new_z = None
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
                            min_dis_2d = 10000000.
                            for i in range(len(points)):
                                diff_x = x - points[i][0]
                                diff_y = y - points[i][1]
                                dis_2d = math.sqrt(diff_x ** 2. + diff_y ** 2.)
                                if dis_2d < rails_max_dis_2d:
                                    if dis_2d < min_dis_2d:
                                        min_dis_2d = dis_2d
                                        new_z = points[i][3]
                                        chang_offset = z_rand_offset
                                        change_rand = (-1)**round(random.random())*(random.random() * z_rand)
                                        change = chang_offset + change_rand
                                        new_z = new_z + change
                if not new_z:
                    output_dxf.write(line)
                    output_dxf.write('\n')
                else:
                    str_new_z = f'{new_z:.3f}'
                    output_dxf.write(str_new_z)
                    output_dxf.write('\n')

    # if cont_line == 100000:
    #     break
input_dxf.close()
output_dxf.close()

