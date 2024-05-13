# authors:
# David Hernandez Lopez, david.hernandez@uclm.es
# Miguel Angel Moreno Hidalgo, miguelangel.moreno@uclm.es

import optparse
import numpy
from osgeo import gdal, osr, ogr
import os
import json
from urllib.parse import unquote
import shutil
from os.path import exists
import datetime
import glob
from math import floor, ceil, sqrt, isnan, modf, trunc
import csv
import re
import cv2 as cv


def copy_shapefile(input_shp, output_shp):
    str_error = ''
    input_base_name = os.path.splitext(os.path.basename(input_shp))[0]
    input_base_path = os.path.dirname(input_shp)
    output_base_path = os.path.dirname(output_shp)
    output_base_name = os.path.splitext(os.path.basename(output_shp))[0]
    for file in os.listdir(input_base_path):
        file_base_name = os.path.splitext(os.path.basename(file))[0]
        if file_base_name == input_base_name:
            file_extension = os.path.splitext(os.path.basename(file))[1]
            output_file = output_base_path + "/" + output_base_name  + file_extension
            output_file = os.path.normcase(output_file)
            input_file = input_base_path + "/" + file
            input_file = os.path.normcase(input_file)
            try:
                shutil.copyfile(input_file, output_file)
            except EnvironmentError as e:
                str_error = "Unable to copy file. %s" % e
                return str_error
    return str_error

class OptionParser(optparse.OptionParser):
    def check_required(self, opt):
        option = self.get_option(opt)
        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            self.error("%s option not supplied" % option)


def julian_date(day, month, year):
    if month <= 2:  # january & february
        year = year - 1.0
        month = month + 12.0
    jd = floor(365.25 * (year + 4716.0)) + floor(30.6001 * (month + 1.0)) + 2.0
    jd = jd - floor(year / 100.0) + floor(floor(year / 100.0) / 4.0)
    jd = jd + day - 1524.5
    # jd = jd + day - 1524.5 + (utc_time)/24.
    mjd = jd - 2400000.5
    return jd, mjd


def julian_date_to_date(jd):
    jd = jd + 0.5
    F, I = modf(jd)
    I = int(I)
    A = trunc((I - 1867216.25)/36524.25)
    if I > 2299160:
        B = I + 1 + A - trunc(A / 4.)
    else:
        B = I
    C = B + 1524
    D = trunc((C - 122.1) / 365.25)
    E = trunc(365.25 * D)
    G = trunc((C - E) / 30.6001)
    day = C - E + F - trunc(30.6001 * G)
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
    return year, month, day


def is_number(n):
    is_number = True
    try:
        num = float(n)
        # check for "nan" floats
        is_number = num == num  # or use `math.isnan(num)`
    except ValueError:
        is_number = False
    return is_number


def sortFunction(e):
    return e['value']

def process_gcc_or_vol(input_shp,
                       crop_minimum_value,
                       kmeans_clusters,
                       percentile_minimum_threshold,
                       gcc,
                       vol,
                       input_field_name,
                       date_format):
    str_error = None
    if kmeans_clusters < 0 and percentile_minimum_threshold < 0:
        str_error = "Function process"
        str_error += "\nkmeans_clusters or percentile_minimum_threshold must be greather than 0"
        return str_error
    elif kmeans_clusters > 0 and percentile_minimum_threshold > 0:
        str_error = "Function process"
        str_error += "\nkmeans_clusters or percentile_minimum_threshold must be greather than 0"
        return str_error
    if not exists(input_shp):
        str_error = "Function process"
        str_error += "\nNot exists file: {}".format(input_shp)
        return str_error
    driver = ogr.GetDriverByName('ESRI Shapefile')
    in_vec_ds = None
    try:
        in_vec_ds = driver.Open(input_shp, 1)  # 0 means read-only. 1 means writeable.
    except ValueError:
        str_error = "Function process"
        str_error += "\nError opening dataset file:\n{}".format(input_shp)
        return str_error
    in_layer = in_vec_ds.GetLayer()
    in_crs = in_layer.GetSpatialRef()
    in_crs_wkt = in_crs.ExportToWkt()
    in_geometry_type = in_layer.GetGeomType()
    if in_geometry_type != ogr.wkbPolygon \
            and in_geometry_type != ogr.wkbMultiPolygon \
            and in_geometry_type != ogr.wkbPolygonM and in_geometry_type != ogr.wkbPolygonZM:
        str_error = "Function process"
        str_error += "\nNot Polygon geometry type in file:\n{}".format(input_shp)
        return str_error
    in_layer_definition = in_layer.GetLayerDefn()
    input_field_id_index = in_layer_definition.GetFieldIndex(input_field_name)
    if input_field_id_index == -1:
        str_error = "Function process"
        str_error += "\nNot exists field: {} in file:\n{}".format(input_field_name, input_shp)
        return str_error
    output_field_name = None
    if gcc:
        output_field_name = '_dg'
    else:
        output_field_name = '_dv'
    if kmeans_clusters > -1:
        output_field_name = output_field_name + 'k'
    else:
        output_field_name = output_field_name + 'p'
    str_date = None
    if '_' in input_field_name:
        str_values = input_field_name.split('_')
        for i in range(len(str_values)):
            str_value = str_values[i]
            is_date = True
            if len(str_value) == 6:
                str_value = '20' + str_value
            try:
                date = datetime.datetime.strptime(str_value, date_format)
            except ValueError as error:
                is_date = False
            if is_date:
                str_date = str(date.strftime('%Y')[2:4]) + str(date.strftime('%m')) + str(date.strftime('%d'))
                break
    if str_date:
        output_field_name = str_date + output_field_name
    output_field_id_index = in_layer_definition.GetFieldIndex(output_field_name)
    if output_field_id_index == -1:
        in_layer.CreateField(ogr.FieldDefn(output_field_name, ogr.OFTInteger))#ogr.OFTReal))
    cont_feature = 0
    input_values = []
    position_in_input_values_by_feature_position = {}
    for feature in in_layer:
        value = feature.GetFieldAsDouble(input_field_id_index)
        if value >= crop_minimum_value:
            input_value = {}
            input_value['position'] = cont_feature
            input_value['value'] = value
            input_values.append(input_value)
            position_in_input_values_by_feature_position[cont_feature] = len(input_values) - 1
        cont_feature = cont_feature + 1
    number_of_features = in_layer.GetFeatureCount()
    if kmeans_clusters > -1:
        input_values_cv = numpy.zeros([len(input_values), 1], dtype=numpy.float32)
        cont_feature_crop = 0
        for input_value in input_values:
            input_values_cv[cont_feature_crop][0] = input_value['value']
            cont_feature_crop = cont_feature_crop + 1
        # Define criteria = ( type, max_iter = 10 , epsilon = 1.0 )
        # criteria = (cv.TERM_CRITERIA_EPS + cv.TERM_CRITERIA_MAX_ITER, 10, 1.0)
        criteria = (cv.TERM_CRITERIA_MAX_ITER, 100, 1.0)
        flags = cv.KMEANS_RANDOM_CENTERS
        compactness, labels, centers = cv.kmeans(input_values_cv, kmeans_clusters,
                                                 None, criteria, 10, flags)
        pos_center_min_value = -1
        center_min_value = 100000000.
        for i in range(6):
            if centers[i] < center_min_value:
                center_min_value = centers[i]
                pos_center_min_value = i
        cont_feature = 0
        for feature in in_layer:
            damaged = 0
            if not cont_feature in position_in_input_values_by_feature_position:
                damaged = -1
            else:
                pos_in_input_values = position_in_input_values_by_feature_position[cont_feature]
                pos_center = labels[position_in_input_values_by_feature_position[cont_feature]][0]
                if pos_center == pos_center_min_value:
                    damaged = 1
            cont_feature = cont_feature + 1
            feature.SetField(output_field_name, damaged)
            in_layer.SetFeature(feature)
    elif percentile_minimum_threshold > 0.:
        input_values.sort(key=sortFunction)
        damage_positions = []
        number_of_damages = 0
        threshold_value = -1
        for i in range(0, len(input_values)):
            damage_positions.append(input_values[i]['position'])
            if number_of_damages / number_of_features > percentile_minimum_threshold:
                threshold_value = input_values[i]['value']
                break
            number_of_damages = number_of_damages + 1
        cont_feature = 0
        for feature in in_layer:
            damaged = 0
            if not cont_feature in position_in_input_values_by_feature_position:
                damaged = -1
            else:
                if cont_feature in damage_positions:
                    damaged = 1
            cont_feature = cont_feature + 1
            feature.SetField(output_field_name, damaged)
            in_layer.SetFeature(feature)
    in_vec_ds = None
    return str_error

def process_ndvi(input_shp,
                 kmeans_clusters,
                 percentile_minimum_threshold,
                 input_orthomosaic,
                 blue_band_position,
                 green_band_position,
                 red_band_position,
                 nir_band_position,
                 factor_to_reflectance,
                 soil_maximum_rgb_reflectance,
                 crop_minimum_value,
                 str_date):
    str_error = None
    if kmeans_clusters < 0 and percentile_minimum_threshold < 0:
        str_error = "Function process"
        str_error += "\nkmeans_clusters or percentile_minimum_threshold must be greather than 0"
        return str_error
    elif kmeans_clusters > 0 and percentile_minimum_threshold > 0:
        str_error = "Function process"
        str_error += "\nkmeans_clusters or percentile_minimum_threshold must be greather than 0"
        return str_error
    if not exists(input_shp):
        str_error = "Function process"
        str_error += "\nNot exists file: {}".format(input_shp)
        return str_error
    if not exists(input_orthomosaic):
        str_error = "Function process"
        str_error += "\nNot exists file: {}".format(input_orthomosaic)
        return str_error
    orthomosaic_ds = None
    try:
        orthomosaic_ds = gdal.Open(input_orthomosaic)
    except ValueError:
        str_error = "Function process"
        str_error += "\nError opening dataset file:\n{}".format(input_orthomosaic)
        return str_error
    # orthomosaic_number_of_bands = orthomosaic_ds.GetRasterCount()
    # if blue_band_position > orthomosaic_number_of_bands:
    #     str_error = "Function process"
    #     str_error += "\nBlue band position is greather than orthomosaic number of bands"
    #     return str_error
    orthomosaic_ds_rb_blue = None
    try:
        orthomosaic_ds_rb_blue = orthomosaic_ds.GetRasterBand(blue_band_position)
    except ValueError:
        str_error = "Function process"
        str_error += "\nError getting BLUE raster band from file:\n{}".format(input_orthomosaic)
        return str_error
    orthomosaic_ds_rb_green = None
    try:
        orthomosaic_ds_rb_green = orthomosaic_ds.GetRasterBand(green_band_position)
    except ValueError:
        str_error = "Function process"
        str_error += "\nError getting GREEN raster band from file:\n{}".format(input_orthomosaic)
        return str_error
    orthomosaic_ds_rb_red = None
    try:
        orthomosaic_ds_rb_red = orthomosaic_ds.GetRasterBand(red_band_position)
    except ValueError:
        str_error = "Function process"
        str_error += "\nError getting RED raster band from file:\n{}".format(input_orthomosaic)
        return str_error
    orthomosaic_ds_rb_nir = None
    try:
        orthomosaic_ds_rb_nir = orthomosaic_ds.GetRasterBand(nir_band_position)
    except ValueError:
        str_error = "Function process"
        str_error += "\nError getting NIR raster band from file:\n{}".format(input_orthomosaic)
        return str_error
    orthomosaic_geotransform = orthomosaic_ds.GetGeoTransform()
    orthomosaic_crs = osr.SpatialReference()
    orthomosaic_crs.ImportFromWkt(orthomosaic_ds.GetProjectionRef())
    orthomosaic_crs_wkt = orthomosaic_crs.ExportToWkt()
    ulx, xres, xskew, uly, yskew, yres = orthomosaic_ds.GetGeoTransform()
    lrx = ulx + (orthomosaic_ds.RasterXSize * xres)
    lry = uly + (orthomosaic_ds.RasterYSize * yres)
    out_ring = ogr.Geometry(ogr.wkbLinearRing)
    out_ring.AddPoint(ulx, uly)
    out_ring.AddPoint(lrx, uly)
    out_ring.AddPoint(lrx, lry)
    out_ring.AddPoint(ulx, lry)
    out_ring.AddPoint(ulx, uly)
    orthomosaic_poly = ogr.Geometry(ogr.wkbPolygon)
    orthomosaic_poly.AddGeometry(out_ring)
    rs_pixel_width = orthomosaic_geotransform[1]
    rs_pixel_height = orthomosaic_geotransform[5]
    orthomosaic_pixel_area = abs(rs_pixel_width) * abs(rs_pixel_height)
    orthomosaic_x_origin = orthomosaic_geotransform[0]
    orthomosaic_y_origin = orthomosaic_geotransform[3]
    orthomosaic_pixel_width = orthomosaic_geotransform[1]
    orthomosaic_pixel_height = orthomosaic_geotransform[5]
    driver = ogr.GetDriverByName('ESRI Shapefile')
    in_vec_ds = None
    try:
        in_vec_ds = driver.Open(input_shp, 1)  # 0 means read-only. 1 means writeable.
    except ValueError:
        str_error = "Function process"
        str_error += "\nError opening dataset file:\n{}".format(input_shp)
        return str_error
    in_layer = in_vec_ds.GetLayer()
    in_crs = in_layer.GetSpatialRef()
    in_crs_wkt = in_crs.ExportToWkt()
    in_geometry_type = in_layer.GetGeomType()
    if in_geometry_type != ogr.wkbPolygon \
            and in_geometry_type != ogr.wkbMultiPolygon \
            and in_geometry_type != ogr.wkbPolygonM and in_geometry_type != ogr.wkbPolygonZM:
        str_error = "Function process"
        str_error += "\nNot Polygon geometry type in file:\n{}".format(input_shp)
        return str_error
    in_layer_definition = in_layer.GetLayerDefn()
    number_of_features = in_layer.GetFeatureCount()
    cont_feature = 0
    input_values = []
    position_in_input_values_by_feature_position = {}
    output_field_name = str_date
    output_field_name = output_field_name + '_' + 'dn'
    if kmeans_clusters > -1:
        output_field_name = output_field_name + 'k'
    else:
        output_field_name = output_field_name + 'p'
    output_field_id_index = in_layer_definition.GetFieldIndex(output_field_name)
    if output_field_id_index == -1:
        in_layer.CreateField(ogr.FieldDefn(output_field_name, ogr.OFTInteger))#ogr.OFTReal))
    for feature in in_layer:
        print('Processing plant: {}, of {}'.format(str(cont_feature + 1),
                                                   str(number_of_features)))
        plot_geometry_full = feature.GetGeometryRef().Clone()
        crs_transform = None
        if in_crs_wkt != orthomosaic_crs_wkt:
            crs_transform = osr.CoordinateTransformation(in_crs, orthomosaic_crs)
        if crs_transform:
            plot_geometry_full.Transform(crs_transform)
        plot_geometry = None
        if orthomosaic_poly.Overlaps(plot_geometry_full):
            plot_geometry = plot_geometry_full.Intersection(orthomosaic_poly)
        if orthomosaic_poly.Contains(plot_geometry_full):
            plot_geometry = plot_geometry_full
        if orthomosaic_poly.Within(plot_geometry_full):
            plot_geometry = orthomosaic_poly
        if not plot_geometry:
            continue
        plot_geometry = plot_geometry_full.Intersection(orthomosaic_poly)
        plot_geometry_area = plot_geometry.GetArea()
        if plot_geometry_area < (3 * orthomosaic_pixel_area):
            continue
        geom_points_x = []
        geom_points_y = []
        geom_type_name = plot_geometry.GetGeometryName().lower()
        if "multipolygon" in geom_type_name:
            for i in range(0, plot_geometry.GetGeometryCount()):
                ring = plot_geometry.GetGeometryRef(i).GetGeometryRef(0)
                numpoints = ring.GetPointCount()
                for p in range(numpoints):
                    fc, sc, tc = ring.GetPoint(p)
                    geom_points_x.append(fc)
                    geom_points_y.append(sc)
        elif "polygon" in geom_type_name:
            ring = plot_geometry.GetGeometryRef(0)
            numpoints = ring.GetPointCount()
            for p in range(numpoints):
                fc, sc, tc = ring.GetPoint(p)
                geom_points_x.append(fc)
                geom_points_y.append(sc)
        else:
            # sys.exit("ERROR: Geometry needs to be either Polygon or Multipolygon")
            continue
        plot_geom_x_min = min(geom_points_x)
        plot_geom_x_max = max(geom_points_x)
        plot_geom_y_min = min(geom_points_y)
        plot_geom_y_max = max(geom_points_y)
        # Specify offset and rows and columns to read
        rs_x_off = int((plot_geom_x_min - orthomosaic_x_origin) / rs_pixel_width)
        rs_y_off = int((orthomosaic_y_origin - plot_geom_y_max) / rs_pixel_width)
        x_ul = orthomosaic_x_origin + rs_x_off * rs_pixel_width
        y_ul = orthomosaic_y_origin - rs_y_off * rs_pixel_width
        rs_x_count = int((plot_geom_x_max - plot_geom_x_min) / rs_pixel_width) + 1
        rs_y_count = int((plot_geom_y_max - plot_geom_y_min) / rs_pixel_width) + 1
        # Create memory target raster
        target_orthomosaic = gdal.GetDriverByName('MEM').Create('', rs_x_count, rs_y_count, 1, gdal.GDT_Byte)
        target_orthomosaic.SetGeoTransform((
            plot_geom_x_min, rs_pixel_width, 0,
            plot_geom_y_max, 0, rs_pixel_height,
        ))
        # Create for target raster the same projection as for the value raster
        raster_srs = osr.SpatialReference()
        raster_srs.ImportFromWkt(orthomosaic_ds.GetProjectionRef())
        target_orthomosaic.SetProjection(raster_srs.ExportToWkt())
        target_orthomosaic.SetProjection(raster_srs.ExportToWkt())
        feature_drv = ogr.GetDriverByName('ESRI Shapefile')
        feature_ds= feature_drv.CreateDataSource("/vsimem/memory_name.shp")
        # geometryType = plot_geometry.getGeometryType()
        feature_layer = feature_ds.CreateLayer("layer", orthomosaic_crs, geom_type=plot_geometry.GetGeometryType())
        featureDefnHeaders = feature_layer.GetLayerDefn()
        out_feature = ogr.Feature(featureDefnHeaders)
        out_feature.SetGeometry(plot_geometry)
        feature_layer.CreateFeature(out_feature)
        feature_ds.FlushCache()
        # Rasterize zone polygon to raster blue
        gdal.RasterizeLayer(target_orthomosaic, [1], feature_layer, burn_values=[1])
        feature_orthomosaic_band_mask = target_orthomosaic.GetRasterBand(1)
        feature_orthomosaic_data_mask = (feature_orthomosaic_band_mask.ReadAsArray(0, 0,
                                                                                  rs_x_count, rs_y_count)
                                         .astype(float))
        # Mask zone of raster blue
        feature_orthomosaic_data_blue = (orthomosaic_ds_rb_blue.ReadAsArray(rs_x_off, rs_y_off, rs_x_count, rs_y_count)
                                    .astype(float))
        feature_dsm_raster_array_blue = numpy.ma.masked_array(feature_orthomosaic_data_blue,
                                                         numpy.logical_not(feature_orthomosaic_data_mask))
        orthomosaic_first_indexes_blue, orthomosaic_second_indexes_blue = feature_dsm_raster_array_blue.nonzero()
        # Mask zone of raster green
        feature_orthomosaic_data_green = (orthomosaic_ds_rb_green.ReadAsArray(rs_x_off, rs_y_off, rs_x_count, rs_y_count)
                                    .astype(float))
        feature_dsm_raster_array_green = numpy.ma.masked_array(feature_orthomosaic_data_green,
                                                         numpy.logical_not(feature_orthomosaic_data_mask))
        orthomosaic_first_indexes_green, orthomosaic_second_indexes_green = feature_dsm_raster_array_green.nonzero()
        # Mask zone of raster red
        feature_orthomosaic_data_red = (orthomosaic_ds_rb_red.ReadAsArray(rs_x_off, rs_y_off, rs_x_count, rs_y_count)
                                    .astype(float))
        feature_dsm_raster_array_red = numpy.ma.masked_array(feature_orthomosaic_data_red,
                                                         numpy.logical_not(feature_orthomosaic_data_mask))
        orthomosaic_first_indexes_red, orthomosaic_second_indexes_red = feature_dsm_raster_array_red.nonzero()
        # Mask zone of raster nir
        feature_orthomosaic_data_nir = (orthomosaic_ds_rb_nir.ReadAsArray(rs_x_off, rs_y_off, rs_x_count, rs_y_count)
                                    .astype(float))
        feature_dsm_raster_array_nir = numpy.ma.masked_array(feature_orthomosaic_data_nir,
                                                         numpy.logical_not(feature_orthomosaic_data_mask))
        orthomosaic_first_indexes_nir, orthomosaic_second_indexes_nir = feature_dsm_raster_array_nir.nonzero()
        ndvi_mean = 0.
        ndvi_min = 2.
        ndvi_max = -2.
        ndvi_number_of_values = 0
        for i in range(len(orthomosaic_first_indexes_blue)):
            fi = orthomosaic_first_indexes_blue[i]
            si = orthomosaic_second_indexes_blue[i]
            if (not fi in orthomosaic_first_indexes_green
                    or not fi in orthomosaic_first_indexes_red
                    or not fi in orthomosaic_first_indexes_nir):
                continue
            if (not si in orthomosaic_second_indexes_green
                    or not si in orthomosaic_second_indexes_red
                    or not si in orthomosaic_second_indexes_nir):
                continue
            blue = feature_dsm_raster_array_blue[fi][si] * factor_to_reflectance
            green = feature_dsm_raster_array_green[fi][si] * factor_to_reflectance
            red = feature_dsm_raster_array_red[fi][si] * factor_to_reflectance
            nir = feature_dsm_raster_array_nir[fi][si] * factor_to_reflectance
            if (blue < soil_maximum_rgb_reflectance
                    and green < soil_maximum_rgb_reflectance
                    and red < soil_maximum_rgb_reflectance):
                continue
            ndvi = (nir - red) / (red + nir)
            if ndvi < crop_minimum_value:
                continue
            ndvi_mean = ndvi_mean + ndvi
            if ndvi < ndvi_min:
                ndvi_min = ndvi
            if ndvi > ndvi_max:
                ndvi_max = ndvi
            ndvi_number_of_values = ndvi_number_of_values + 1
        if ndvi_number_of_values > 0:
            ndvi_mean = ndvi_mean / ndvi_number_of_values
            input_value = {}
            input_value['position'] = cont_feature
            input_value['value'] = ndvi_mean
            input_values.append(input_value)
            position_in_input_values_by_feature_position[cont_feature] = len(input_values) - 1
        else:
            ndvi_mean = -1
        cont_feature = cont_feature + 1
    if kmeans_clusters > -1:
        input_values_cv = numpy.zeros([len(input_values), 1], dtype=numpy.float32)
        cont_feature_crop = 0
        for input_value in input_values:
            input_values_cv[cont_feature_crop][0] = input_value['value']
            cont_feature_crop = cont_feature_crop + 1
        # Define criteria = ( type, max_iter = 10 , epsilon = 1.0 )
        # criteria = (cv.TERM_CRITERIA_EPS + cv.TERM_CRITERIA_MAX_ITER, 10, 1.0)
        criteria = (cv.TERM_CRITERIA_MAX_ITER, 100, 1.0)
        flags = cv.KMEANS_RANDOM_CENTERS
        compactness, labels, centers = cv.kmeans(input_values_cv, kmeans_clusters,
                                                 None, criteria, 10, flags)
        pos_center_min_value = -1
        center_min_value = 100000000.
        for i in range(6):
            if centers[i] < center_min_value:
                center_min_value = centers[i]
                pos_center_min_value = i
        cont_feature = 0
        for feature in in_layer:
            damaged = 0
            if not cont_feature in position_in_input_values_by_feature_position:
                damaged = -1
            else:
                pos_in_input_values = position_in_input_values_by_feature_position[cont_feature]
                pos_center = labels[position_in_input_values_by_feature_position[cont_feature]][0]
                if pos_center == pos_center_min_value:
                    damaged = 1
            cont_feature = cont_feature + 1
            feature.SetField(output_field_name, damaged)
            in_layer.SetFeature(feature)
    else:
        input_values.sort(key=sortFunction)
        damage_positions = []
        number_of_damages = 0
        threshold_value = -1
        for i in range(0, len(input_values)):
            damage_positions.append(input_values[i]['position'])
            if number_of_damages / number_of_features > percentile_minimum_threshold:
                threshold_value = input_values[i]['value']
                break
            number_of_damages = number_of_damages + 1
        cont_feature = 0
        for feature in in_layer:
            damaged = 0
            if not cont_feature in position_in_input_values_by_feature_position:
                damaged = -1
            else:
                if cont_feature in damage_positions:
                    damaged = 1
            cont_feature = cont_feature + 1
            feature.SetField(output_field_name, damaged)
            in_layer.SetFeature(feature)
    in_vec_ds = None
    return str_error

def main():
    # ==================
    # parse command line
    # ==================
    usage = "usage: %prog [options] "
    parser = OptionParser(usage=usage)
    parser.add_option("--input_crops_frames_shp", dest="input_crops_frames_shp", action="store", type="string",
                      help="Crops frames shapefile", default=None)
    parser.add_option("--method_data", dest="method_data", action="store", type="string",
                      help="Method input data type: gcc, vol or ndvi", default=None)
    parser.add_option("--method_segmentation", dest="method_segmentation", action="store", type="string",
                      help="Method segmentation: kmeans or percentile", default=None)
    parser.add_option("--crop_minimum_value", dest="crop_minimum_value", action="store", type="string",
                      help="Crop minimum value: gcc (per unit), vol (dm3) or NDVI (per unit)", default=None)
    parser.add_option("--kmeans_clusters", dest="kmeans_clusters", action="store", type="string",
                      help="Number of cluster for kmeans segmentation", default=None)
    parser.add_option("--percentile_minimum_threshold", dest="percentile_minimum_threshold", action="store",
                      type="string", help="Minimum value (per unit) for percentile segmentation", default=None)
    parser.add_option("--input_orthomosaic", dest="input_orthomosaic", action="store", type="string",
                      help="Input multispectral orthomosaic, for ndvi method", default=None)
    parser.add_option("--blue_band_position", dest="blue_band_position", action="store", type="int",
                      help="Blue band position in multispectral orthomosaic, for ndvi method", default=None)
    parser.add_option("--green_band_position", dest="green_band_position", action="store", type="int",
                      help="Green band position in multispectral orthomosaic, for ndvi method", default=None)
    parser.add_option("--red_band_position", dest="red_band_position", action="store", type="int",
                      help="Red band position in multispectral orthomosaic, for ndvi method", default=None)
    parser.add_option("--nir_band_position", dest="nir_band_position", action="store", type="int",
                      help="Nir band position in multispectral orthomosaic, for ndvi method", default=None)
    parser.add_option("--factor_to_reflectance", dest="factor_to_reflectance", action="store",
                      type="string", help="Factor for convert band values to reflectance (per unit), for ndvi method",
                      default=None)
    parser.add_option("--soil_maximum_rgb_reflectance", dest="soil_maximum_rgb_reflectance", action="store",
                      type="string", help="Soil maximum RGB reflectance (per unit) from orthomosaic, for ndvi method",
                      default=None)
    parser.add_option("--output_shp", dest="output_shp", action="store", type="string",
                      help="Output shapefile or none for use input shapefile", default=None)
    parser.add_option("--date_from_orthomosaic_file", dest="date_from_orthomosaic_file", action="store",
                      type="int", help="Read date from orthomosaic file name: 1-yes, 0-No, for ndvi method",
                      default=None)
    parser.add_option("--orthomosaic_file_string_separator", dest="orthomosaic_file_string_separator",
                      action="store", type="string", help="Orthomosaic file string separator, for ndvi method",
                      default=None)
    parser.add_option("--orthomosaic_file_date_string_position", dest="orthomosaic_file_date_string_position",
                      action="store", type="int",
                      help="Orthomosaic file date string position, for ndvi method", default=None)
    parser.add_option("--date_format", dest="date_format", action="store", type="string",
                      help="Date format (%Y%m%d, ...)", default=None)
    parser.add_option("--date", dest="date", action="store", type="string",
                      help="None or date value no from orthomosaic files, for ndvi method", default=None)
    parser.add_option("--input_field", dest="input_field", action="store", type="string",
                      help="Input field name of GCC or VOL, for gcc or vol methods", default=None)
    (options, args) = parser.parse_args()
    if not options.input_crops_frames_shp:
        parser.print_help()
        return
    input_shp = options.input_crops_frames_shp
    if not exists(input_shp):
        print("Error:\nInput crops shapefile does not exists:\n{}".format(input_shp))
        return
    if not options.output_shp:
        parser.print_help()
        return
    use_input_shp = True
    output_shp = options.output_shp
    if output_shp != 'none':
        use_input_shp = False
    if not options.method_segmentation:
        parser.print_help()
        return
    method_segmentation = options.method_segmentation
    kmeans_clusters = -1
    percentile_minimum_threshold = -1.
    if method_segmentation == 'kmeans':
        if not options.kmeans_clusters:
            parser.print_help()
            return
        str_kmeans_clusters = options.kmeans_clusters
        if not is_number(str_kmeans_clusters):
            print("Error:\nInvalid value for kmeans clusters: {}".
                  format(str_kmeans_clusters))
            return
        kmeans_clusters = int(str_kmeans_clusters)
        if kmeans_clusters < 2 or kmeans_clusters > 20:
            print("Error:\nInvalid value for kmeans clusters: {}".
                  format(str_kmeans_clusters))
            return
    elif method_segmentation == 'percentile':
        if not options.percentile_minimum_threshold:
            parser.print_help()
            return
        str_percentile_minimum_threshold = options.percentile_minimum_threshold
        if not is_number(str_percentile_minimum_threshold):
            print("Error:\nInvalid value for percentile minimum threshold: {}".
                  format(str_percentile_minimum_threshold))
            return
        percentile_minimum_threshold = float(str_percentile_minimum_threshold)
        if percentile_minimum_threshold < 0 or percentile_minimum_threshold > 1:
            print("Error:\nInvalid value for percentile minimum threshold: {}".
                  format(str_percentile_minimum_threshold))
            return
    if kmeans_clusters < 0 and percentile_minimum_threshold < 0:
        print("Error:\nMethod segmentation must be: kmeans or percentile")
        return
    if not options.method_data:
        parser.print_help()
        return
    method_data = options.method_data
    ndvi = False
    gcc = False
    vol = False
    if method_data == 'ndvi':
        ndvi = True
    elif method_data == 'gcc':
        gcc = True
    elif method_data == 'vol':
        vol = True
    if not ndvi and not gcc and not vol:
        print("Error:\nMethod data must be: ndvi, gcc or vol")
        return
    if not options.crop_minimum_value:
        parser.print_help()
        return
    str_crop_minimum_value = options.crop_minimum_value
    if not is_number(str_crop_minimum_value):
        print("Error:\nInvalid value for crop minimum value: {}".
              format(str_crop_minimum_value))
        return
    crop_minimum_value = float(str_crop_minimum_value)
    if ndvi:
        if crop_minimum_value < 0 or crop_minimum_value > 1:
            print("Error:\nInvalid value for crop minimum NDVI value: {}".
                  format(str_crop_minimum_value))
            return
    elif gcc:
        if crop_minimum_value < 0 or crop_minimum_value > 10:
            print("Error:\nInvalid value for crop minimum GCC value: {}".
                  format(str_crop_minimum_value))
            return
    elif vol:
        if crop_minimum_value < 0 or crop_minimum_value > 10000:
            print("Error:\nInvalid value for crop minimum VOL value: {}".
                  format(str_crop_minimum_value))
            return
    if not options.date_format:
        parser.print_help()
        return
    date_format = options.date_format.strip()
    input_field = None
    if not use_input_shp:
        str_error = copy_shapefile(input_shp, output_shp)
        if str_error:
            print("Error:\n{}".format(str_error))
            return
        input_shp = output_shp
    if gcc or vol:
        if not options.input_field:
            parser.print_help()
            return
        input_field = options.input_field
        str_error = process_gcc_or_vol(input_shp,
                                       crop_minimum_value,
                                       kmeans_clusters,
                                       percentile_minimum_threshold,
                                       gcc,
                                       vol,
                                       input_field,
                                       date_format)
        if str_error:
            print("Error:\n{}".format(str_error))
            return
        print("... Process finished")
    elif ndvi:
        if not options.input_orthomosaic:
            parser.print_help()
            return
        input_orthomosaic = options.input_orthomosaic
        if not exists(input_orthomosaic):
            print("Error:\nInput orthomosaic does not exists:\n{}".format(input_orthomosaic))
            return
        if not options.blue_band_position:
            parser.print_help()
            return
        str_blue_band_position = options.blue_band_position
        if not is_number(str_blue_band_position):
            print("Error:\nInvalid value for blue band position: {}".
                  format(str_blue_band_position))
            return
        blue_band_position = int(str_blue_band_position)
        if blue_band_position < 1 or blue_band_position > 8:
            print("Error:\nInvalid value for blue band position: {}".
                  format(str_blue_band_position))
            return
        str_green_band_position = options.green_band_position
        if not is_number(str_green_band_position):
            print("Error:\nInvalid value for green band position: {}".
                  format(str_green_band_position))
            return
        green_band_position = int(str_green_band_position)
        if green_band_position < 1 or green_band_position > 8:
            print("Error:\nInvalid value for green band position: {}".
                  format(str_green_band_position))
            return
        if green_band_position == blue_band_position:
            print("Error:\nBand positions must be different")
            return
        str_red_band_position = options.red_band_position
        if not is_number(str_red_band_position):
            print("Error:\nInvalid value for green band position: {}".
                  format(str_red_band_position))
            return
        red_band_position = int(str_red_band_position)
        if red_band_position < 1 or red_band_position > 8:
            print("Error:\nInvalid value for green band position: {}".
                  format(str_red_band_position))
            return
        if red_band_position == blue_band_position or red_band_position == green_band_position:
            print("Error:\nBand positions must be different")
            return
        str_nir_band_position = options.nir_band_position
        if not is_number(str_nir_band_position):
            print("Error:\nInvalid value for green band position: {}".
                  format(str_nir_band_position))
            return
        nir_band_position = int(str_nir_band_position)
        if nir_band_position < 1 or nir_band_position > 8:
            print("Error:\nInvalid value for green band position: {}".
                  format(str_nir_band_position))
            return
        if nir_band_position == blue_band_position or nir_band_position == green_band_position \
                or nir_band_position == red_band_position:
            print("Error:\nBand positions must be different")
            return
        if not options.soil_maximum_rgb_reflectance:
            parser.print_help()
            return
        str_factor_to_reflectance = options.factor_to_reflectance
        if not is_number(str_factor_to_reflectance):
            print("Error:\nInvalid value for factor to reflectance: {}".
                  format(str_factor_to_reflectance))
            return
        factor_to_reflectance = float(str_factor_to_reflectance)
        if factor_to_reflectance < 0 or factor_to_reflectance > 10000.:
            print("Error:\nInvalid value for factor to reflectance: {}".
                  format(str_factor_to_reflectance))
            return
        str_soil_maximum_rgb_reflectance = options.soil_maximum_rgb_reflectance
        if not is_number(str_soil_maximum_rgb_reflectance):
            print("Error:\nInvalid value for percentile minimum threshold: {}".
                  format(str_soil_maximum_rgb_reflectance))
            return
        soil_maximum_rgb_reflectance = float(str_soil_maximum_rgb_reflectance)
        if soil_maximum_rgb_reflectance < 0 or soil_maximum_rgb_reflectance > 1:
            print("Error:\nInvalid value for soil maximum RGB reflectance: {}".
                  format(str_soil_maximum_rgb_reflectance))
            return
        if options.date_from_orthomosaic_file == None:
            parser.print_help()
            return
        date_from_orthomosaic = False
        if options.date_from_orthomosaic_file == 1:
            date_from_orthomosaic = True
        date = None
        if not date_from_orthomosaic:
            if options.date == 'none':
                print("Error:\nDate must be a value if not read from orthomosaic file name")
                return
            str_date = options.date
            is_date = True
            if len(options.date) == 6:
                str_date = '20' + str_date
            try:
                date = datetime.datetime.strptime(str_date, date_format)
            except ValueError as error:
                is_date = False
            if not is_date:
                print("Error:\nInvalid string date from orthomosaic name: {} and format: {}".
                      format(options.date, date_format))
                return
        else:
            if not options.orthomosaic_file_string_separator:
                parser.print_help()
                return
            orthomosaic_file_string_separator = options.orthomosaic_file_string_separator
            if not options.orthomosaic_file_date_string_position:
                parser.print_help()
                return
            orthomosaic_file_date_string_position = options.orthomosaic_file_date_string_position
            orthomosaic_file_name_without_path = os.path.splitext(os.path.basename(input_orthomosaic))[0]
            orthomosaic_file_name_values = orthomosaic_file_name_without_path.split(orthomosaic_file_string_separator)
            if (orthomosaic_file_date_string_position < 0
                    or orthomosaic_file_date_string_position > len(orthomosaic_file_name_values)):
                print("Error:\nInvalid value for orthomosaic files date string position: {}".
                      format(str(orthomosaic_file_date_string_position)))
                return
            str_date = orthomosaic_file_name_values[orthomosaic_file_date_string_position - 1]
            is_date = True
            if len(str_date) == 6:
                str_date = '20' + str_date
            try:
                date = datetime.datetime.strptime(str_date, date_format)
            except ValueError as error:
                is_date = False
            if not is_date:
                print("Error:\nInvalid string date from orthomosaic name: {} and format: {}".
                      format(orthomosaic_file_name_values[orthomosaic_file_date_string_position - 1],
                             date_format))
                return
        str_date = str(date.strftime('%Y')[2:4]) + str(date.strftime('%m')) + str(date.strftime('%d'))
        str_error = process_ndvi(input_shp,
                                 kmeans_clusters,
                                 percentile_minimum_threshold,
                                 input_orthomosaic,
                                 blue_band_position,
                                 green_band_position,
                                 red_band_position,
                                 nir_band_position,
                                 factor_to_reflectance,
                                 soil_maximum_rgb_reflectance,
                                 crop_minimum_value,
                                 str_date)
        if str_error:
            print("Error:\n{}".format(str_error))
            return
        print("... Process finished")


if __name__ == '__main__':
    main()
