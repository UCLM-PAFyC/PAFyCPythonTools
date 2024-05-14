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
import math


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
            output_file = output_base_path + "/" + output_base_name + file_extension
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
    A = trunc((I - 1867216.25) / 36524.25)
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


def process_cwsith(input_shp,
                   kmeans_clusters,
                   percentile_maximum_threshold,
                   input_orthomosaic,
                   temperature_air,
                   relative_humidity,
                   upper_line_coef_a,
                   upper_line_coef_b,
                   lower_line_coef_a,
                   lower_line_coef_b,
                   factor_to_temperature,
                   str_date):
    str_error = None
    if kmeans_clusters < 0 and percentile_maximum_threshold < 0:
        str_error = "Function process"
        str_error += "\nkmeans_clusters or percentile_maximum_threshold must be greather than 0"
        return str_error
    elif kmeans_clusters > 0 and percentile_maximum_threshold > 0:
        str_error = "Function process"
        str_error += "\nkmeans_clusters or percentile_maximum_threshold must be greather than 0"
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
    # if orthomosaic_number_of_bands != 1:
    #     str_error = "Function process"
    #     str_error += "\nOrthomosaic number of bands must be 1"
    #     return str_error
    orthomosaic_ds_rb = None
    try:
        orthomosaic_ds_rb = orthomosaic_ds.GetRasterBand(1)
    except ValueError:
        str_error = "Function process"
        str_error += "\nError getting BLUE raster band from file:\n{}".format(input_orthomosaic)
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
    output_field_name_tm = str_date
    output_field_name_tm = output_field_name_tm + '_' + 'tm'
    if kmeans_clusters > -1:
        output_field_name_tm = output_field_name_tm + 'k'
    else:
        output_field_name_tm = output_field_name_tm + 'p'
    output_field_tm_id_index = in_layer_definition.GetFieldIndex(output_field_name_tm)
    if output_field_tm_id_index == -1:
        in_layer.CreateField(ogr.FieldDefn(output_field_name_tm, ogr.OFTReal))#ogr.OFTReal))
    output_field_name_ts = str_date
    output_field_name_ts = output_field_name_ts + '_' + 'ts'
    if kmeans_clusters > -1:
        output_field_name_ts = output_field_name_ts + 'k'
    else:
        output_field_name_ts = output_field_name_ts + 'p'
    output_field_ts_id_index = in_layer_definition.GetFieldIndex(output_field_name_ts)
    if output_field_ts_id_index == -1:
        in_layer.CreateField(ogr.FieldDefn(output_field_name_ts, ogr.OFTReal))#ogr.OFTReal))
    output_field_name_nvs = str_date
    output_field_name_nvs = output_field_name_nvs + '_' + 'nvs'
    output_field_nvs_id_index = in_layer_definition.GetFieldIndex(output_field_name_nvs)
    if output_field_nvs_id_index == -1:
        in_layer.CreateField(ogr.FieldDefn(output_field_name_nvs, ogr.OFTInteger))#ogr.OFTReal))
    output_field_name_c = str_date
    output_field_name_c = output_field_name_c + '_' + 'c'
    if kmeans_clusters > -1:
        output_field_name_c = output_field_name_c + 'k'
    else:
        output_field_name_c = output_field_name_c + 'p'
    output_field_c_id_index = in_layer_definition.GetFieldIndex(output_field_name_c)
    if output_field_c_id_index == -1:
        in_layer.CreateField(ogr.FieldDefn(output_field_name_c, ogr.OFTReal))#ogr.OFTReal))
    output_field_name_cx = str_date
    output_field_name_cx = output_field_name_cx + '_' + 'cx'
    if kmeans_clusters > -1:
        output_field_name_cx = output_field_name_cx + 'k'
    else:
        output_field_name_cx = output_field_name_cx + 'p'
    output_field_cx_id_index = in_layer_definition.GetFieldIndex(output_field_name_cx)
    if output_field_cx_id_index == -1:
        in_layer.CreateField(ogr.FieldDefn(output_field_name_cx, ogr.OFTReal))#ogr.OFTReal))
    output_field_name_cn = str_date
    output_field_name_cn = output_field_name_cn + '_' + 'cn'
    if kmeans_clusters > -1:
        output_field_name_cn = output_field_name_cn + 'k'
    else:
        output_field_name_cn = output_field_name_cn + 'p'
    output_field_cn_id_index = in_layer_definition.GetFieldIndex(output_field_name_cn)
    if output_field_cn_id_index == -1:
        in_layer.CreateField(ogr.FieldDefn(output_field_name_cn, ogr.OFTReal))#ogr.OFTReal))
    es = (0.611 * math.exp((17.27 * temperature_air) / (237.3 + temperature_air)))
    ea = es * relative_humidity / 100.
    VPD = es - ea
    dTul = upper_line_coef_a * VPD + upper_line_coef_b
    dTll = lower_line_coef_a * VPD + lower_line_coef_b
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
        feature_orthomosaic_data = (orthomosaic_ds_rb.ReadAsArray(rs_x_off, rs_y_off, rs_x_count, rs_y_count)
                                    .astype(float))
        feature_raster_array = numpy.ma.masked_array(feature_orthomosaic_data,
                                                         numpy.logical_not(feature_orthomosaic_data_mask))
        orthomosaic_first_indexes, orthomosaic_second_indexes = feature_raster_array.nonzero()
        temperature_values = []
        crop_temperature_values = []
        for i in range(len(orthomosaic_first_indexes)):
            fi = orthomosaic_first_indexes[i]
            si = orthomosaic_second_indexes[i]
            temperature = feature_raster_array[fi][si] * factor_to_temperature
            temperature_values.append(temperature)
        if kmeans_clusters > -1:
            input_values_cv = numpy.zeros([len(temperature_values), 1], dtype=numpy.float32)
            cont_value = 0
            for temperature_value in temperature_values:
                input_values_cv[cont_value][0] = temperature_values[cont_value]
                cont_value = cont_value + 1
            # Define criteria = ( type, max_iter = 10 , epsilon = 1.0 )
            # criteria = (cv.TERM_CRITERIA_EPS + cv.TERM_CRITERIA_MAX_ITER, 10, 1.0)
            criteria = (cv.TERM_CRITERIA_MAX_ITER, 100, 1.0)
            flags = cv.KMEANS_RANDOM_CENTERS
            compactness, labels, centers = cv.kmeans(input_values_cv, kmeans_clusters,
                                                     None, criteria, 10, flags)
            pos_center_min_value = -1
            center_min_value = 100000000.
            # for i in range(6):
            for i in range(kmeans_clusters):
                if centers[i] < center_min_value:
                    center_min_value = centers[i]
                    pos_center_min_value = i
            cont_value = 0
            for temperature_value in temperature_values:
                if labels[cont_value] == pos_center_min_value:
                    crop_temperature_values.append(temperature_value)
                cont_value = cont_value + 1
        else:
            temperature_values.sort(key=sortFunction)
            threshold_value = -1
            cont_value = 0
            for i in range(0, len(temperature_values)):
                crop_temperature_values.append(temperature_values[cont_value])
                cont_value = cont_value + 1
                if cont_value / len(temperature_values) > percentile_maximum_threshold:
                    threshold_value = temperature_values[cont_value-1]
                    break
        crop_temperature_mean = 0
        for crop_temperature_value in crop_temperature_values:
            crop_temperature_mean = crop_temperature_mean + crop_temperature_value
        crop_temperature_mean = crop_temperature_mean / len(crop_temperature_values)
        crop_temperature_std = 0
        if len(crop_temperature_values) > 1:
            for crop_temperature_value in crop_temperature_values:
                crop_temperature_std = crop_temperature_std + pow(crop_temperature_value - crop_temperature_mean, 2.)
            crop_temperature_std = sqrt(crop_temperature_std / (len(crop_temperature_values) - 1))
        else:
            crop_temperature_std = -1.
        dTx = crop_temperature_mean - temperature_air
        cwsi = (dTx - dTll) / (dTul - dTll)
        dTx_min = crop_temperature_mean - crop_temperature_std - temperature_air
        cwsi_min = (dTx_min - dTll) / (dTul - dTll)
        dTx_max = crop_temperature_mean + crop_temperature_std - temperature_air
        cwsi_max = (dTx_max - dTll) / (dTul - dTll)
        cont_feature = cont_feature + 1
        feature.SetField(output_field_name_tm, crop_temperature_mean)
        feature.SetField(output_field_name_ts, crop_temperature_std)
        feature.SetField(output_field_name_nvs, len(crop_temperature_values))
        feature.SetField(output_field_name_c, cwsi)
        feature.SetField(output_field_name_cx, cwsi_max)
        feature.SetField(output_field_name_cn, cwsi_min)
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
    parser.add_option("--method_segmentation", dest="method_segmentation", action="store", type="string",
                      help="Method segmentation: kmeans or percentile", default=None)
    parser.add_option("--kmeans_clusters", dest="kmeans_clusters", action="store", type="string",
                      help="Number of cluster for kmeans segmentation", default=None)
    parser.add_option("--percentile_maximum_threshold", dest="percentile_maximum_threshold", action="store",
                      type="string", help="Maximum value (per unit) for percentile segmentation", default=None)
    parser.add_option("--input_orthomosaic", dest="input_orthomosaic", action="store", type="string",
                      help="Input thermal orthomosaic", default=None)
    parser.add_option("--temperature", dest="temperature", action="store", type="string",
                      help="Air temperature (celsius degrees)", default=None)
    parser.add_option("--relative_humidity", dest="relative_humidity", action="store", type="string",
                      help="Relative humidity (percentage)", default=None)
    parser.add_option("--upper_line_coef_a", dest="upper_line_coef_a", action="store", type="string",
                      help="Upper line coefficient A (slope)", default=None)
    parser.add_option("--upper_line_coef_b", dest="upper_line_coef_b", action="store", type="string",
                      help="Upper line coefficient B (y value for x equal to 0)", default=None)
    parser.add_option("--lower_line_coef_a", dest="lower_line_coef_a", action="store", type="string",
                      help="Lower line coefficient A (slope)", default=None)
    parser.add_option("--lower_line_coef_b", dest="lower_line_coef_b", action="store", type="string",
                      help="Lower line coefficient B (y value for x equal to 0)", default=None)
    parser.add_option("--factor_to_temperature", dest="factor_to_temperature", action="store", type="string",
                      help="Factor for convert band values in thermal raster to temperature (celsius degrees)")
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
    percentile_maximum_threshold = -1.
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
        if not options.percentile_maximum_threshold:
            parser.print_help()
            return
        str_percentile_maximum_threshold = options.percentile_maximum_threshold
        if not is_number(str_percentile_maximum_threshold):
            print("Error:\nInvalid value for percentile maximum threshold: {}".
                  format(str_percentile_maximum_threshold))
            return
        percentile_maximum_threshold = float(str_percentile_maximum_threshold)
        if percentile_maximum_threshold < 0 or percentile_maximum_threshold > 1:
            print("Error:\nInvalid value for percentile maximum threshold: {}".
                  format(str_percentile_maximum_threshold))
            return
    if kmeans_clusters < 0 and percentile_maximum_threshold < 0:
        print("Error:\nMethod segmentation must be: kmeans or percentile")
        return
    input_orthomosaic = options.input_orthomosaic
    if not exists(input_orthomosaic):
        print("Error:\nInput orthomosaic does not exists:\n{}".format(input_orthomosaic))
        return
    if not options.temperature:
        parser.print_help()
        return
    str_temperature = options.temperature
    if not is_number(str_temperature):
        print("Error:\nInvalid value for temperature value: {}".
              format(str_temperature))
        return
    temperature_air = float(str_temperature)
    if not options.relative_humidity:
        parser.print_help()
        return
    str_relative_humidity = options.relative_humidity
    if not is_number(str_relative_humidity):
        print("Error:\nInvalid value for relative humidity value: {}".
              format(str_relative_humidity))
        return
    relative_humidity = float(str_relative_humidity)
    if not options.upper_line_coef_a:
        parser.print_help()
        return
    str_upper_line_coef_a = options.upper_line_coef_a
    if not is_number(str_upper_line_coef_a):
        print("Error:\nInvalid value for upper line coefficient A value: {}".
              format(str_upper_line_coef_a))
        return
    upper_line_coef_a = float(str_upper_line_coef_a)
    if not options.upper_line_coef_b:
        parser.print_help()
        return
    str_upper_line_coef_b = options.upper_line_coef_b
    if not is_number(str_upper_line_coef_b):
        print("Error:\nInvalid value for upper line coefficient B value: {}".
              format(str_upper_line_coef_b))
        return
    upper_line_coef_b = float(str_upper_line_coef_b)
    if not options.lower_line_coef_a:
        parser.print_help()
        return
    str_lower_line_coef_a = options.lower_line_coef_a
    if not is_number(str_lower_line_coef_a):
        print("Error:\nInvalid value for lower line coefficient A value: {}".
              format(str_lower_line_coef_a))
        return
    lower_line_coef_a = float(str_lower_line_coef_a)
    if not options.lower_line_coef_b:
        parser.print_help()
        return
    str_lower_line_coef_b = options.lower_line_coef_b
    if not is_number(str_lower_line_coef_b):
        print("Error:\nInvalid value for lower line coefficient B value: {}".
              format(str_lower_line_coef_b))
        return
    lower_line_coef_b = float(str_lower_line_coef_b)
    if not options.factor_to_temperature:
        parser.print_help()
        return
    str_factor_to_temperature = options.factor_to_temperature
    if not is_number(str_factor_to_temperature):
        print("Error:\nInvalid value for factor to temperature: {}".
              format(str_factor_to_temperature))
        return
    factor_to_temperature = float(str_factor_to_temperature)
    if options.date_from_orthomosaic_file == None:
        parser.print_help()
        return
    date_from_orthomosaic = False
    if options.date_from_orthomosaic_file == 1:
        date_from_orthomosaic = True
    date = None
    if not options.date_format:
        parser.print_help()
        return
    date_format = options.date_format.strip()
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
    input_field = None
    if not use_input_shp:
        str_error = copy_shapefile(input_shp, output_shp)
        if str_error:
            print("Error:\n{}".format(str_error))
            return
        input_shp = output_shp
    process_cwsith(input_shp,
                   kmeans_clusters,
                   percentile_maximum_threshold,
                   input_orthomosaic,
                   temperature_air,
                   relative_humidity,
                   upper_line_coef_a,
                   upper_line_coef_b,
                   lower_line_coef_a,
                   lower_line_coef_b,
                   factor_to_temperature,
                   str_date)
    if str_error:
        print("Error:\n{}".format(str_error))
        return
    print("... Process finished")


if __name__ == '__main__':
    main()
