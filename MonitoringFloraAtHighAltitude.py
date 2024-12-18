# authors:
# David Hernandez Lopez, david.hernandez@uclm.es
# Miguel Angel Moreno Hidalgo, miguelangel.moreno@uclm.es

# import optparse
import argparse
import numpy as np
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
from pathlib import Path


class GdalErrorHandler(object):
    def __init__(self):
        self.err_level = gdal.CE_None
        self.err_no = 0
        self.err_msg = ''

    def handler(self, err_level, err_no, err_msg):
        self.err_level = err_level
        self.err_no = err_no
        self.err_msg = err_msg


def is_number(n):
    is_number = True
    try:
        num = float(n)
        # check for "nan" floats
        is_number = num == num  # or use `math.isnan(num)`
    except ValueError:
        is_number = False
    return is_number


def process(input_orthomosaic,
            input_shp,
            factor_to_reflectance,
            bands_to_use,
            red_band_number,
            nir_band_number,
            minimum_ndvi,
            maximum_ndvi,
            minimum_explained_variance,
            only_one_principal_component,
            max_number_of_kmeans_clusters,
            minimum_classification_area,
            output_path):
    str_error = None

    return str_error



def clip_raster(input_raster,
                input_shp,
                no_data_value,
                output_path):
    str_error = None
    output_raster_suffix = ''
    if not output_path:
        output_path = os.path.dirname(os.path.abspath(input_raster))
    raster_base_name = os.path.basename(input_raster).split('.')[0]
    output_raster = output_path + '\\' + raster_base_name
    # output_raster = os.path.splitext(input_raster)[0]
    output_raster = output_raster + "_rois"
    output_raster = output_raster + os.path.splitext(input_raster)[1]
    if not exists(input_shp):
        str_error = "Function clip_raster"
        str_error += "\nNot exists file: {}".format(input_shp)
        return str_error, output_raster
    if not exists(input_raster):
        str_error = "Function clip_raster"
        str_error += "\nNot exists file: {}".format(input_raster)
        return str_error, output_raster
    # return str_error, output_raster
    if exists(output_raster):
        try:
            os.remove(output_raster)
        except FileNotFoundError as e:
            str_error = ("Removing output raster:\n{}".format(output_raster))
            str_error = str_error + ("\nError\t" + e.strerror)
            return str_error, output_raster
        except OSError as e:
            str_error = ("Removing output raster:\n{}".format(output_raster))
            str_error = str_error + ("\nError\t" + e.strerror)
            return str_error, output_raster
    try:
        input_raster_ds = gdal.Open(input_raster)
    except Exception as e:
        assert err.err_level == gdal.CE_Failure, (
                'The handler error level should now be at failure')
        assert err.err_msg == e.args[0], 'raised exception should contain the message'
        str_error = "Function clip_raster"
        str_error = ('Handled warning: level={}, no={}, msg={}'.format(
                err.err_level, err.err_no, err.err_msg))
        return str_error, output_raster
    raster_no_data_value = input_raster_ds.GetRasterBand(1).GetNoDataValue()
    if not raster_no_data_value:
        datatype = gdal.GetDataTypeName(input_raster_ds.GetRasterBand(1).DataType)
        raster_no_data_value = no_data_value
    orthomosaic_geotransform = input_raster_ds.GetGeoTransform()
    gsd_x = abs(orthomosaic_geotransform[1])
    gsd_y = abs(orthomosaic_geotransform[5])
    gdalwarp_str_options = " -cutline " + input_shp
    gdalwarp_str_options += " -crop_to_cutline -dstnodata " + str(no_data_value)
    gdalwarp_str_options += " -tr "
    gdalwarp_str_options += "{:.12f}".format(gsd_x)
    gdalwarp_str_options += " "
    gdalwarp_str_options += "{:.12f}".format(gsd_y)
    gdalwarp_str_options += " -co COMPRESS=LZW"
    print('Clipping raster ...', flush=True)
    try:
        output_raster_ds = gdal.Warp(output_raster, input_raster_ds, options = gdalwarp_str_options)
    except Exception as e:
        assert err.err_level == gdal.CE_Failure, (
                'The handler error level should now be at failure')
        assert err.err_msg == e.args[0], 'raised exception should contain the message'
        str_error = "Function clip_raster"
        str_error = ('Handled warning: level={}, no={}, msg={}'.format(
                err.err_level, err.err_no, err.err_msg))
        return str_error, output_raster
    print('   ... Process finished', flush=True)
    return str_error, output_raster


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_orthomosaic", help="Input orthomosaic", type=str)
    parser.add_argument("--no_data_value", type=int,
                        help="Raster no data value, if not defined in file")
    parser.add_argument("--input_rois_shp", help="Input input rois shapefile, if exists", type=str)
    parser.add_argument("--factor_to_reflectance", type=float,
                        help="Multiplicative factor for convert raster values to reflectance")
    parser.add_argument("--bands_to_use", nargs="+", type=int, help="Bands to use, starting 1")
    parser.add_argument("--red_band_number", type=int, help="Red band number, starting 1")
    parser.add_argument("--nir_band_number", type=int, help="Nir band number, starting 1")
    parser.add_argument("--minimum_ndvi", type=float, help="Minimmum NDVI, in range [-1,-1]")
    parser.add_argument("--maximum_ndvi", type=float, help="Minimmum NDVI, in range [-1,-1]")
    parser.add_argument("--minimum_explained_variance", type=float,
                        help="Minimmum explained variance by PCA components, in range [0,1]")
    parser.add_argument("--only_one_principal_component", type=int,
                        help="Use only one principal compponent. Not if compute from explained variance")
    parser.add_argument("--max_number_of_kmeans_clusters", type=int,
                        help="Maximum number of clusters in Kmeans classification process")
    parser.add_argument("--minimum_classification_area", type=float,
                        help="Minimum classification area, in meters")
    parser.add_argument("--output_path", type=str,
                        help="Output path or empty for multispectral orthomosaic path")
    args = parser.parse_args()
    if not args.input_orthomosaic:
        parser.print_help()
        return
    input_orthomosaic = args.input_orthomosaic
    if not exists(input_orthomosaic):
        print("Error:\nInput orthomosaic does not exists:\n{}".format(input_orthomosaic))
        return
    if not args.no_data_value:
        parser.print_help()
    no_data_value = args.no_data_value
    input_rois_shp = None
    if args.input_rois_shp:
        input_rois_shp = args.input_rois_shp
        if input_rois_shp:
            if not exists(input_rois_shp):
                print("Error:\nInput ROIs shapefile does not exists:\n{}".format(input_rois_shp))
                return
    if not args.factor_to_reflectance:
        parser.print_help()
    factor_to_reflectance = args.factor_to_reflectance
    if not args.bands_to_use:
        parser.print_help()
        return
    bands_to_use = args.bands_to_use
    if not args.red_band_number:
        parser.print_help()
        return
    red_band_number = args.red_band_number
    if not args.nir_band_number:
        parser.print_help()
        return
    nir_band_number = args.nir_band_number
    if not args.minimum_ndvi:
        parser.print_help()
        return
    minimum_ndvi = args.minimum_ndvi
    if not args.maximum_ndvi:
        parser.print_help()
        return
    maximum_ndvi = args.maximum_ndvi
    if not args.minimum_explained_variance:
        parser.print_help()
        return
    minimum_explained_variance = args.minimum_explained_variance
    if not args.max_number_of_kmeans_clusters:
        parser.print_help()
        return
    if not args.only_one_principal_component:
        parser.print_help()
        return
    int_only_one_principal_component = args.only_one_principal_component
    if int_only_one_principal_component < 0 or args.only_one_principal_component > 1:
        print("Error:\nParameter only_one_principal_component must be 0 or 1")
    only_one_principal_component = False
    if int_only_one_principal_component == 1:
        only_one_principal_component = True
    max_number_of_kmeans_clusters = args.max_number_of_kmeans_clusters
    if not args.minimum_classification_area:
        parser.print_help()
    minimum_classification_area = args.minimum_classification_area
    output_path = args.output_path
    if output_path:
        if not exists(output_path):
            print("Error:\nOutput path does not exists:\n{}".format(output_path))
            return
    if input_rois_shp:
        str_error, input_orthomosaic_rois = clip_raster(input_orthomosaic,
                                                        input_rois_shp,
                                                        no_data_value,
                                                        output_path)
        if str_error:
            print("Error:\n{}".format(str_error))
            return
        input_orthomosaic = input_orthomosaic_rois
        input_rois_shp = None
    str_error = process(input_orthomosaic,
                        input_rois_shp,
                        factor_to_reflectance,
                        bands_to_use,
                        red_band_number,
                        nir_band_number,
                        minimum_ndvi,
                        maximum_ndvi,
                        minimum_explained_variance,
                        only_one_principal_component,
                        max_number_of_kmeans_clusters,
                        minimum_classification_area,
                        output_path)
    if str_error:
        print("Error:\n{}".format(str_error))
        return
    print("... Process finished", flush=True)


if __name__ == '__main__':
    err = GdalErrorHandler()
    gdal.PushErrorHandler(err.handler)
    gdal.UseExceptions()  # Exceptions will get raised on anything >= gdal.CE_Failure
    assert err.err_level == gdal.CE_None, 'the error level starts at 0'
    main()
