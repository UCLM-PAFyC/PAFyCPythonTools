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

# class OptionParser(optparse.OptionParser):
#     def check_required(self, opt):
#         option = self.get_option(opt)
#         # Assumes the option's 'default' is set to None!
#         if getattr(self.values, option.dest) is None:
#             self.error("%s option not supplied" % option)


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
            factor_to_reflectance,
            bands_to_use,
            red_band_number,
            nir_band_number,
            minimum_ndvi,
            maximum_ndvi,
            minimum_explained_variance,
            max_number_of_kmeans_clusters):
    str_error = None
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
    orthomosaic_number_of_bands = orthomosaic_ds.RasterCount
    if red_band_number > orthomosaic_number_of_bands:
        str_error = "Function process"
        str_error += "\nRed band number is greather than orthomosaic number of bands"
        return str_error
    if nir_band_number > orthomosaic_number_of_bands:
        str_error = "Function process"
        str_error += "\nNir band number is greather than orthomosaic number of bands"
        return str_error
    values_by_band = {}
    xSize =orthomosaic_ds.RasterXSize
    ySize = orthomosaic_ds.RasterYSize
    geotransform = orthomosaic_ds.GetGeoTransform()
    projection = orthomosaic_ds.GetProjection()
    for band_number in bands_to_use:
        if band_number > orthomosaic_number_of_bands or band_number < 1:
            str_error = "Function process"
            str_error += ("\nBand to use number: {} is out of valid number for orthomosaic number of bands"
                          .format(str(band_number)))
            return str_error
        band = orthomosaic_ds.GetRasterBand(band_number)
        values_by_band[band_number] = band.ReadAsArray() # rows, columns
        yo = 1
    if not red_band_number in values_by_band:
        band = orthomosaic_ds.GetRasterBand(red_band_number)
        values_by_band[red_band_number] = band.ReadAsArray() # rows, columns
    if not nir_band_number in values_by_band:
        band = orthomosaic_ds.GetRasterBand(nir_band_number)
        values_by_band[nir_band_number] = band.ReadAsArray() # rows, columns
    rows, columns = values_by_band[red_band_number].shape
    valid_indexes = []
    fileformat = "GTiff"
    driver = gdal.GetDriverByName(fileformat)
    output_path = os.path.dirname(os.path.abspath(input_orthomosaic))
    var_path = Path(input_orthomosaic)
    base_name = var_path.stem
    file_ext = '.tif'
    dst_filename_ndvi = (output_path + '/' + base_name + '_ndvi' + file_ext)
    dst_filename_ndvi = os.path.normpath(dst_filename_ndvi)
    dst_ds_ndvi = driver.Create(dst_filename_ndvi, xsize=xSize, ysize=ySize, bands=1, eType=gdal.GDT_Float32)
    dst_ds_ndvi.SetGeoTransform(geotransform)
    dst_ds_ndvi.SetProjection(projection)
    np_raster_ndvi = np.zeros((ySize, xSize), dtype=np.float32)
    dst_filename_ndvi_mask = (output_path + '/' + base_name + '_ndvi_mask' + file_ext)
    dst_filename_ndvi_mask = os.path.normpath(dst_filename_ndvi_mask)
    dst_ds_ndvi_mask = driver.Create(dst_filename_ndvi_mask, xsize=xSize, ysize=ySize, bands=1, eType=gdal.GDT_Byte)
    dst_ds_ndvi_mask.SetGeoTransform(geotransform)
    dst_ds_ndvi_mask.SetProjection(projection)
    np_raster_ndvi_mask = np.zeros((ySize, xSize), dtype=np.uint8)
    # for row in range(100):#range(rows):
    for row in range(rows):
        # for column in range(100):#range(columns):
        for column in range(columns):
            red_value = values_by_band[red_band_number][row][column] * factor_to_reflectance
            nir_value = values_by_band[nir_band_number][row][column] * factor_to_reflectance
            ndvi_value = (nir_value - red_value) / (nir_value + red_value)
            if ndvi_value >= minimum_ndvi and ndvi_value <= maximum_ndvi:
                valid_index = [row, column]
                valid_indexes.append(valid_index)
                np_raster_ndvi_mask[row][column] = 1
            np_raster_ndvi[row][column] = ndvi_value
    # dst_ds_ndvi.GetRasterBand(1).SetNoDataValue(0)
    dst_ds_ndvi.GetRasterBand(1).WriteArray(np_raster_ndvi)
    dst_ds_ndvi_mask.GetRasterBand(1).WriteArray(np_raster_ndvi_mask)
    dst_ds_ndvi = None
    dst_ds_ndvi_mask = None
    data = np.zeros((len(valid_indexes), len(bands_to_use)))
    for i in range(len(valid_indexes)):
        row = valid_indexes[i][0]
        column = valid_indexes[i][1]
        for j in range(len(bands_to_use)):
            data[i][j] = values_by_band[bands_to_use[j]][row][column] * factor_to_reflectance
    values_by_band = None
    standardized_data = (data - data.mean(axis=0)) / data.std(axis=0)
    data = None
    covariance_matrix = np.cov(standardized_data, ddof=1, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
    # np.argsort can only provide lowest to highest; use [::-1] to reverse the list
    order_of_importance = np.argsort(eigenvalues)[::-1]
    # utilize the sort order to sort eigenvalues and eigenvectors
    sorted_eigenvalues = eigenvalues[order_of_importance]
    sorted_eigenvectors = eigenvectors[:, order_of_importance]  # sort the columns
    # use sorted_eigenvalues to ensure the explained variances correspond to the eigenvectors
    explained_variance = sorted_eigenvalues / np.sum(sorted_eigenvalues)
    explained_variance_by_selected_components = 0
    number_of_pca_components = 0
    while explained_variance_by_selected_components < minimum_explained_variance:
        explained_variance_by_selected_components = (explained_variance_by_selected_components
                                                     + explained_variance[number_of_pca_components])
        number_of_pca_components = number_of_pca_components + 1
    # number_of_pca_components = 2
    reduced_data = np.matmul(standardized_data, sorted_eigenvectors[:, :number_of_pca_components])  # transform the original data
    standardized_data = None
    criteria = (cv.TERM_CRITERIA_MAX_ITER, 100, 1.0)
    flags = cv.KMEANS_RANDOM_CENTERS
    input_values_cv = np.zeros([len(valid_indexes), number_of_pca_components], dtype=np.float32)
    for i in range(len(valid_indexes)):
        for j in range(number_of_pca_components):
            input_values_cv[i][j] = reduced_data[i][j]
    mse = {}
    for kmeans_clusters in range(1, max_number_of_kmeans_clusters + 1):
        compactness, labels, centers = cv.kmeans(input_values_cv, kmeans_clusters,
                                                 None, criteria, 10, flags)
        mse[kmeans_clusters] = compactness / len(valid_indexes)
        str_mse = str(round(mse[kmeans_clusters] * 100))
        dst_filename = (output_path + '/' + base_name + '_' + 'npc_'
                        + str(number_of_pca_components) + '_nckms_' + str(kmeans_clusters)
                        + '_mse_' + str_mse + file_ext)
        dst_filename = os.path.normpath(dst_filename)
        dst_ds = driver.Create(dst_filename, xsize=xSize, ysize=ySize, bands=1, eType=gdal.GDT_Byte)
        dst_ds.SetGeoTransform(geotransform)
        dst_ds.SetProjection(projection)
        np_raster = np.zeros((ySize, xSize), dtype=np.uint8)
        for i in range(len(valid_indexes)):
            row = valid_indexes[i][0]
            column = valid_indexes[i][1]
            np_raster[row][column] = labels[i] + 1
        dst_ds.GetRasterBand(1).SetNoDataValue(0)
        dst_ds.GetRasterBand(1).WriteArray(np_raster)
        dst_ds = None

    return str_error


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_orthomosaic", help="Input orthomosaic",type=str)
    parser.add_argument("--factor_to_reflectance", type=float,
                        help="Multiplicative factor for convert raster values to reflectance")
    parser.add_argument("--bands_to_use", nargs="+", type=int, help="Bands to use, starting 1")
    parser.add_argument("--red_band_number", type=int, help="Red band number, starting 1")
    parser.add_argument("--nir_band_number", type=int, help="Nir band number, starting 1")
    parser.add_argument("--minimum_ndvi", type=float, help="Minimmum NDVI, in range [-1,-1]")
    parser.add_argument("--maximum_ndvi", type=float, help="Minimmum NDVI, in range [-1,-1]")
    parser.add_argument("--minimum_explained_variance", type=float,
                        help="Minimmum explained variance by PCA components, in range [0,1]")
    parser.add_argument("--max_number_of_kmeans_clusters", type=int,
                        help="Maximum number of clusters in Kmeans classification process")
    args = parser.parse_args()
    if not args.input_orthomosaic:
        parser.print_help()
        return
    input_orthomosaic = args.input_orthomosaic
    if not exists(input_orthomosaic):
        print("Error:\nInput orthomosaic does not exists:\n{}".format(input_orthomosaic))
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
    max_number_of_kmeans_clusters = args.max_number_of_kmeans_clusters

    str_error = process(input_orthomosaic,
                        factor_to_reflectance,
                        bands_to_use,
                        red_band_number,
                        nir_band_number,
                        minimum_ndvi,
                        maximum_ndvi,
                        minimum_explained_variance,
                        max_number_of_kmeans_clusters)
    if str_error:
        print("Error:\n{}".format(str_error))
        return
    print("... Process finished")


if __name__ == '__main__':
    main()

