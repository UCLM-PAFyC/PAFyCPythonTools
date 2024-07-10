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
            minimum_classification_area):
    str_error = None
    if input_shp:
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
    orthomosaic_number_of_bands = orthomosaic_ds.RasterCount
    if red_band_number > orthomosaic_number_of_bands:
        str_error = "Function process"
        str_error += "\nRed band number is greather than orthomosaic number of bands"
        return str_error
    if nir_band_number > orthomosaic_number_of_bands:
        str_error = "Function process"
        str_error += "\nNir band number is greather than orthomosaic number of bands"
        return str_error
    xSize = orthomosaic_ds.RasterXSize
    ySize = orthomosaic_ds.RasterYSize
    orthomosaic_geotransform = orthomosaic_ds.GetGeoTransform()
    gsd_x = abs(orthomosaic_geotransform[1])
    gsd_y = abs(orthomosaic_geotransform[5])
    gsd = gsd_x
    median_filter_number_of_pixels = ceil(minimum_classification_area / gsd)
    if median_filter_number_of_pixels % 2 == 0:
        median_filter_number_of_pixels = median_filter_number_of_pixels + 1
    median_filter_position = round((median_filter_number_of_pixels ** 2 + 1) / 2)
    projection = orthomosaic_ds.GetProjection()
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
    values_by_band = {}
    no_data_values_by_band = {}
    for band_number in bands_to_use:
        if band_number > orthomosaic_number_of_bands or band_number < 1:
            str_error = "Function process"
            str_error += ("\nBand to use number: {} is out of valid number for orthomosaic number of bands"
                          .format(str(band_number)))
            return str_error
        band = orthomosaic_ds.GetRasterBand(band_number)
        no_data_value = band.GetNoDataValue()
        no_data_values_by_band[band_number] = no_data_value
        values_by_band[band_number] = band.ReadAsArray()  # rows, columns
    if not red_band_number in values_by_band:
        band = orthomosaic_ds.GetRasterBand(red_band_number)
        values_by_band[red_band_number] = band.ReadAsArray()  # rows, columns
    if not nir_band_number in values_by_band:
        band = orthomosaic_ds.GetRasterBand(nir_band_number)
        values_by_band[nir_band_number] = band.ReadAsArray()  # rows, columns
    valid_indexes = []
    invalid_indexes = []
    if input_shp:
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
        features_indexes = np.full((ySize * xSize, 2), -1, dtype=np.compat.long)
        index_pos = -1
        for feature in in_layer:
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
            feature_ds = feature_drv.CreateDataSource("/vsimem/memory_name.shp")
            # geometryType = plot_geometry.getGeometryType()
            feature_layer = feature_ds.CreateLayer("layer", orthomosaic_crs,
                                                   geom_type=plot_geometry.GetGeometryType())
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
            band = orthomosaic_ds.GetRasterBand(1)
            feature_orthomosaic_data = (band.ReadAsArray(rs_x_off, rs_y_off, rs_x_count, rs_y_count)
                                        .astype(float))
            feature_raster_array = np.ma.masked_array(feature_orthomosaic_data,
                                                      np.logical_not(feature_orthomosaic_data_mask))
            orthomosaic_first_indexes, orthomosaic_second_indexes = feature_raster_array.nonzero()
            for i in range(len(orthomosaic_first_indexes)):
                fi = orthomosaic_first_indexes[i]
                si = orthomosaic_second_indexes[i]
                index_pos = index_pos + 1
                features_indexes[index_pos][0] = fi + rs_y_off
                features_indexes[index_pos][1] = si + rs_x_off
        for i in range(index_pos + 1):
            row = features_indexes[i][0]
            column = features_indexes[i][1]
            red_value = values_by_band[red_band_number][row][column]
            nir_value = values_by_band[nir_band_number][row][column]
            index = [row, column]
            if no_data_values_by_band[red_band_number]:
                if red_value == no_data_values_by_band[red_band_number]:
                    invalid_indexes.append(index)
                    continue
            if no_data_values_by_band[nir_band_number]:
                if nir_value == no_data_values_by_band[nir_band_number]:
                    invalid_indexes.append(index)
                    continue
            red_value = red_value * factor_to_reflectance
            nir_value = nir_value * factor_to_reflectance
            ndvi_value = (nir_value - red_value) / (nir_value + red_value)
            if ndvi_value >= minimum_ndvi and ndvi_value <= maximum_ndvi:
                valid_indexes.append(index)
                # np_raster_ndvi_mask[row][column] = 1
            else:
                invalid_indexes.append(index)
            # np_raster_ndvi[row][column] = ndvi_value
    else:
        rows, columns = values_by_band[red_band_number].shape
        # for row in range(400):#range(rows):
        for row in range(rows):
            # for column in range(400):#range(columns):
            for column in range(columns):
                red_value = values_by_band[red_band_number][row][column]
                nir_value = values_by_band[nir_band_number][row][column]
                index = [row, column]
                if no_data_values_by_band[red_band_number]:
                    if red_value == no_data_values_by_band[red_band_number]:
                        invalid_indexes.append(index)
                        continue
                if no_data_values_by_band[nir_band_number]:
                    if nir_value == no_data_values_by_band[nir_band_number]:
                        invalid_indexes.append(index)
                        continue
                red_value = red_value * factor_to_reflectance
                nir_value = nir_value * factor_to_reflectance
                ndvi_value = (nir_value - red_value) / (nir_value + red_value)
                if ndvi_value >= minimum_ndvi and ndvi_value <= maximum_ndvi:
                    valid_indexes.append(index)
                    # np_raster_ndvi_mask[row][column] = 1
                else:
                    invalid_indexes.append(index)
                # np_raster_ndvi[row][column] = ndvi_value

    # fileformat = "GTiff"
    # driver = gdal.GetDriverByName(fileformat)
    # output_path = os.path.dirname(os.path.abspath(input_orthomosaic))
    # var_path = Path(input_orthomosaic)
    # base_name = var_path.stem
    # file_ext = '.tif'
    # dst_filename_ndvi = (output_path + '/' + base_name + '_ndvi' + file_ext)
    # dst_filename_ndvi = os.path.normpath(dst_filename_ndvi)
    # dst_ds_ndvi = driver.Create(dst_filename_ndvi, xsize=xSize, ysize=ySize, bands=1, eType=gdal.GDT_Float32)
    # dst_ds_ndvi.SetGeoTransform(orthomosaic_geotransform)
    # dst_ds_ndvi.SetProjection(projection)
    # np_raster_ndvi = np.zeros((ySize, xSize), dtype=np.float32)
    # dst_filename_ndvi_mask = (output_path + '/' + base_name + '_ndvi_mask' + file_ext)
    # dst_filename_ndvi_mask = os.path.normpath(dst_filename_ndvi_mask)
    # dst_ds_ndvi_mask = driver.Create(dst_filename_ndvi_mask, xsize=xSize, ysize=ySize, bands=1, eType=gdal.GDT_Byte)
    # dst_ds_ndvi_mask.SetGeoTransform(orthomosaic_geotransform)
    # dst_ds_ndvi_mask.SetProjection(projection)
    # np_raster_ndvi_mask = np.zeros((ySize, xSize), dtype=np.uint8)
    # # for row in range(400):#range(rows):
    # for row in range(rows):
    #     # for column in range(400):#range(columns):
    #     for column in range(columns):
    #         red_value = values_by_band[red_band_number][row][column] * factor_to_reflectance
    #         nir_value = values_by_band[nir_band_number][row][column] * factor_to_reflectance
    #         ndvi_value = (nir_value - red_value) / (nir_value + red_value)
    #         index = [row, column]
    #         if ndvi_value >= minimum_ndvi and ndvi_value <= maximum_ndvi:
    #             valid_indexes.append(index)
    #             np_raster_ndvi_mask[row][column] = 1
    #         else:
    #             invalid_indexes.append(index)
    #         np_raster_ndvi[row][column] = ndvi_value
    # # dst_ds_ndvi.GetRasterBand(1).SetNoDataValue(0)
    # dst_ds_ndvi.GetRasterBand(1).WriteArray(np_raster_ndvi)
    # dst_ds_ndvi_mask.GetRasterBand(1).WriteArray(np_raster_ndvi_mask)
    # dst_ds_ndvi = None
    # dst_ds_ndvi_mask = None

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
    if number_of_pca_components > 1 and only_one_principal_component:
        number_of_pca_components = 1
    reduced_data = np.matmul(standardized_data,
                             sorted_eigenvectors[:, :number_of_pca_components])  # transform the original data
    standardized_data = None
    criteria = (cv.TERM_CRITERIA_MAX_ITER, 100, 1.0)
    flags = cv.KMEANS_RANDOM_CENTERS
    input_values_cv = np.zeros([len(valid_indexes), number_of_pca_components], dtype=np.float32)
    for i in range(len(valid_indexes)):
        for j in range(number_of_pca_components):
            input_values_cv[i][j] = reduced_data[i][j]
    mse = {}
    fileformat = "GTiff"
    driver = gdal.GetDriverByName(fileformat)
    output_path = os.path.dirname(os.path.abspath(input_orthomosaic))
    var_path = Path(input_orthomosaic)
    base_name = var_path.stem
    file_ext = '.tif'
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
        dst_ds.SetGeoTransform(orthomosaic_geotransform)
        dst_ds.SetProjection(projection)
        # np_raster = np.zeros((ySize, xSize), dtype=np.uint8)
        np_raster = np.full((ySize, xSize), 0, dtype=np.uint8)
        for i in range(len(valid_indexes)):
            row = valid_indexes[i][0]
            column = valid_indexes[i][1]
            np_raster[row][column] = labels[i] + 1
        np_raster = cv.medianBlur(np_raster, median_filter_number_of_pixels)
        for i in range(len(invalid_indexes)):
            row = invalid_indexes[i][0]
            column = invalid_indexes[i][1]
            np_raster[row][column] = 0
        dst_ds.GetRasterBand(1).SetNoDataValue(0)
        dst_ds.GetRasterBand(1).WriteArray(np_raster)
        dst_ds = None
    return str_error


def clip_raster(input_raster,
                input_shp,
                no_data_value):
    str_error = None
    output_raster_suffix = ''
    output_raster = os.path.splitext(input_raster)[0]
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
    output_raster_suffix = ''
    output_raster = os.path.splitext(input_raster)[0]
    output_raster = output_raster + "_rois"
    output_raster = output_raster + os.path.splitext(input_raster)[1]
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
    if input_rois_shp:
        str_error, input_orthomosaic_rois = clip_raster(input_orthomosaic,
                                                        input_rois_shp,
                                                        no_data_value)
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
                        minimum_classification_area)
    if str_error:
        print("Error:\n{}".format(str_error))
        return
    print("... Process finished")


if __name__ == '__main__':
    err = GdalErrorHandler()
    gdal.PushErrorHandler(err.handler)
    gdal.UseExceptions()  # Exceptions will get raised on anything >= gdal.CE_Failure
    assert err.err_level == gdal.CE_None, 'the error level starts at 0'
    main()
