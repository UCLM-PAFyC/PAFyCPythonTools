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
            input_rois_shp,
            factor_to_reflectance,
            bands_to_use,
            red_band_number,
            nir_band_number,
            minimum_ndvi,
            minimum_nir_reflectance,
            grid_spacing,
            minimum_explained_variance,
            only_one_principal_component,
            weight_factor_by_cluster,
            input_dsm,
            input_dtm,
            crop_minimum_height,
            output_path):
    str_error = None
    str_error = None
    if input_rois_shp:
        if not exists(input_rois_shp):
            str_error = "Function process"
            str_error += "\nNot exists file: {}".format(input_rois_shp)
            return str_error
    if not exists(input_orthomosaic):
        str_error = "Function process"
        str_error += "\nNot exists file: {}".format(input_orthomosaic)
        return str_error
    dsm_ds = None
    dsm_rb = None
    dsm_crs = None
    dsm_gsd_x = None
    dsm_gsd_y = None
    dsm_xSize = None
    dsm_ySize = None
    dsm_ulx = None
    dsm_uly = None
    dsm_lrx = None
    dsm_lry = None
    dsm_geotransform = None
    dsm_data = None
    dtm_ds = None
    dtm_rb = None
    dtm_crs = None
    dtm_gsd_x = None
    dtm_gsd_y = None
    dtm_xSize = None
    dtm_ySize = None
    dtm_ulx = None
    dtm_uly = None
    dtm_lrx = None
    dtm_lry = None
    dtm_geotransform = None
    dtm_data = None
    if crop_minimum_height > 0.001 or crop_minimum_height < -0.001:
        if not exists(input_dsm):
            str_error = "Function process"
            str_error += "\nNot exists file: {}".format(input_dsm)
            return str_error
        if not exists(input_dtm):
            str_error = "Function process"
            str_error += "\nNot exists file: {}".format(input_dtm)
            return str_error
        try:
            dsm_ds = gdal.Open(input_dsm)
        except ValueError:
            str_error = "Function process"
            str_error += "\nError opening dataset file:\n{}".format(input_dsm)
            return str_error
        try:
            dsm_rb = dsm_ds.GetRasterBand(1)
        except ValueError:
            str_error = "Function process"
            str_error += "\nError getting raster band from file:\n{}".format(input_dsm)
            return str_error
        try:
            dtm_ds = gdal.Open(input_dtm)
        except ValueError:
            str_error = "Function process"
            str_error += "\nError opening dataset file:\n{}".format(input_dtm)
            return str_error
        try:
            dtm_rb = dtm_ds.GetRasterBand(1)
        except ValueError:
            str_error = "Function process"
            str_error += "\nError getting raster band from file:\n{}".format(input_dtm)
            return str_error
        dsm_crs = osr.SpatialReference()
        dsm_crs.ImportFromWkt(dsm_ds.GetProjectionRef())
        dsm_crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        dsm_crs_wkt = dsm_crs.ExportToWkt()
        dsm_xSize = dsm_ds.RasterXSize
        dsm_ySize = dsm_ds.RasterYSize
        dsm_geotransform = dsm_ds.GetGeoTransform()
        dsm_gsd_x = abs(dsm_geotransform[1])
        dsm_gsd_y = abs(dsm_geotransform[5])
        dsm_ulx, dsm_xres, dsm_xskew, dsm_uly, dsm_yskew, dsm_yres = dsm_ds.GetGeoTransform()
        dsm_lrx = dsm_ulx + (dsm_ds.RasterXSize * dsm_xres)
        dsm_lry = dsm_uly + (dsm_ds.RasterYSize * dsm_yres)
        dsm_data = dsm_rb.ReadAsArray()
        dtm_crs = osr.SpatialReference()
        dtm_crs.ImportFromWkt(dtm_ds.GetProjectionRef())
        dtm_crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        dtm_crs_wkt = dtm_crs.ExportToWkt()
        dtm_xSize = dtm_ds.RasterXSize
        dtm_ySize = dtm_ds.RasterYSize
        dtm_geotransform = dtm_ds.GetGeoTransform()
        dtm_gsd_x = abs(dtm_geotransform[1])
        dtm_gsd_y = abs(dtm_geotransform[5])
        dtm_ulx, dtm_xres, dtm_xskew, dtm_uly, dtm_yskew, dtm_yres = dtm_ds.GetGeoTransform()
        dtm_lrx = dtm_ulx + (dtm_ds.RasterXSize * dtm_xres)
        dtm_lry = dtm_uly + (dtm_ds.RasterYSize * dtm_yres)
        dtm_data = dtm_rb.ReadAsArray()
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
    for band_number in bands_to_use:
        if band_number > orthomosaic_number_of_bands or band_number < 1:
            str_error = "Function process"
            str_error += ("\nBand to use number: {} is out of valid number for orthomosaic number of bands"
                          .format(str(band_number)))
            return str_error
    xSize = orthomosaic_ds.RasterXSize
    ySize = orthomosaic_ds.RasterYSize
    orthomosaic_geotransform = orthomosaic_ds.GetGeoTransform()
    gsd_x = abs(orthomosaic_geotransform[1])
    gsd_y = abs(orthomosaic_geotransform[5])
    gsd = gsd_x
    projection = orthomosaic_ds.GetProjection()
    orthomosaic_crs = osr.SpatialReference()
    transform_orthomosaic_to_dsm = None
    transform_orthomosaic_to_dtm = None
    orthomosaic_crs.ImportFromWkt(orthomosaic_ds.GetProjectionRef())
    orthomosaic_crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    if crop_minimum_height > 0.0:
        transform_orthomosaic_to_dsm = osr.CoordinateTransformation(orthomosaic_crs, dsm_crs)
        transform_orthomosaic_to_dtm = osr.CoordinateTransformation(orthomosaic_crs, dtm_crs)
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
    band_red = orthomosaic_ds.GetRasterBand(red_band_number)
    band_red_no_data_value = band_red.GetNoDataValue()
    band_red_data = band_red.ReadAsArray()
    band_redmasked_data = np.ma.masked_where(band_red_data == band_red_no_data_value, band_red_data)
    band_redindexes_without_no_data_value = band_redmasked_data.nonzero()
    columns_by_row = {}
    for i in range(len(band_redindexes_without_no_data_value[0])):
        row = band_redindexes_without_no_data_value[0][i]
        column = band_redindexes_without_no_data_value[1][i]
        if not row in columns_by_row:
            columns_by_row[row] = []
        columns_by_row[row].append(column)
    band_redmasked_data = None
    band_redindexes_without_no_data_value = None
    band_nir = orthomosaic_ds.GetRasterBand(nir_band_number)
    band_nir_no_data_value = band_nir.GetNoDataValue()
    band_nir_data = band_nir.ReadAsArray()
    valid_indexes = []
    invalid_indexes = []
    if crop_minimum_height > 0.0:
        print('Filtering by NIR Reflectance, NDVI value and crop height...', flush=True)
    else:
        print('Filtering by NIR Reflectance and NDVI value ...', flush=True)
    ndvi_valid_values_by_row = {}
    for row in columns_by_row:
        for j in range(len(columns_by_row[row])):
            column = columns_by_row[row][j]
            index = [row, column]
            nir_value = band_nir_data[row][column] * factor_to_reflectance
            if nir_value < minimum_nir_reflectance:
                invalid_indexes.append(index)
                continue
            red_value = band_red_data[row][column] * factor_to_reflectance
            ndvi_value = (nir_value - red_value) / (nir_value + red_value)
            if ndvi_value < minimum_ndvi:
                invalid_indexes.append(index)
                continue
            if crop_minimum_height > 0.001 or crop_minimum_height < -0.001:
                x_coord = column * gsd_x + ulx + (gsd_x / 2.)  # add half the cell size
                y_coord = uly - row * gsd_y - (gsd_y / 2.)  # to centre the point
                x_coord_dsm = x_coord
                y_coord_dsm = y_coord
                x_coord_dsm, y_coord_dsm, _ = transform_orthomosaic_to_dsm.TransformPoint(x_coord_dsm, y_coord_dsm)
                dsm_column = math.floor((x_coord_dsm - dsm_ulx) / dsm_gsd_x)
                dsm_row = math.floor((dsm_uly - y_coord_dsm) / dsm_gsd_y)
                if dsm_column < 0 or dsm_column > dsm_xSize:
                    invalid_indexes.append(index)
                    continue
                if dsm_row < 0 or dsm_row > dsm_ySize:
                    invalid_indexes.append(index)
                    continue
                dsm_height = dsm_data[dsm_row][dsm_column]
                x_coord_dtm = x_coord
                y_coord_dtm = y_coord
                x_coord_dtm, y_coord_dtm, _ = transform_orthomosaic_to_dtm.TransformPoint(x_coord_dtm, y_coord_dtm)
                dtm_column = math.floor((x_coord_dtm - dtm_ulx) / dtm_gsd_x)
                dtm_row = math.floor((dtm_uly - y_coord_dtm) / dtm_gsd_y)
                if dtm_column < 0 or dtm_column > dtm_xSize:
                    invalid_indexes.append(index)
                    continue
                if dtm_row < 0 or dtm_row > dtm_ySize:
                    continue
                dtm_height = dtm_data[dtm_row][dtm_column]
                crop_height = dsm_height - dtm_height
                if crop_minimum_height > 0.001 and crop_height < crop_minimum_height:
                    invalid_indexes.append(index)
                    continue
                elif crop_minimum_height < 0.001 and crop_height > (-1. * crop_minimum_height):
                    invalid_indexes.append(index)
                    continue
            valid_indexes.append(index)
            if not row in ndvi_valid_values_by_row:
                ndvi_valid_values_by_row[row] = {}
            ndvi_valid_values_by_row[row][column] = ndvi_value
    print('   ... Process finished', flush=True)
    if crop_minimum_height > 0.0:
        dsm_data = None
        dtm_data = None
    data = np.zeros((len(valid_indexes), len(bands_to_use)))
    invalid_positions_in_valid_indexes = []
    for j in range(len(bands_to_use)):
        band_number = bands_to_use[j]
        print('Getting data for band {} ...'.format(band_number), flush=True)
        band_data = None
        band_no_data_value = None
        if band_number == red_band_number:
            band_data = band_red_data
            band_no_data_value = band_red_no_data_value
        elif band_number == nir_band_number:
            band_data = band_nir_data
            band_no_data_value = band_nir_no_data_value
        else:
            band = orthomosaic_ds.GetRasterBand(band_number)
            band_no_data_value = band.GetNoDataValue()
            band_data = band.ReadAsArray()  # rows, columns
        for i in range(len(valid_indexes)):
            if i in invalid_positions_in_valid_indexes:
                continue
            row = valid_indexes[i][0]
            column = valid_indexes[i][1]
            value = band_data[row][column]
            if value == band_no_data_value:
                index = [row, column]
                invalid_positions_in_valid_indexes.append(i)
            data[i][j] = value * factor_to_reflectance
        print('   ... Process finished', flush=True)
    band_red_data = None
    band_nir_data = None
    number_of_kmeans_clusters = len(weight_factor_by_cluster)
    if len(invalid_positions_in_valid_indexes) > 0:
        print('Removing pixels with no data value in some band ...', flush=True)
        valid_indexes_with_out_outliers = []
        data_without_outliers = np.zeros((len(valid_indexes) - len(invalid_positions_in_valid_indexes),
                                          len(bands_to_use)))
        pos = -1
        for i in range(len(valid_indexes)):
            if i in invalid_positions_in_valid_indexes:
                continue
            pos = pos + 1
            row = valid_indexes[i][0]
            column = valid_indexes[i][1]
            for j in range(len(bands_to_use)):
                data_without_outliers[pos][j] = data[i][j]
            index = [row, column]
            valid_indexes_with_out_outliers.append(index)
        data = None
        valid_indexes = None
        data = data_without_outliers
        valid_indexes = valid_indexes_with_out_outliers
        data_without_outliers = None
        valid_indexes_with_out_outliers = None
        print('   ... Process finished', flush=True)
    print('Computing principal components and transforming data values to new base ...', flush=True)
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
    print('   ... Process finished', flush=True)
    print('Computing results for {} clusters ...'.format(str(number_of_kmeans_clusters)), flush=True)
    compactness, labels, centers = cv.kmeans(input_values_cv, number_of_kmeans_clusters,
                                             None, criteria, 10, flags)
    ndvi_by_cluster = {}
    for i in range(len(valid_indexes)):
        row = valid_indexes[i][0]
        column = valid_indexes[i][1]
        n_cluster = labels[i].item() + 1
        if row in ndvi_valid_values_by_row:
            if column in ndvi_valid_values_by_row[row]:
                if not n_cluster in ndvi_by_cluster:
                    ndvi_by_cluster[n_cluster] = {}
                    ndvi_by_cluster[n_cluster]['values'] = []
                    ndvi_by_cluster[n_cluster]['values_rows'] = []
                    ndvi_by_cluster[n_cluster]['values_columns'] = []
                    ndvi_by_cluster[n_cluster]['mean'] = 0.
                    ndvi_by_cluster[n_cluster]['std'] = 0.
                ndvi_value = ndvi_valid_values_by_row[row][column]
                ndvi_by_cluster[n_cluster]['values'].append(ndvi_value)
                ndvi_by_cluster[n_cluster]['values_rows'].append(row)
                ndvi_by_cluster[n_cluster]['values_columns'].append(column)
                ndvi_by_cluster[n_cluster]['mean'] = ndvi_by_cluster[n_cluster]['mean'] + ndvi_value
    output_grid_by_column_by_row = {}
    for n_cluster in ndvi_by_cluster:
        number_of_values = len(ndvi_by_cluster[n_cluster]['values'])
        ndvi_by_cluster[n_cluster]['mean'] = ndvi_by_cluster[n_cluster]['mean'] / float(number_of_values)
        # if number_of_values > 1:
        ndvi_mean_value = ndvi_by_cluster[n_cluster]['mean']
        for i_ndvi_value in range(len(ndvi_by_cluster[n_cluster]['values'])):
            ndvi_value = ndvi_by_cluster[n_cluster]['values'][i_ndvi_value]
            ndvi_column = ndvi_by_cluster[n_cluster]['values_columns'][i_ndvi_value]
            ndvi_row = ndvi_by_cluster[n_cluster]['values_rows'][i_ndvi_value]
            x_coord = ndvi_column * gsd_x + ulx + (gsd_x / 2.)  # add half the cell size
            y_coord = uly - ndvi_row * gsd_y - (gsd_y / 2.)  # to centre the point
            grid_column = math.floor((x_coord - math.floor(ulx)) / grid_spacing)
            # grid_row = math.floor((x_coord - math.ceil(uly)) / grid_spacing)
            grid_row = math.floor((math.ceil(uly) - y_coord) / grid_spacing)
            if not grid_column in output_grid_by_column_by_row:
                output_grid_by_column_by_row[grid_column] = {}
            if not grid_row in output_grid_by_column_by_row[grid_column]:
                output_grid_by_column_by_row[grid_column][grid_row] = {}
            if not n_cluster in output_grid_by_column_by_row[grid_column][grid_row]:
                output_grid_by_column_by_row[grid_column][grid_row][n_cluster] = 0
            output_grid_by_column_by_row[grid_column][grid_row][n_cluster] = (
                    output_grid_by_column_by_row[grid_column][grid_row][n_cluster] + 1)
            ndvi_diff_value = ndvi_value - ndvi_mean_value
            ndvi_by_cluster[n_cluster]['std'] = ndvi_by_cluster[n_cluster]['std'] + ndvi_diff_value ** 2.
        ndvi_by_cluster[n_cluster]['std'] = math.sqrt(ndvi_by_cluster[n_cluster]['std'] / number_of_values)
    cluster_descending_order_position = []
    position_in_vector_descending_order_by_cluster = {}
    for n_cluster in ndvi_by_cluster:
        max_mean_ndvi_value = 0.
        n_cluster_max_mean_ndvi_value = -1
        for n_cluster_bis in ndvi_by_cluster:
            if n_cluster_bis in cluster_descending_order_position:
                continue
            if ndvi_by_cluster[n_cluster_bis]['mean'] > max_mean_ndvi_value:
                max_mean_ndvi_value = ndvi_by_cluster[n_cluster_bis]['mean']
                n_cluster_max_mean_ndvi_value = n_cluster_bis
        position_in_vector_descending_order_by_cluster[n_cluster_max_mean_ndvi_value] \
            = len(cluster_descending_order_position)
        cluster_descending_order_position.append(n_cluster_max_mean_ndvi_value)
    fileformat = "GTiff"
    driver = gdal.GetDriverByName(fileformat)
    if not output_path:
        output_path = os.path.dirname(os.path.abspath(input_orthomosaic))
    var_path = Path(input_orthomosaic)
    base_name = var_path.stem
    file_ext = '.tif'
    for i_n_cluster in range(len(cluster_descending_order_position)):
        n_cluster = cluster_descending_order_position[i_n_cluster]
        str_mean = str(int(round(ndvi_by_cluster[n_cluster]['mean'] * 100)))
        str_std = str(int(round(ndvi_by_cluster[n_cluster]['std'] * 100)))
        dst_filename = (output_path + '/' + base_name + '_nckms_' + str(i_n_cluster+1)
                        + '_mean_' + str_mean + '_std_' + str_std + file_ext)
        dst_filename = os.path.normpath(dst_filename)
        dst_ds = driver.Create(dst_filename, xsize=xSize, ysize=ySize, bands=1, eType=gdal.GDT_Byte)
        dst_ds.SetGeoTransform(orthomosaic_geotransform)
        dst_ds.SetProjection(projection)
        # np_raster = np.zeros((ySize, xSize), dtype=np.uint8)
        np_raster = np.full((ySize, xSize), 255, dtype=np.uint8)
        for i_ndvi_value in range(len(ndvi_by_cluster[n_cluster]['values'])):
            ndvi_value = ndvi_by_cluster[n_cluster]['values'][i_ndvi_value]
            ndvi_row = ndvi_by_cluster[n_cluster]['values_rows'][i_ndvi_value]
            ndvi_column = ndvi_by_cluster[n_cluster]['values_columns'][i_ndvi_value]
            np_raster[ndvi_row][ndvi_column] = int(round(ndvi_value * 100))
        dst_ds.GetRasterBand(1).SetNoDataValue(0)
        dst_ds.GetRasterBand(1).WriteArray(np_raster)
        dst_ds = None
    outShapefile = (output_path + '/' + base_name + '_grid.shp')
    outShapefile = os.path.normpath(outShapefile)
    outDriver = ogr.GetDriverByName("ESRI Shapefile")
    # Remove output shapefile if it already exists
    if os.path.exists(outShapefile):
        outDriver.DeleteDataSource(outShapefile)
    # Create the output shapefile
    outDataSource = outDriver.CreateDataSource(outShapefile)
    outLayer = outDataSource.CreateLayer("grid", orthomosaic_crs, ogr.wkbPolygon)
    idField = ogr.FieldDefn("id", ogr.OFTInteger)
    outLayer.CreateField(idField)
    fractionCoverField = ogr.FieldDefn("frac_cover", ogr.OFTReal)
    fractionCoverField.SetPrecision(2)
    outLayer.CreateField(fractionCoverField)
    indexField = ogr.FieldDefn("index",ogr.OFTReal)
    indexField.SetPrecision(2)
    outLayer.CreateField(indexField)
    for i_n_cluster in range(len(cluster_descending_order_position)):
        cluster_percentage_field_name = "cl_fc_" + str(i_n_cluster + 1)
        cluster_percentage_Field = ogr.FieldDefn(cluster_percentage_field_name, ogr.OFTReal)
        cluster_percentage_Field.SetPrecision(2)
        outLayer.CreateField(cluster_percentage_Field)
        cluster_index_field_name = "cl_idx_" + str(i_n_cluster + 1)
        cluster_index_Field = ogr.FieldDefn(cluster_index_field_name, ogr.OFTReal)
        cluster_index_Field.SetPrecision(2)
        outLayer.CreateField(cluster_index_Field)
    feature_count = 0
    number_of_pixels_in_grid = round((grid_spacing * grid_spacing) / (gsd_x * gsd_y))
    for grid_column in output_grid_by_column_by_row:
        for grid_row in output_grid_by_column_by_row[grid_column]:
            featureDefn = outLayer.GetLayerDefn()
            feature = ogr.Feature(featureDefn)
            ring = ogr.Geometry(ogr.wkbLinearRing)
            grid_ul_x = math.floor(ulx) + grid_column * grid_spacing
            grid_ul_y = math.ceil(uly) - grid_row * grid_spacing
            grid_ur_x = math.floor(ulx) + (grid_column + 1) * grid_spacing
            grid_ur_y = math.ceil(uly) - grid_row * grid_spacing
            grid_lr_x = math.floor(ulx) + (grid_column + 1) * grid_spacing
            grid_lr_y = math.ceil(uly) - (grid_row + 1) * grid_spacing
            grid_ll_x = math.floor(ulx) + grid_column * grid_spacing
            grid_ll_y = math.ceil(uly) - (grid_row + 1) * grid_spacing
            ring.AddPoint(grid_ul_x, grid_ul_y)
            ring.AddPoint(grid_ur_x, grid_ur_y)
            ring.AddPoint(grid_lr_x, grid_lr_y)
            ring.AddPoint(grid_ll_x, grid_ll_y)
            ring.AddPoint(grid_ul_x, grid_ul_y)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            feature.SetGeometry(poly)
            feature_count = feature_count + 1
            feature.SetField("id", feature_count)
            number_of_pixels_in_clusters = 0
            for n_cluster in output_grid_by_column_by_row[grid_column][grid_row]:
                ordered_cluster = position_in_vector_descending_order_by_cluster[n_cluster]
                number_of_pixels_in_clusters = (number_of_pixels_in_clusters
                                                + output_grid_by_column_by_row[grid_column][grid_row][n_cluster])
            # fraction_cover = (number_of_pixels_in_clusters * xSize * ySize) / (grid_spacing * grid_spacing) * 100.
            fraction_cover = ( number_of_pixels_in_clusters / number_of_pixels_in_grid ) * 100.
            if fraction_cover > 100.:
                fraction_cover = 100.
            feature.SetField("frac_cover", fraction_cover)
            index = 0.
            for n_cluster in output_grid_by_column_by_row[grid_column][grid_row]:
                ordered_cluster = position_in_vector_descending_order_by_cluster[n_cluster]
                # percentage_pixels_in_cluster = (output_grid_by_column_by_row[grid_column][n_cluster]
                #                                 / number_of_pixels_in_clusters * 100.)
                percentage_pixels_in_cluster = float((output_grid_by_column_by_row[grid_column][grid_row][n_cluster]
                                                / number_of_pixels_in_grid ) * 100.)
                if percentage_pixels_in_cluster > 100.:
                    percentage_pixels_in_cluster = 100.
                index_in_cluster = float(weight_factor_by_cluster[ordered_cluster] * percentage_pixels_in_cluster)
                index = index + index_in_cluster
                cluster_percentage_field_name = "cl_fc_" + str(ordered_cluster + 1)
                cluster_index_field_name = "cl_idx_" + str(ordered_cluster + 1)
                feature.SetField(cluster_percentage_field_name, percentage_pixels_in_cluster)
                feature.SetField(cluster_index_field_name, index_in_cluster)
            feature.SetField("index", index)
            outLayer.CreateFeature(feature)
            feature = None
    outDataSource = None
    print('   ... Process finished', flush=True)
    return str_error


def clip_raster(input_raster,
                input_shp,
                no_data_value,
                output_path,
                remove_existing):
    str_error = None
    output_raster = ""
    if not exists(input_shp):
        str_error = "Function clip_raster"
        str_error += "\nNot exists file: {}".format(input_shp)
        return str_error, output_raster
    if not exists(input_raster):
        str_error = "Function clip_raster"
        str_error += "\nNot exists file: {}".format(input_raster)
        return str_error, output_raster
    shp_var_path = Path(input_shp)
    shp_base_name = shp_var_path.stem
    output_raster_suffix = ''
    if not output_path:
        output_path = os.path.dirname(os.path.abspath(input_raster))
    raster_base_name = os.path.basename(input_raster).split('.')[0]
    output_raster = output_path + '\\' + raster_base_name
    # output_raster = os.path.splitext(input_raster)[0]
    output_raster = output_raster + "_rois"
    output_raster = output_raster + "_"
    output_raster = output_raster + shp_base_name
    output_raster = output_raster + os.path.splitext(input_raster)[1]
    # return str_error, output_raster
    if exists(output_raster):
        if not remove_existing:
            return str_error, output_raster
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
        output_raster_ds = gdal.Warp(output_raster, input_raster_ds, options=gdalwarp_str_options)
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
    parser.add_argument("--minimum_nir_reflectance", type=float, help="Minimmum NIR reflectance, in range [0,-1]")
    parser.add_argument("--minimum_explained_variance", type=float,
                        help="Minimmum explained variance by PCA components, in range [0,1]")
    parser.add_argument("--only_one_principal_component", type=int,
                        help="Use only one principal compponent. Not if compute from explained variance")
    parser.add_argument("--grid_spacing", type=float,
                        help="Grid spacing, in meters")
    parser.add_argument("--weight_factor_by_cluster", nargs="+", type=float,
                        help="Weight factor by cluster")
    parser.add_argument("--input_dsm", dest="input_dsm", action="store", type=str,
                        help="DSM geotiff, or '' for no use it, crop_minimum_height == 0.0", default=None)
    parser.add_argument("--input_dtm", dest="input_dtm", action="store", type=str,
                        help="DTM geotiff, or '' for no use it, crop_minimum_height == 0.0", default=None)
    parser.add_argument("--crop_minimum_height, or 0.0 for no use it", dest="crop_minimum_height", action="store", type=float,
                        help="Crop minimum height, in meters, positive for vegetation over value and negative for vegetation below value", default=None)
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
    if not args.minimum_nir_reflectance:
        parser.print_help()
        return
    minimum_nir_reflectance = args.minimum_nir_reflectance
    if not args.grid_spacing:
        parser.print_help()
        return
    grid_spacing = args.grid_spacing
    if not args.minimum_explained_variance:
        parser.print_help()
        return
    minimum_explained_variance = args.minimum_explained_variance
    if not args.only_one_principal_component:
        parser.print_help()
        return
    int_only_one_principal_component = args.only_one_principal_component
    if int_only_one_principal_component < 0 or args.only_one_principal_component > 1:
        print("Error:\nParameter only_one_principal_component must be 0 or 1")
    only_one_principal_component = False
    if int_only_one_principal_component == 1:
        only_one_principal_component = True
    if not args.weight_factor_by_cluster:
        parser.print_help()
    weight_factor_by_cluster = args.weight_factor_by_cluster
    if args.crop_minimum_height == None:
        parser.print_help()
        return
    crop_minimum_height = args.crop_minimum_height
    input_dsm = ''
    input_dtm = ''
    if crop_minimum_height > 0.001 or crop_minimum_height < -0.001:
        if not args.input_dsm:
            parser.print_help()
            return
        if not args.input_dtm:
            parser.print_help()
            return
        input_dsm = args.input_dsm
        if not exists(input_dsm):
            print("Error:\nInput DSM does not exists:\n{}".format(input_dsm))
            return
        input_dtm = args.input_dtm
        if not exists(input_dtm):
            print("Error:\nInput DTM does not exists:\n{}".format(input_dtm))
            return
    output_path = args.output_path
    if output_path:
        if not exists(output_path):
            print("Error:\nOutput path does not exists:\n{}".format(output_path))
            return
    if input_rois_shp:
        remove_existing = False
        str_error, input_orthomosaic_rois = clip_raster(input_orthomosaic,
                                                        input_rois_shp,
                                                        no_data_value,
                                                        output_path,
                                                        remove_existing)
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
                        minimum_nir_reflectance,
                        grid_spacing,
                        minimum_explained_variance,
                        only_one_principal_component,
                        weight_factor_by_cluster,
                        input_dsm,
                        input_dtm,
                        crop_minimum_height,
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
