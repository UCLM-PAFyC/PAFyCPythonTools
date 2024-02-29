# authors:
# David Hernandez Lopez, david.hernandez@uclm.es
# Miguel Angel Moreno Hidalgo, miguelangel.moreno@uclm.es
# Diego Guerrero Sevilla, diego.guerrero@uclm.es

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


def process(input_shp,
            input_dsm,
            input_dtm,
            use_input_shp,
            output_shp,
            crop_minimum_height,
            compute_GCC,
            compute_volume,
            str_date):
    str_error = ''
    if not exists(input_dsm):
        str_error = "Function process"
        str_error += "\nNot exists file: {}".format(input_dsm)
        return str_error
    if not exists(input_dtm):
        str_error = "Function process"
        str_error += "\nNot exists file: {}".format(input_dtm)
        return str_error
    if not exists(input_shp):
        str_error = "Function process"
        str_error += "\nNot exists file: {}".format(input_shp)
        return str_error
    dsm_ds = None
    try:
        dsm_ds = gdal.Open(input_dsm)
    except ValueError:
        str_error = "Function process"
        str_error += "\nError opening dataset file:\n{}".format(input_dsm)
        return str_error
    dsm_rb = None
    try:
        dsm_rb = dsm_ds.GetRasterBand(1)
    except ValueError:
        str_error = "Function process"
        str_error += "\nError getting raster band from file:\n{}".format(input_dsm)
        return str_error
    dsm_geotransform = dsm_ds.GetGeoTransform()
    dtm_ds = None
    try:
        dtm_ds = gdal.Open(input_dtm)
    except ValueError:
        str_error = "Function process"
        str_error += "\nError opening dataset file:\n{}".format(input_dtm)
        return str_error
    dtm_rb = None
    try:
        dtm_rb = dtm_ds.GetRasterBand(1)
    except ValueError:
        str_error = "Function process"
        str_error += "\nError getting raster band from file:\n{}".format(input_dtm)
        return str_error
    dtm_geotransform = dtm_ds.GetGeoTransform()
    for i in range(len(dsm_geotransform)):
        if abs(dsm_geotransform[i]-dtm_geotransform[i]) > 1.0:
            str_error = "Function process"
            str_error += "\nDifferent georreferencing in DTM and DSM"
            return str_error
    dsm_crs = osr.SpatialReference()
    dsm_crs.ImportFromWkt(dsm_ds.GetProjectionRef())
    dtm_crs = osr.SpatialReference()
    dtm_crs.ImportFromWkt(dtm_ds.GetProjectionRef())
    dsm_crs_wkt = dsm_crs.ExportToWkt()
    dtm_crs_wkt = dtm_crs.ExportToWkt()
    if dsm_crs_wkt != dtm_crs_wkt:
        str_error = "Function process"
        str_error += "\nDifferent crs in DTM and DSM"
        return str_error
    ulx, xres, xskew, uly, yskew, yres = dsm_ds.GetGeoTransform()
    lrx = ulx + (dsm_ds.RasterXSize * xres)
    lry = uly + (dsm_ds.RasterYSize * yres)
    out_ring = ogr.Geometry(ogr.wkbLinearRing)
    out_ring.AddPoint(ulx, uly)
    out_ring.AddPoint(lrx, uly)
    out_ring.AddPoint(lrx, lry)
    out_ring.AddPoint(ulx, lry)
    out_ring.AddPoint(ulx, uly)
    dems_poly = ogr.Geometry(ogr.wkbPolygon)
    dems_poly.AddGeometry(out_ring)
    rs_pixel_width = dsm_geotransform[1]
    rs_pixel_height = dsm_geotransform[5]
    dems_pixel_area = abs(rs_pixel_width) * abs(rs_pixel_height)
    dems_x_origin = dsm_geotransform[0]
    dems_y_origin = dsm_geotransform[3]
    dems_pixel_width = dsm_geotransform[1]
    dems_pixel_height = dsm_geotransform[5]
    driver = ogr.GetDriverByName('ESRI Shapefile')
    in_vec_ds = None
    try:
        if not use_input_shp:
            in_vec_ds = driver.Open(input_shp, 0)  # 0 means read-only. 1 means writeable.
        else:
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
    out_vec_ds = None
    out_layer = None
    field_name_gcc = str_date + "_gcc"
    field_name_vol = str_date + "_vol"
    in_layer_definition = in_layer.GetLayerDefn()
    if use_input_shp:
        if compute_GCC:
            field_id_index = in_layer_definition.GetFieldIndex(field_name_gcc)
            if field_id_index == -1:
                in_layer.CreateField(ogr.FieldDefn(field_name_gcc, ogr.OFTReal))
        if compute_volume:
            field_id_index = in_layer_definition.GetFieldIndex(field_name_vol)
            if field_id_index == -1:
                in_layer.CreateField(ogr.FieldDefn(field_name_vol, ogr.OFTReal))
    else:
        out_ds = driver.CreateDataSource(output_shp)
        out_layer = out_ds.CreateLayer(output_shp.split(".")[0], in_crs, in_geometry_type)
        for i in range(0, in_layer_definition.GetFieldCount()):
            fieldDefn = in_layer_definition.GetFieldDefn(i)
            out_layer.CreateField(fieldDefn)
        if compute_GCC:
            out_layer.CreateField(ogr.FieldDefn(field_name_gcc, ogr.OFTReal))
        if compute_volume:
            out_layer.CreateField(ogr.FieldDefn(field_name_vol, ogr.OFTReal))
    number_of_features = in_layer.GetFeatureCount()
    cont_feature = 0
    for feature in in_layer:
        cont_feature = cont_feature + 1
        print('Processing plant: {}, of {}'.format(str(cont_feature),
                                                   str(number_of_features)))
        plot_geometry_full = feature.GetGeometryRef().Clone()
        crs_transform = None
        if in_crs_wkt != dsm_crs_wkt:
            crs_transform = osr.CoordinateTransformation(in_crs, dsm_crs)
        if crs_transform:
            plot_geometry_full.Transform(crs_transform)
        plot_geometry = None
        if dems_poly.Overlaps(plot_geometry_full) or dems_poly.Intersects(plot_geometry_full):
            plot_geometry = plot_geometry_full.Intersection(dems_poly)
        if dems_poly.Contains(plot_geometry_full):
            plot_geometry = plot_geometry_full
        if dems_poly.Within(plot_geometry_full):
            plot_geometry = dems_poly
        if not plot_geometry:
            continue
        plot_geometry = plot_geometry_full.Intersection(dems_poly)
        plot_geometry_area = plot_geometry.GetArea()
        if plot_geometry_area < (3 * dems_pixel_area):
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
        rs_x_off = int((plot_geom_x_min - dems_x_origin) / rs_pixel_width)
        rs_y_off = int((dems_y_origin - plot_geom_y_max) / rs_pixel_width)
        x_ul = dems_x_origin + rs_x_off * rs_pixel_width
        y_ul = dems_y_origin - rs_y_off * rs_pixel_width
        rs_x_count = int((plot_geom_x_max - plot_geom_x_min) / rs_pixel_width) + 1
        rs_y_count = int((plot_geom_y_max - plot_geom_y_min) / rs_pixel_width) + 1
        # Create memory target raster
        target_dsm = gdal.GetDriverByName('MEM').Create('', rs_x_count, rs_y_count, 1, gdal.GDT_Byte)
        target_dsm.SetGeoTransform((
            plot_geom_x_min, rs_pixel_width, 0,
            plot_geom_y_max, 0, rs_pixel_height,
        ))
        target_dtm = gdal.GetDriverByName('MEM').Create('', rs_x_count, rs_y_count, 1, gdal.GDT_Byte)
        target_dtm.SetGeoTransform((
            plot_geom_x_min, rs_pixel_width, 0,
            plot_geom_y_max, 0, rs_pixel_height,
        ))
        # Create for target raster the same projection as for the value raster
        raster_srs = osr.SpatialReference()
        raster_srs.ImportFromWkt(dsm_ds.GetProjectionRef())
        target_dsm.SetProjection(raster_srs.ExportToWkt())
        target_dtm.SetProjection(raster_srs.ExportToWkt())
        feature_drv = ogr.GetDriverByName('ESRI Shapefile')
        feature_ds= feature_drv.CreateDataSource("/vsimem/memory_name.shp")
        # geometryType = plot_geometry.getGeometryType()
        feature_layer = feature_ds.CreateLayer("layer", dsm_crs, geom_type=plot_geometry.GetGeometryType())
        featureDefnHeaders = feature_layer.GetLayerDefn()
        out_feature = ogr.Feature(featureDefnHeaders)
        out_feature.SetGeometry(plot_geometry)
        feature_layer.CreateFeature(out_feature)
        feature_ds.FlushCache()
        # Rasterize zone polygon to raster dsm
        gdal.RasterizeLayer(target_dsm, [1], feature_layer, burn_values=[1])
        feature_dsm_data = dsm_rb.ReadAsArray(rs_x_off, rs_y_off, rs_x_count, rs_y_count).astype(float)
        feature_dsm_band_mask = target_dsm.GetRasterBand(1)
        feature_dsm_data_mask = feature_dsm_band_mask.ReadAsArray(0, 0, rs_x_count, rs_y_count).astype(float)
        # Mask zone of raster dsm
        feature_dsm_raster_array = numpy.ma.masked_array(feature_dsm_data, numpy.logical_not(feature_dsm_data_mask))
        dsm_first_indexes, dsm_second_indexes = feature_dsm_raster_array.nonzero()
        # Rasterize zone polygon to raster dtm
        gdal.RasterizeLayer(target_dtm, [1], feature_layer, burn_values=[1])
        feature_dtm_data = dtm_rb.ReadAsArray(rs_x_off, rs_y_off, rs_x_count, rs_y_count).astype(float)
        feature_dtm_band_mask = target_dtm.GetRasterBand(1)
        feature_dtm_data_mask = feature_dtm_band_mask.ReadAsArray(0, 0, rs_x_count, rs_y_count).astype(float)
        # Mask zone of raster dsm
        feature_dtm_raster_array = numpy.ma.masked_array(feature_dtm_data, numpy.logical_not(feature_dtm_data_mask))
        dtm_first_indexes, dtm_second_indexes = feature_dtm_raster_array.nonzero()
        values = []
        mean_value = 0.
        min_value = 10000000.
        max_value = -10000000.
        number_of_pixels_in_frame = 0
        number_of_plant_pixels_in_frame = 0
        volume = 0
        for i in range(len(dsm_first_indexes)):
            fi = dsm_first_indexes[i]
            si = dsm_second_indexes[i]
            number_of_pixels_in_frame = number_of_pixels_in_frame + 1
            if not fi in dtm_first_indexes:
                continue
            if not si in dtm_second_indexes:
                continue
            # check mask
            dsm_height = feature_dsm_raster_array[dsm_first_indexes[i]][dsm_second_indexes[i]]
            dtm_height = feature_dtm_raster_array[dsm_first_indexes[i]][dsm_second_indexes[i]]
            value = dsm_height - dtm_height
            if value >= crop_minimum_height:
                values.append(value)
                if value < min_value:
                    min_value = value
                if value > max_value:
                    max_value = value
                mean_value = mean_value + value
                number_of_plant_pixels_in_frame = number_of_plant_pixels_in_frame + 1
                volume = volume + dems_pixel_area * value
        gcc = (float(number_of_plant_pixels_in_frame)) / (float(number_of_pixels_in_frame))
        if compute_GCC:
            feature.SetField(field_name_gcc, gcc)
        if compute_volume:
            feature.SetField(field_name_vol, volume * 1000.0)
        if not use_input_shp:
            out_layer.CreateFeature(feature)
        else:
            in_layer.SetFeature(feature)
    in_vec_ds = None
    out_vec_ds = None
    return str_error


def main():
    # ==================
    # parse command line
    # ==================
    usage = "usage: %prog [options] "
    parser = OptionParser(usage=usage)
    parser.add_option("--input_crops_frames_shp", dest="input_crops_frames_shp", action="store", type="string",
                      help="Crops frames shapefile", default=None)
    parser.add_option("--input_dsm", dest="input_dsm", action="store", type="string",
                      help="DSM geotiff", default=None)
    parser.add_option("--input_dtm", dest="input_dtm", action="store", type="string",
                      help="DTM geotiff", default=None)
    parser.add_option("--crop_minimum_height", dest="crop_minimum_height", action="store", type="string",
                      help="Crop minimum height, in meters", default=None)
    parser.add_option("--output_shp", dest="output_shp", action="store", type="string",
                      help="Output shapefile or none for use input shapefile", default=None)
    parser.add_option("--compute_GCC", dest="compute_GCC", action="store", type="int",
                      help="Compute GCC: 1-yes, 0-No", default=None)
    parser.add_option("--compute_volume", dest="compute_volume", action="store", type="int",
                      help="Compute volume: 1-yes, 0-No", default=None)
    parser.add_option("--date_from_dem_files", dest="date_from_dem_files", action="store", type="int",
                      help="Read date from DEM files: 1-yes, 0-No", default=None)
    parser.add_option("--dem_files_string_separator", dest="dem_files_string_separator", action="store", type="string",
                      help="DEM files string separator", default=None)
    parser.add_option("--dem_files_date_string_position", dest="dem_files_date_string_position", action="store", type="int",
                      help="DEM files date string position", default=None)
    parser.add_option("--date_format", dest="date_format", action="store", type="string",
                      help="Date format (%Y%m%d, ...)", default=None)
    parser.add_option("--date", dest="date", action="store", type="string",
                      help="Date value no from DEM files", default=None)
    (options, args) = parser.parse_args()
    if not options.input_crops_frames_shp \
        or not options.input_dsm \
        or not options.input_dtm \
        or not options.crop_minimum_height \
        or not options.output_shp \
        or not options.compute_GCC \
        or not options.compute_volume \
        or not options.date_from_dem_files \
        or not options.dem_files_string_separator \
        or not options.dem_files_date_string_position \
        or not options.date_format \
        or not options.date :
        parser.print_help()
        return
    input_crops_shapefile = options.input_crops_frames_shp
    if not exists(input_crops_shapefile):
        print("Error:\nInput crops shapefile does not exists:\n{}".format(input_crops_shapefile))
        return
    input_dsm = options.input_dsm
    if not exists(input_dsm):
        print("Error:\nInput DSM does not exists:\n{}".format(input_dsm))
        return
    input_dtm = options.input_dtm
    if not exists(input_dtm):
        print("Error:\nInput DTM does not exists:\n{}".format(input_dtm))
        return
    crop_minimum_height = options.crop_minimum_height
    if not is_number(crop_minimum_height):
        print("Error:\nInvalid value for crops minimum height: {}".
              format(crop_minimum_height))
        return
    crop_minimum_height = float(crop_minimum_height)
    compute_GCC = False
    if options.compute_GCC == 1:
        compute_GCC = True
    compute_volume = False
    if options.compute_volume == 1:
        compute_volume = True
    if not compute_GCC and not compute_volume:
        print("Error:\nNothing to do, no GCC or volume computation")
        return
    use_input_shp = True
    output_shp = options.output_shp
    if output_shp != 'none':
        use_input_shp = False
    date_from_dem_files = False
    if options.date_from_dem_files == 1:
        date_from_dem_files = True
    date_format = options.date_format.strip()
    if not date_from_dem_files:
        if options.date == 'none':
            print("Error:\nDate must be a value if not read from dems file name")
            return
        date = datetime.datetime.strptime(options.date, date_format)
    else:
        dem_files_string_separator = options.dem_files_string_separator
        dem_files_date_string_position = options.dem_files_date_string_position
        dsm_file_name_without_path = os.path.basename(input_dsm)
        dsm_file_name_values = dsm_file_name_without_path.split(dem_files_string_separator)
        if dem_files_date_string_position < 0 or dem_files_date_string_position > len(dsm_file_name_values):
            print("Error:\nInvalid value for dsm files date string position: {}".
                  format(str(dem_files_date_string_position)))
        str_date = dsm_file_name_values[dem_files_date_string_position-1]
        date_dsm = datetime.datetime.strptime(str_date, date_format)
        dtm_file_name_without_path = os.path.basename(input_dtm)
        dtm_file_name_values = dtm_file_name_without_path.split(dem_files_string_separator)
        if dem_files_date_string_position < 0 or dem_files_date_string_position > len(dtm_file_name_values):
            print("Error:\nInvalid value for dtm files date string position: {}".
                  format(str(dem_files_date_string_position)))
        str_date = dtm_file_name_values[dem_files_date_string_position-1]
        date_dtm = datetime.datetime.strptime(str_date, date_format)
        if date_dsm != date_dtm:
            print("Error:\nDifferents dates from DSM and DTM")
            return
        date = date_dsm
    str_date = str(date.strftime('%Y')[2:4]) + str(date.strftime('%m')) + str(date.strftime('%d'))
    str_error = process(input_crops_shapefile,
                        input_dsm,
                        input_dtm,
                        use_input_shp,
                        output_shp,
                        crop_minimum_height,
                        compute_GCC,
                        compute_volume,
                        str_date)
    if str_error:
        print("Error:\n{}".format(str_error))
        return
    print("... Process finished")


if __name__ == '__main__':
    main()
