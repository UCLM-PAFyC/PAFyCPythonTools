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

def process_cwsith(input_shp,
                   kmeans_clusters,
                   percentile_minimum_threshold,
                   input_orthomosaic,
                   temperature,
                   relative_humidity,
                   upper_line_coef_a,
                   upper_line_coef_b,
                   lower_line_coef_a,
                   lower_line_coef_b,
                   factor_to_temperature,
                   str_date):
    str_error = None

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
    temperature = float(str_temperature)
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
    str_error = process_cwsith(input_shp,
                               kmeans_clusters,
                               percentile_minimum_threshold,
                               input_orthomosaic,
                               temperature,
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
