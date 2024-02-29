echo off
REM set PROCESS_PATH=D:/PAFyCToolsGui/20220426_Tarazona_Vid_A6000/output
set OSGEO4W_ROOT=C:/Program Files/QGIS 3.22.12
set TOOL=D:/PAFyCToolsGui/20220426_Tarazona_Vid_A6000/output/CropsCharacterizationFromPhotogrammetricGeomaticProducts.py
REM cd /d "%PROCESS_PATH%"
call "%OSGEO4W_ROOT%\bin\o4w_env.bat"
REM ndvi, kmeans
REM python E:\python\PAFyCPythonTools\CropMonitoringFromPhotogrammetricGeomaticProducts.py --input_crops_frames_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_initial.shp" --method_data ndvi --method_segmentation  kmeans --kmeans_clusters 6 --input_orthomosaic "D:\PAFyCToolsGui\20230630_Tarazona_MULTI_Vid\ORT-303-RFL-Tarazona_Vid_20230630_25830_4cm.tif" --blue_band_position 1 --green_band_position 2 --red_band_position 3 --nir_band_position 5 --soil_maximum_rgb_reflectance 0.08  --crop_minimum_value 0.5 --factor_to_reflectance 1.0 --output_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_initial_ndvi_kmeans.shp" --date_from_orthomosaic_file=1 --orthomosaic_file_string_separator="_" --orthomosaic_file_date_string_position=3 --date_format="%%Y%%m%%d" --date=none

REM ndvi, percentile
python E:\python\PAFyCPythonTools\CropMonitoringFromPhotogrammetricGeomaticProducts.py --input_crops_frames_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_initial.shp" --method_data ndvi --method_segmentation  percentile --percentile_minimum_threshold 0.05 --input_orthomosaic "D:\PAFyCToolsGui\20230630_Tarazona_MULTI_Vid\ORT-303-RFL-Tarazona_Vid_20230630_25830_4cm.tif" --blue_band_position 1 --green_band_position 2 --red_band_position 3 --nir_band_position 5 --soil_maximum_rgb_reflectance 0.08  --crop_minimum_value 0.5 --factor_to_reflectance 1.0 --output_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_initial_ndvi_percentile.shp" --date_from_orthomosaic_file=1 --orthomosaic_file_string_separator="_" --orthomosaic_file_date_string_position=3 --date_format="%%Y%%m%%d" --date=none

REM gcc_or_vol, gcc, kmeans
REM python E:\python\PAFyCPythonTools\CropMonitoringFromPhotogrammetricGeomaticProducts.py --input_crops_frames_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_gcc_vol_allDates.shp" --method_data gcc --crop_minimum_value 0.02 --method_segmentation  kmeans --kmeans_clusters 6 --input_field "220905_gcc" --output_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_gcc_vol_allDates_gcc_kmeans.shp" --date_format="%Y%m%d"

REM gcc_or_vol, gcc, percentile
REM python E:\python\PAFyCPythonTools\CropMonitoringFromPhotogrammetricGeomaticProducts.py --input_crops_frames_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_gcc_vol_allDates.shp" --method_data gcc --crop_minimum_value 0.02 --method_segmentation  percentile --percentile_minimum_threshold 0.05 --input_field "220905_gcc" --output_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_gcc_vol_allDates_gcc_percentile.shp" --date_format="%Y%m%d"

REM gcc_or_vol, vol, kmeans
REM python E:\python\PAFyCPythonTools\CropMonitoringFromPhotogrammetricGeomaticProducts.py --input_crops_frames_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_gcc_vol_allDates.shp" --method_data vol --crop_minimum_value 0.1 --method_segmentation  kmeans --kmeans_clusters 6 --input_field "220812_vol" --output_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_gcc_vol_allDates_vol_kmeans.shp" --date_format="%Y%m%d"

REM gcc_or_vol, vol, percentile
REM python E:\python\PAFyCPythonTools\CropMonitoringFromPhotogrammetricGeomaticProducts.py --input_crops_frames_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_gcc_vol_allDates.shp" --method_data vol --crop_minimum_value 0.1 --method_segmentation  percentile --percentile_minimum_threshold 0.05 --input_field "220812_vol" --output_shp "D:\PAFyCToolsGui\20220426_Tarazona_Vid_A6000\output\vines_frames_gcc_vol_allDates_vol_percentile.shp" --date_format="%Y%m%d"
