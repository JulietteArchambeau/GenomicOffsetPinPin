# Script to calculate the raster of the topographic variables
# -----------------------------------------------------------
#
# 07/08/2020
# Using saga version: 2.3.1


cd /home/juliette/Documents/GenomicOffset/GenomicOffsetPinPin/data/Topography/   # working directory



####### 1/ Download the SRTM tiles
####### --------------------------

# Get the GeoTIFF from: http://dwtkns.com/srtm/
# Another possibility (from the official website: http://srtm.csi.cgiar.org/)

# These GeoTIFF are in the WGS84 geographic coordinate system (long/lat).
# In order to use calculate the topographic variables in SAGA, we have to project these tiles in a projected coordinate system, e.g. UTM projections,
# which relies on the metre as unit of measure. 

# Resolution: 90m.




####### 2/ Convert TIF in SGRD
####### --------------------------

# saga uses grids in format .sgrd
# So we have to convert the SRTM tifs to grids (.sgrd).


# Tiles in the WGS 84 / UTM zone 29N (https://epsg.io/32629)
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_35_04.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_35_04.tif 
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_35_05.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_35_05.tif 
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_35_06.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_35_06.tif 

# Tiles in the WGS 84 / UTM zone 30N
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_36_03.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_36_03.tif 
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_36_04.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_36_04.tif 
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_36_05.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_36_05.tif 
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_36_06.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_36_06.tif 

# Tiles in the WGS 84 / UTM zone 31N
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_37_03.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_37_03.tif 
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_37_04.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_37_04.tif 
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_37_05.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_37_05.tif 

# Tiles in the WGS 84 / UTM zone 32N
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_38_04.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_38_04.tif 

# Tiles in the WGS 84 / UTM zone 33N
saga_cmd io_gdal 0 -GRIDS SrtmWGS84/SGRD/srtm_39_04.sgrd -FILES SrtmWGS84/OriginalTifsTiles/srtm_39_04.tif 





####### 3/ Mosaicking
####### -------------

# Creating a mosaic of all tiles
# At this step, the files are still in the WGS84 geographic coordinate system.

saga_cmd grid_tools 3 -GRIDS SrtmWGS84/SGRD/srtm_35_04.sgrd\;SrtmWGS84/SGRD/srtm_35_05.sgrd\;SrtmWGS84/SGRD/srtm_35_06.sgrd\;SrtmWGS84/SGRD/srtm_36_03.sgrd\;SrtmWGS84/SGRD/srtm_36_04.sgrd\;SrtmWGS84/SGRD/srtm_36_05.sgrd\;SrtmWGS84/SGRD/srtm_36_06.sgrd\;SrtmWGS84/SGRD/srtm_37_03.sgrd\;SrtmWGS84/SGRD/srtm_37_04.sgrd\;SrtmWGS84/SGRD/srtm_37_05.sgrd\;SrtmWGS84/SGRD/srtm_38_04.sgrd\;SrtmWGS84/SGRD/srtm_39_04.sgrd \
					  -TYPE 7 \
					  -OVERLAP 1 \
					  -BLEND_DIST 10.000000 \
					  -TARGET_OUT_GRID SrtmWGS84/SGRD/srtm_mosaic.sgrd



# save the file in tif for further analyses:
gdal_translate -of GTiff SrtmWGS84/SGRD/srtm_mosaic.sdat SrtmWGS84/MosaicTif/srtm_mosaic.tif




####### 4/ Reprojecting from WGS 84 to UTM 31N
####### --------------------------------------


# Following the SAGA manual: https://sagatutorials.files.wordpress.com/2016/02/saga_manual_english_cdu_june-2017.pdf
# Section 7.1 Coordinate transform.
# I choose to project in UTM 31N as it is the middle zone: https://epsg.io/32631

saga_cmd pj_proj4 0 -CRS_METHOD 0 -CRS_PROJ4 "+proj=longlat +datum=WGS84 +no_defs" -GRIDS SrtmWGS84/SGRD/srtm_mosaic.sgrd

saga_cmd pj_proj4 4 -SOURCE SrtmWGS84/SGRD/srtm_mosaic.sgrd \
					-CRS_METHOD 0 \
					-TARGET_USER_SIZE 90 \
					-CRS_PROJ4 "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs" \
					-TARGET_GRID SrtmUTM/srtm_mosaic_UTM31N.sgrd





####### 5/ Calculating the slope, aspect and curvature
####### ----------------------------------------------


# Units for slope and aspect are in radians (default option)

saga_cmd ta_morphometry 0 -ELEVATION SrtmUTM/srtm_mosaic_UTM31N.sgrd \
						  -SLOPE SlopeAspectCurvature/gridsUTM31N/slope_mosaic_UTM31N.sgrd \
						  -ASPECT SlopeAspectCurvature/gridsUTM31N/aspect_mosaic_UTM31N.sgrd \
						  -C_GENE SlopeAspectCurvature/gridsUTM31N/curvature_mosaic_UTM31N.sgrd

# save the file in tif for further analyses:			
gdal_translate -of GTiff SlopeAspectCurvature/gridsUTM31N/slope_mosaic_UTM31N.sdat SlopeAspectCurvature/TifsUTM31N/slope_mosaic_UTM31N.tif
gdal_translate -of GTiff SlopeAspectCurvature/gridsUTM31N/aspect_mosaic_UTM31N.sdat SlopeAspectCurvature/TifsUTM31N/aspect_mosaic_UTM31N.tif
gdal_translate -of GTiff SlopeAspectCurvature/gridsUTM31N/curvature_mosaic_UTM31N.sdat SlopeAspectCurvature/TifsUTM31N/curvature_mosaic_UTM31N.tif


			  
	
####### 6/ Calculating the topographic wetness index
####### --------------------------------------------


saga_cmd ta_preprocessor 5 -ELEV SrtmUTM/srtm_mosaic_UTM31N.sgrd -FILLED TWI/demfilledsinks_UTM31N.sgrd # bassins-versants
saga_cmd ta_hydrology 15 -DEM TWI/demfilledsinks_UTM31N.sgrd -TWI TWI/TWI_UTM31N.srgd # TWI
gdal_translate -of GTiff TWI/TWI_UTM31N.sdat TWI/TWI_UTM31N.tif # SAGA format -> GEOTif tif



# More info
# ---------

# https://www.researchgate.net/post/How_do_I_mosaic_images_of_two_different_UTM_zones_for_change_detection
# " One thing to note is that reprojecting the image data from one UTM zone into a different UTM zone for an image that already has projection information stored for it, such as a geotiff, will not impact the image pixel values. This is because it will simply convert the header coordinate information over to the new coordinate system origin and then each pixel x,y location is based on its row,column offset and pixel size distance from origin. This process is different than georectification where you would be applying a geometric transformation to each pixel location to warp it to match control points in some other coordinate system. The image would then need to be resampled and hence the pixel values change which could mess up your change detection scheme."


####### 7/ Reprojecting from UTM 31N to WGS 84
####### --------------------------------------


# I need to reproject in WGS 84 to include in a same stack the climatic, soil and topographic variables for GF and GDM analyses.


# Curvature
saga_cmd pj_proj4 0 -CRS_METHOD 0 -CRS_PROJ4 "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs" -GRIDS SlopeAspectCurvature/gridsUTM31N/curvature_mosaic_UTM31N.sgrd
saga_cmd pj_proj4 4 -SOURCE SlopeAspectCurvature/gridsUTM31N/curvature_mosaic_UTM31N.sgrd \
					-CRS_METHOD 0 \
					-CRS_PROJ4 "+proj=longlat +datum=WGS84 +no_defs" \
					-TARGET_GRID SlopeAspectCurvature/TifsWGS84/curvature_mosaic_WGS84.sgrd			
gdal_translate -of GTiff SlopeAspectCurvature/TifsWGS84/curvature_mosaic_WGS84.sdat SlopeAspectCurvature/TifsWGS84/curvature_mosaic_WGS84.tif


# Aspect
saga_cmd pj_proj4 0 -CRS_METHOD 0 -CRS_PROJ4 "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs" -GRIDS SlopeAspectCurvature/gridsUTM31N/aspect_mosaic_UTM31N.sgrd
saga_cmd pj_proj4 4 -SOURCE SlopeAspectCurvature/gridsUTM31N/aspect_mosaic_UTM31N.sgrd \
					-CRS_METHOD 0 \
					-CRS_PROJ4 "+proj=longlat +datum=WGS84 +no_defs" \
					-TARGET_GRID SlopeAspectCurvature/TifsWGS84/aspect_mosaic_WGS84.sgrd		
gdal_translate -of GTiff SlopeAspectCurvature/TifsWGS84/aspect_mosaic_WGS84.sdat SlopeAspectCurvature/TifsWGS84/aspect_mosaic_WGS84.tif


# Slope
saga_cmd pj_proj4 0 -CRS_METHOD 0 -CRS_PROJ4 "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs" -GRIDS SlopeAspectCurvature/gridsUTM31N/slope_mosaic_UTM31N.sgrd
saga_cmd pj_proj4 4 -SOURCE SlopeAspectCurvature/gridsUTM31N/slope_mosaic_UTM31N.sgrd \
					-CRS_METHOD 0 \
					-CRS_PROJ4 "+proj=longlat +datum=WGS84 +no_defs" \
					-TARGET_GRID SlopeAspectCurvature/TifsWGS84/slope_mosaic_WGS84.sgrd
gdal_translate -of GTiff SlopeAspectCurvature/TifsWGS84/slope_mosaic_WGS84.sdat SlopeAspectCurvature/TifsWGS84/slope_mosaic_WGS84.tif



####### 8/ Calculating the Topographic Ruggedness Index (TRI)
####### -----------------------------------------------------

# Calculated based on: Riley, S.J., De Gloria, S.D., Elliot, R. (1999): A Terrain Ruggedness that Quantifies Topographic Heterogeneity. Intermountain Journal of Science, Vol.5, No.1-4, pp.23-27.

saga_cmd ta_morphometry 16 -DEM SrtmUTM/srtm_mosaic_UTM31N.sgrd -TRI TRI/gridsUTM31N/TRI.sgrd

# Reproject in WGS84
saga_cmd pj_proj4 0 -CRS_METHOD 0 -CRS_PROJ4 "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs" -GRIDS TRI/gridsUTM31N/TRI.sgrd
saga_cmd pj_proj4 4 -SOURCE TRI/gridsUTM31N/TRI.sgrd \
					-CRS_METHOD 0 \
					-CRS_PROJ4 "+proj=longlat +datum=WGS84 +no_defs" \
					-TARGET_GRID TRI/gridsWGS84/TRI_WGS84.sgrd		
gdal_translate -of GTiff TRI/gridsWGS84/TRI_WGS84.sdat TRI/TifsWGS84/TRI_WGS84.tif
