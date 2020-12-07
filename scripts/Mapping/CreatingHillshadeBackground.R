####################################
# Hillshade background for the maps
####################################


# https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1-0-and-derived-products/slope?tab=download
# Need to create an account

# we took the hillshade map from here:
# https://www.eea.europa.eu/data-and-maps/data/digital-elevation-model-of-europe

background <- raster("data/maps/Hillshade/EuropeanData/hillshade1x1.tif")
background <- projectRaster(background, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
e <- extent(-10, 14, 31, 50)
background <- crop(background, e)

writeRaster(background, filename="data/maps/Hillshade/background.grd")
