
# Calculate the slope, TPI, aspect, roughness and TRI index with the package 'raster' 
# Using data from SRTM.


library(raster)

# SRTM raster with the tiles merged in SAGA
srtm.rast <- raster("data/Topography/SrtmWGS84/MosaicTif/srtm_mosaic.tif")

x <- terrain(srtm.rast, opt=c("slope"), neighbors=8)
writeRaster(x, filename="data/Topography/TopoWithR/Slope_90m.grd", format="raster", overwrite=TRUE)

x <- terrain(srtm.rast, opt=c("tpi"), neighbors=8)
writeRaster(x, filename="data/Topography/TopoWithR/TPI_90m.grd",format="raster",overwrite=TRUE)

x <- terrain(srtm.rast, opt=c("aspect"), neighbors=8)
writeRaster(x, filename="data/Topography/TopoWithR/Aspect_90m.grd",format="raster", overwrite=TRUE)

x <- terrain(srtm.rast, opt=c("roughness"), neighbors=8)
writeRaster(x, filename="data/Topography/TopoWithR/Roughness_90m.grd",format="raster", overwrite=TRUE)

x <- terrain(srtm.rast, opt=c("tri"), neighbors=8)
writeRaster(x, filename="data/Topography/TopoWithR/TRI_90m.grd",format="raster",overwrite=TRUE)

