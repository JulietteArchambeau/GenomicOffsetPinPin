#
# Preparing data for GDM and GF analyses
# --------------------------------------"


library(raster)

# Maritime pine distribution from EUFORGEN (shapefile)
shp  <- shapefile('../../Pinpin_Clonapin/maps/pinus_pinaster_distribution/Pinus_pinaster_EUFORGEN.shp')


# --------------------------------------------------------------------------------------------------------"

# >>> CURRENT DATA FOR MAPS ####

#  Climatic rasters - Set 1 ####
grids <- list.files("data/climate/Tiff_1901_2009/" , pattern = "*1901-2009_with_34_sampling_points.tif$")
grids <- grids[grepl("bio5_|bio14_|bio12_|bio1_",grids)==T]
# grids <- c(grids[[1]],grids[[4]],grids[[2]],grids[[3]])
rast <- raster::stack(paste0("data/climate/Tiff_1901_2009/", grids))
names(rast) <- c("bio1","bio12","bio14","bio5")
rast <- mask(rast,shp, updatevalue=NA)


# SRTM altitude raster ####
rast.alt <- raster("data/Topography/SrtmWGS84/MosaicTif/srtm_mosaic.tif")
rast.alt <- crop(rast.alt,extent(rast)) # attribute the extent of climatic raster to altitude rasters
rast.alt <- raster::resample(rast.alt,rast) # same number of cells between the altitude and climatic rasters
names(rast.alt) <- c("altitude")
rast.alt <- mask(rast.alt,shp,updatevalue=NA) # keep only cells within maritime pine distribution 

# Soil rasters ####
stack.soil <- stack("data/soil/STU_EU_T_SAND_WGS84.tif","data/soil/STU_EU_DEPTH_ROOTS_WGS84.tif")
names(stack.soil) <- c("sand_top","depth_roots")
stack.soil <- raster::crop(stack.soil,extent(rast)) # attribute the extent of climatic raster to soil raster
stack.soil <- raster::resample(stack.soil,rast) # same number of cells between the soil and climatic rasters
stack.soil <- mask(stack.soil,shp,updatevalue=NA) # keep only cells within maritime pine distribution 

# >> Save with altitude
stackall <- stack(rast.alt,rast,stack.soil)
writeRaster(stackall, filename="data/StacksEnvVars/StackAltSoilClimSet1.grd", format="raster",overwrite=TRUE)

# >> Save without altitude
stackall <- stack(rast,stack.soil)
writeRaster(stackall, filename="data/StacksEnvVars/StackSoilClimSet1.grd", format="raster",overwrite=TRUE)


#--------------------------------------------------------------------------"

# Climatic rasters - Set 2 ####
grids <- list.files("data/climate/Tiff_1901_2009/" , pattern = "*1901-2009_with_34_sampling_points.tif$")
grids <- grids[grepl("bio12_|bio1_|ppet.min_|bio6_",grids)==T]
rast <- raster::stack(paste0("data/climate/Tiff_1901_2009/", grids))
names(rast) <- c("bio1","bio12","bio6","ppet.min")
rast <- mask(rast,shp, updatevalue=NA)

# >> Save with altitude
stackall <- stack(rast.alt,rast,stack.soil)
writeRaster(stackall, filename="data/StacksEnvVars/StackAltSoilClimSet2.grd", format="raster",overwrite=TRUE)

# >> Save without altitude
stackall <- stack(rast,stack.soil)
writeRaster(stackall, filename="data/StacksEnvVars/StackSoilClimSet2.grd", format="raster",overwrite=TRUE)


#--------------------------------------------------------------------------"

# Climatic rasters - Set 3 ####
grids <- list.files("data/climate/Tiff_1901_2009/" , pattern = "*1901-2009_with_34_sampling_points.tif$")
grids <- grids[grepl("bio14_|bio1_|ppet.mean_",grids)==T]
rast <- raster::stack(paste0("data/climate/Tiff_1901_2009/", grids))
names(rast) <- c("bio1","bio14","ppet.mean")
rast <- mask(rast,shp, updatevalue=NA)

# >> Save with altitude
stackall <- stack(rast.alt,rast,stack.soil)
writeRaster(stackall, filename="data/StacksEnvVars/StackAltSoilClimSet3.grd", format="raster",overwrite=TRUE)


#--------------------------------------------------------------------------"

# Climatic rasters - Set 4 ####
grids <- list.files("data/climate/Tiff_1901_2009/" , pattern = "*1901-2009_with_34_sampling_points.tif$")
grids <- grids[grepl("bio14_|bio1_|ppet.mean_|ppet.min",grids)==T]
rast <- raster::stack(paste0("data/climate/Tiff_1901_2009/", grids))
names(rast) <- c("bio1","bio14","ppet.mean","ppet.min")
rast <- mask(rast,shp, updatevalue=NA)

# >> Save without altitude
stackall <- stack(rast,stack.soil)
writeRaster(stackall,"data/StacksEnvVars/StackSoilClimSet4.grd", format="raster",overwrite=TRUE)


#--------------------------------------------------------------------------"

# Climatic rasters - Set 5 ####
grids <- list.files("data/climate/Tiff_1901_2009/" , pattern = "*1901-2009_with_34_sampling_points.tif$")
grids <- grids[grepl("bio13_|bio14_|bio5_|bio6_",grids)==T]
rast <- raster::stack(paste0("data/climate/Tiff_1901_2009/", grids))
names(rast) <- c("bio13","bio14","bio5","bio6")
rast <- mask(rast,shp, updatevalue=NA)


# >> Save without altitude
stackall <- stack(rast,stack.soil)
writeRaster(stackall,"data/StacksEnvVars/StackSoilClimSet5.grd", format="raster",overwrite=TRUE)



# --------------------------------------------------------------------------------------------------------"

# >>> FUTURE DATA ####


rast.to.crop <- raster("data/climate/Tiff_1901_2009/map_bio1_1901-2009_with_34_sampling_points.tif")

ssp245_2041_2060 <- stack("data/FutureClimate/share/spatial03/worldclim/cmip6/7_fut/2.5m/BCC-CSM2-MR/ssp245/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2041-2060.tif")

# I understand that these variables here are the WorldClim bioclimatic variables averaged over the period 2041/2060.

select.layers <- names(ssp245_2041_2060)[grepl("2060.13|2060.14|2060.5|2060.6",names(ssp245_2041_2060))==T]
rast <- raster::subset(ssp245_2041_2060,select.layers)
names(rast) <- c("bio5","bio6","bio13","bio14")
rast <- rast[[c("bio13","bio14","bio5","bio6")]]
rast <- raster::crop(rast,extent(rast.to.crop))
rast <- mask(rast,shp, updatevalue=NA)


# Soil rasters ####
stack.soil <- stack("data/soil/STU_EU_T_SAND_WGS84.tif","data/soil/STU_EU_DEPTH_ROOTS_WGS84.tif")
names(stack.soil) <- c("sand_top","depth_roots")
stack.soil <- raster::crop(stack.soil,extent(rast)) # attribute the extent of climatic raster to soil raster
stack.soil <- raster::resample(stack.soil,rast) # same number of cells between the soil and climatic rasters
stack.soil <- mask(stack.soil,shp,updatevalue=NA) # keep only cells within maritime pine distribution 

# >> Save without altitude
stackall <- stack(rast,stack.soil)
writeRaster(stackall, filename="data/StacksEnvVars/FutClimStackSoilClimSet5.grd", format="raster",overwrite=TRUE)


# --------------------------------------------------------------------------------------------------------"

# >>> CURRENT DATA FOR GENOMIC OFFSET CALCULATION ####

rast.to.crop <- stack("data/StacksEnvVars/FutClimStackSoilClimSet5.grd")

# Climatic rasters - Set 5 ####
grids <- list.files("data/climate/Tiff_1901_2009/" , pattern = "*1901-2009_with_34_sampling_points.tif$")
grids <- grids[grepl("bio13_|bio14_|bio5_|bio6_",grids)==T]
rast <- raster::stack(paste0("data/climate/Tiff_1901_2009/", grids))
names(rast) <- c("bio13","bio14","bio5","bio6")
rast <- raster::crop(rast,extent(rast.to.crop))
rast <- raster::resample(rast,rast.to.crop)
rast <- mask(rast,shp, updatevalue=NA)

# Soil rasters ####
stack.soil <- stack("data/soil/STU_EU_T_SAND_WGS84.tif","data/soil/STU_EU_DEPTH_ROOTS_WGS84.tif")
names(stack.soil) <- c("sand_top","depth_roots")
stack.soil <- raster::crop(stack.soil,extent(rast)) # attribute the extent of climatic raster to soil raster
stack.soil <- raster::resample(stack.soil,rast) # same number of cells between the soil and climatic rasters
stack.soil <- mask(stack.soil,shp,updatevalue=NA) # keep only cells within maritime pine distribution 

# >> Save without altitude
stackall <- stack(rast,stack.soil)
writeRaster(stackall, filename="data/StacksEnvVars/StackSoilClimSet5LowRes.grd", format="raster",overwrite=TRUE)


