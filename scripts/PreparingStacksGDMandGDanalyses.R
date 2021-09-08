#
# Preparing data for GDM and GF analyses
# --------------------------------------"


library(raster)
library(stringr)

# Maritime pine distribution from EUFORGEN (shapefile)
PinpinDistri  <- shapefile('../../Pinpin_Clonapin/maps/pinus_pinaster_distribution/Pinus_pinaster_EUFORGEN.shp')

b <- as(extent(5, 10, 33,37.5), 'SpatialPolygons')
crs(b) <- crs(PinpinDistri)
PinpinDistri <- gDifference(PinpinDistri,b)
plot(PinpinDistri)

#----------------------------------------------------------------------------------------------------------- #

#  I) CURRENT CLIMATE                                                                                        ####
#     ---------------                                                                                        #


#  A) Set AvgSand                                                                                            ####
#     -----------                                                                                            #


#  Climatic variables ####
# > from set AvgSand and AvgWater
grids <- list.files("data/climate/CurrentClimate/EuMedClim_Tiff_1901_2009/" , pattern = "*1901-2009_with_34_sampling_points.tif$")
grids <- grids[grepl("bio1_|bio12_|bio14_|bio2_",grids)==T]
grids <- c(grids[[1]],grids[[4]],grids[[2]],grids[[3]])
rast.clim <- raster::stack(paste0("data/climate/CurrentClimate/EuMedClim_Tiff_1901_2009/", grids))
names(rast.clim) <- c("bio1","bio2","bio12","bio14")
rast.clim <- mask(rast.clim,PinpinDistri, updatevalue=NA)


# TRI raster ####
# > The models are fitted with TRI values from raster at 90m resolution, and the 
# projections will be on raster at 1km resolution (to have the same resolution as WorldClim variables).
rast.tri <- raster("data/Topography/TopoWithR/TRI_WGS84_1km.grd")
rast.tri <- crop(rast.tri,extent(rast.clim)) # attribute the extent of climatic raster to TRI rasters
rast.tri <- raster::resample(rast.tri,rast.clim) # same number of cells between the TRI and climatic rasters
names(rast.tri) <- c("TRI")
rast.tri <- mask(rast.tri,PinpinDistri,updatevalue=NA) # keep only cells within maritime pine distribution 

# Soil rasters ####
rast.soil.sand <- stack("data/soil/STU_EU_T_SAND_WGS84.tif","data/soil/STU_EU_DEPTH_ROOTS_WGS84.tif")
names(rast.soil.sand) <- c("sand_top","depth_roots")
rast.soil.sand <- raster::crop(rast.soil.sand,extent(rast.clim)) # attribute the extent of climatic raster to soil raster
rast.soil.sand <- raster::resample(rast.soil.sand,rast.clim) # same number of cells between the soil and climatic rasters
rast.soil.sand <- mask(rast.soil.sand,PinpinDistri,updatevalue=NA) # keep only cells within maritime pine distribution 

# >> Merge and save rasters
stackall <- stack(rast.clim,rast.soil.sand,rast.tri)
writeRaster(stackall, filename="data/StacksEnvVars/StackAvgSand.grd", format="raster",overwrite=TRUE)


#----------------------------------------------------------------------------------------------------------- #
#  B) Set AvgWater                                                                                           ####
#     -----------                                                                                            #

# Soil rasters ####
rast.soil.water <- stack("data/soil/STU_EU_T_TAWC_WGS84.tif","data/soil/STU_EU_DEPTH_ROOTS_WGS84.tif")
names(rast.soil.water) <- c("water_top","depth_roots")
rast.soil.water <- raster::crop(rast.soil.water,extent(rast.clim)) # attribute the extent of climatic raster to soil raster
rast.soil.water <- raster::resample(rast.soil.water,rast.clim) # same number of cells between the soil and climatic rasters
rast.soil.water <- mask(rast.soil.water,PinpinDistri,updatevalue=NA) # keep only cells within maritime pine distribution 

# >> Merge and save rasters
stackall <- stack(rast.clim,rast.soil.water,rast.tri)
writeRaster(stackall, filename="data/StacksEnvVars/StackAvgWater.grd", format="raster",overwrite=TRUE)


#----------------------------------------------------------------------------------------------------------- #
#  C) Set ExtSand                                                                                            ####
#     -----------                                                                                            #
  
#  Climatic variables ####
# > from set ExtSand and ExtWater
grids <- list.files("data/climate/CurrentClimate/EuMedClim_Tiff_1901_2009/" , pattern = "*1901-2009_with_34_sampling_points.tif$")
grids <- grids[grepl("bio13_|bio14_|bio5_|bio6_",grids)==T]
rast.clim <- raster::stack(paste0("data/climate/CurrentClimate/EuMedClim_Tiff_1901_2009/", grids))
names(rast.clim) <- c("bio13","bio14","bio5","bio6")
rast.clim <- mask(rast.clim,PinpinDistri, updatevalue=NA)

# >> Merge and save rasters
stackall <- stack(rast.clim,rast.soil.sand,rast.tri)
writeRaster(stackall, filename="data/StacksEnvVars/StackExtSand.grd", format="raster",overwrite=TRUE)



#----------------------------------------------------------------------------------------------------------- #
#  D) Set ExtWater                                                                                           ####
#     -----------                                                                                            #


# >> Merge and save rasters
stackall <- stack(rast.clim,rast.soil.water,rast.tri)
writeRaster(stackall, filename="data/StacksEnvVars/StackExtWater.grd", format="raster",overwrite=TRUE)



##############################################################################################################################"

#  II) FUTURE CLIMATE                                                                                       ####
#     --------------                                                                                        #


#  A) Set AvgWater                                                                                           ####
#     -----------                                                                                            #

rast.to.crop <- raster("data/climate/CurrentClimate/EuMedClim_Tiff_1901_2009/map_bio1_1901-2009_with_34_sampling_points.tif")

GCMs <- list.files("data/climate/FutureClimate/share/spatial03/worldclim/cmip6/7_fut/2.5m/")


# 2041-2060 - SSP370 - All GCMs                                                                               ####

for(i in GCMs){
  
  rast.clim <- stack(paste0("data/climate/FutureClimate/share/spatial03/worldclim/cmip6/7_fut/2.5m/",i,"/ssp370/wc2.1_2.5m_bioc_",i,"_ssp370_2041-2060.tif"))
  
  # I understand that these 19 variables here are the WorldClim bioclimatic variables averaged over the period 2041/2060.
  
  select.layers <- names(rast.clim)[grepl("2060.1$|2060.2|2060.12|2060.14",names(rast.clim))==T]
  rast.clim <- raster::subset(rast.clim,select.layers)
  names(rast.clim) <- c("bio1","bio2","bio12","bio14")
  rast.clim <- rast.clim[[c("bio1","bio2","bio12","bio14")]]
  rast.clim <- raster::crop(rast.clim,extent(rast.to.crop))
  rast.clim <- mask(rast.clim,PinpinDistri, updatevalue=NA)
  
  
  # Soil rasters
  rast.soil <- stack("data/soil/STU_EU_T_TAWC_WGS84.tif","data/soil/STU_EU_DEPTH_ROOTS_WGS84.tif")
  names(rast.soil) <- c("water_top","depth_roots")
  rast.soil <- raster::crop(rast.soil,extent(rast.clim)) # attribute the extent of climatic raster to soil raster
  rast.soil <- raster::resample(rast.soil,rast.clim) # same number of cells between the soil and climatic rasters
  rast.soil <- mask(rast.soil,PinpinDistri,updatevalue=NA) # keep only cells within maritime pine distribution 
  
  # TRI raster
  rast.tri <- raster("data/Topography/TopoWithR/TRI_WGS84_1km.grd")
  rast.tri <- crop(rast.tri,extent(rast.clim)) # attribute the extent of climatic raster to TRI rasters
  rast.tri <- raster::resample(rast.tri,rast.clim) # same number of cells between the TRI and climatic rasters
  names(rast.tri) <- c("TRI")
  rast.tri <- mask(rast.tri,PinpinDistri,updatevalue=NA) # keep only cells within maritime pine distribution 
  
  # >> Save without altitude
  stackall <- stack(rast.clim,rast.soil,rast.tri)
  
  writeRaster(stackall, 
              filename=paste0("data/StacksEnvVars/FutClimStacks/2041to2060/SSP370/StackAvgWater_",str_remove_all(i,"-"),".grd"), format="raster",overwrite=TRUE)
  
}


# 2041-2060 - SSP585 - All GCMs ####

GCMs <- GCMs[!GCMs=="GFDL-ESM4"]

for(i in GCMs){
  
  rast.clim <- stack(paste0("data/FutureClimate/share/spatial03/worldclim/cmip6/7_fut/2.5m/",i,"/ssp585/wc2.1_2.5m_bioc_",i,"_ssp585_2041-2060.tif"))
  
  select.layers <- names(rast.clim)[grepl("2060.1$|2060.2|2060.12|2060.14",names(rast.clim))==T]
  rast.clim <- raster::subset(rast.clim,select.layers)
  names(rast.clim) <- c("bio1","bio2","bio12","bio14")
  rast.clim <- rast.clim[[c("bio1","bio2","bio12","bio14")]]
  rast.clim <- raster::crop(rast.clim,extent(rast.to.crop))
  rast.clim <- mask(rast.clim,PinpinDistri, updatevalue=NA)
  
  
  # Soil rasters
  rast.soil <- stack("data/soil/STU_EU_T_TAWC_WGS84.tif","data/soil/STU_EU_DEPTH_ROOTS_WGS84.tif")
  names(rast.soil) <- c("water_top","depth_roots")
  rast.soil <- raster::crop(rast.soil,extent(rast.clim)) # attribute the extent of climatic raster to soil raster
  rast.soil <- raster::resample(rast.soil,rast.clim) # same number of cells between the soil and climatic rasters
  rast.soil <- mask(rast.soil,PinpinDistri,updatevalue=NA) # keep only cells within maritime pine distribution 
  
  # TRI raster
  rast.tri <- raster("data/Topography/TopoWithR/TRI_WGS84_1km.grd")
  rast.tri <- crop(rast.tri,extent(rast.clim)) # attribute the extent of climatic raster to TRI rasters
  rast.tri <- raster::resample(rast.tri,rast.clim) # same number of cells between the TRI and climatic rasters
  names(rast.tri) <- c("TRI")
  rast.tri <- mask(rast.tri,PinpinDistri,updatevalue=NA) # keep only cells within maritime pine distribution 
  
  # >> Save without altitude
  stackall <- stack(rast.clim,rast.soil,rast.tri)
  
  writeRaster(stackall, 
              filename=paste0("data/StacksEnvVars/FutClimStacks/2041to2060/SSP585/StackAvgWater_",str_remove_all(i,"-"),".grd"), format="raster",overwrite=TRUE)
  
}



##############################################################################################################################"

#  III) CURRENT DATA FOR GENOMIC OFFSET CALCULATION                                                          ####
#       -------------------------------------------                                                          #


#  A) Set AvgWater                                                                                           ####
#     -----------                                                                                            #



rast.to.crop <- stack("data/StacksEnvVars/FutClimStacks/2041to2060/SSP585/StackAvgWater_BCCCSM2MR.grd")

# Climatic rasters ####
grids <- list.files("data/climate/CurrentClimate/EuMedClim_Tiff_1901_2009/" , pattern = "*1901-2009_with_34_sampling_points.tif$")
grids <- grids[grepl("bio1_|bio2_|bio12_|bio14_",grids)==T]
grids <- c(grids[[1]],grids[[4]],grids[[2]],grids[[3]])
rast.clim <- raster::stack(paste0("data/climate/CurrentClimate/EuMedClim_Tiff_1901_2009/", grids))
names(rast.clim) <- c("bio1","bio2","bio12","bio14")
rast.clim <- raster::crop(rast.clim,extent(rast.to.crop))
rast.clim <- raster::resample(rast.clim,rast.to.crop)
rast.clim <- mask(rast.clim,PinpinDistri, updatevalue=NA)

# Soil rasters ####
rast.soil <- stack("data/soil/STU_EU_T_TAWC_WGS84.tif","data/soil/STU_EU_DEPTH_ROOTS_WGS84.tif")
names(rast.soil) <- c("water_top","depth_roots")
rast.soil <- raster::crop(rast.soil,extent(rast.clim)) # attribute the extent of climatic raster to soil raster
rast.soil <- raster::resample(rast.soil,rast.clim) # same number of cells between the soil and climatic rasters
rast.soil <- mask(rast.soil,PinpinDistri,updatevalue=NA) # keep only cells within maritime pine distribution 

# TRI raster ####
rast.tri <- raster("data/Topography/TopoWithR/TRI_WGS84_1km.grd")
rast.tri <- crop(rast.tri,extent(rast.clim)) # attribute the extent of climatic raster to TRI rasters
rast.tri <- raster::resample(rast.tri,rast.clim) # same number of cells between the TRI and climatic rasters
names(rast.tri) <- c("TRI")
rast.tri <- mask(rast.tri,PinpinDistri,updatevalue=NA) # keep only cells within maritime pine distribution 


stackall <- stack(rast.clim,rast.soil,rast.tri)
writeRaster(stackall, filename="data/StacksEnvVars/CurStacksResWorldClim/StackAvgWaterResWorldClim.grd", 
            format="raster",
            overwrite=TRUE)

