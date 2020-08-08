# Extracting values in the buffer zones from the raster (tiff of each climtic variable)
ExtractValuesFromBuffers <- function(x){
  
  if(crs(clim.var.raster)@projargs == "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"){
    ext <- raster::extract(clim.var.raster,df[,c("longitude","latitude")],buffer=x) 
  } else {
    d <- df[,c("longitude","latitude")]
    coordinates(d) <- c("longitude", "latitude")
    proj4string(d) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
    d <- spTransform(d, crs(clim.var.raster))
    ext <- raster::extract(clim.var.raster,d,buffer=x) 
  }
  
  ext <- lapply(ext, function(x) x[!is.na(x)])
  ext<- unlist(map(ext, var))
}






