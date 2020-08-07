# Extracting values in the buffer zones from the raster (tiff of each climtic variable)
ExtractValuesFromBuffers <- function(x){
  ext <- raster::extract(clim.var.raster,df[,c("longitude","latitude")],buffer=x) 
  ext <- lapply(ext, function(x) x[!is.na(x)])
  ext<- unlist(map(ext, var))
}