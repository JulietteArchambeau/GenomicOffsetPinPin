CreateRasterBurnedAreas <- function(start_date = NULL,
                                    end_date = NULL){
  
  if (is.null(start_date)) stop("Invalid start_date")
  if (is.null(end_date)) stop("Invalid valid end_date")
  
  base_url <- "ftp://fuoco.geog.umd.edu/gfed4"
  
  
  seq_of_dates <- seq.Date(from = as.Date(start_date),
                           to = as.Date(end_date),
                           by = "month")
  my_dates <- substr(x = gsub("-", "", as.character(seq_of_dates)),
                     start = 1, stop = 6)
  # Assemble file names
  fnms <- paste0(base_url, "/monthly/GFED4.0_MQ_", my_dates, "_BA.hdf")
  
  
  
  # Loop through dates to populate the stack
  for (i in seq_along(fnms)) {
    
    message(paste("Downloading", fnms[i]))
    # Download the file
    x <- try(RCurl::getBinaryURL(fnms[i],
                                 userpwd = "fire:burnt",
                                 ftp.use.epsv = FALSE,
                                 .opts = list(timeout = 1000,
                                              connecttimeout = 1000,
                                              maxredirs = 20)),
             silent = FALSE)
    
    Sys.sleep(3) # add by Charlie Pauvert
    
    if (class(x) == "try-error") {
      
      stop("Server currently unavailable, please try again later.")
      
    } else {
      
      
      # Download the file
      hdf_file_path <- tempfile(fileext = ".hdf")
      writeBin(x, con = hdf_file_path)
      
      # Get subdataset names
      sds <- gdalUtils::get_subdatasets(hdf_file_path)[1]
      message(paste("Importing subdataset:", sds, collapse = "\n"))
      
      # Translate subdataset from hdf file to tiff, this is needed because no
      # cross-platform direct translation to nc is available
      temp_tif_file <- tempfile(fileext = ".tif")
      gdalUtils::gdal_translate(src_dataset = sds[1],
                                dst_dataset = temp_tif_file)
      
      if (i == 1){
        s <- raster::raster(temp_tif_file)
      }else{
        s <- raster::stack(s, raster::raster(temp_tif_file))
      }
      
    }
    
  }
  
  # writeRaster(s,paste0("GFED4_",
  #                      substr(x = gsub("-", "",start_date),start = 1, stop = 6),"_",
  #                      substr(x = gsub("-", "",end_date),start = 1, stop = 6),
  #                      ".grd"), format="raster",overwrite=TRUE)  
  return(s)
  
  
}
