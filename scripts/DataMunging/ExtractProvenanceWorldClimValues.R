

# Add to the DF "AllDataPhenoClimSoil" the provenance climatic values from WorldClim 
# ================================================================================== #

library(raster)    # CRAN v3.3-13 
library(tidyverse) # CRAN v1.3.0 


# Load the data used in paper one with EumedClim climatic data and soil data.
data <- readRDS(file="data/AllDataPhenoClimSoil.RDS")

# create a df of the provenance coordinates
xy <- unique(data[,c("prov","longitude_prov","latitude_prov")])
xy <- xy[!(xy$prov=="ROD"),]
colnames(xy) <- c("prov","longitude","latitude")
xy

# Path to WorldClim data
path="data/climate/CurrentClimate/WorldClim30sec_1970_2000/"

# List the WorldClim files (the bioclimatic variables)
myFiles <- list.files(path,pattern=".tif")

vars <- str_sub(myFiles,11,-5) %>% str_remove("_") %>% str_c("WC")

for (i in 1:length(myFiles)){
  rast <- raster(paste0(path,myFiles[i]))
  xy[,vars[i]] <- raster::extract(rast,xy[,c("longitude","latitude")])}

