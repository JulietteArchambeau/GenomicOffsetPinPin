################################################
# Mask of the maritime distribution for the maps
################################################

library(raster)

# Maritime pine distribution from EUFORGEN (shapefile)
PinpinDistri  <- shapefile('../../Pinpin_Clonapin/maps/pinus_pinaster_distribution/Pinus_pinaster_EUFORGEN.shp')

# As I haven't downloaded SRTM data for north-eastern Morocco, I'm going to remove the polygons in this region.
b <- as(extent(5, 10, 33,37.5), 'SpatialPolygons')
crs(b) <- crs(PinpinDistri)
PinpinDistri <- rgeos::gDifference(PinpinDistri,b)

# To see the mask:
# plot(PinpinDistri)

# saving the mask:
shapefile(PinpinDistri, filename='data/maps/MaskPinpinDistri/PinpinDistriEuforgen.shp', overwrite=TRUE)


# We are going to expand the mask with the NFI plots
#####################################################"

# Load the data from Alexandre Changenet with the plot coordinates
data <- readRDS(file="data/IFN/alexdata/dfplotPINPINA0.8R.M.rds")

# Rm plots where there is no maritime pine
data <- data[!(data$treeNbrJ.IMall==0),]

# From geographic coordinates to spatial points:
NFIplots <- SpatialPoints(coords=data[,c("longitude","latitude")],proj4string=crs(PinpinDistri))

# Buffer with a radius of 10km
NFIplots <- buffer(NFIplots,width=10000)

# I don't know why byt buffering changes the projection, so we fix it:
NFIplots <- spTransform(NFIplots,crs(PinpinDistri))

# Country boundaries (to avoid that parts of the buffer are in the ocean...)
borders <- shapefile("../../Pinpin_Clonapin/maps/countries_boundaries/ne_50m_admin_0_countries.shp")

# Intersection btw countries and buffer
NFIplots <- rgeos::gIntersection(NFIplots,borders)

# Union between NFI plots and maritime pine distribution from Euforgen
mask <- rgeos::gUnion(PinpinDistri,NFIplots)

plot(mask,col="green")
plot(borders,add=T)


# saving the mask:
shapefile(mask, filename='data/maps/MaskPinpinDistri/PinpinDistriEUforgen_NFIplotsBuffer10km.shp', overwrite=TRUE)
