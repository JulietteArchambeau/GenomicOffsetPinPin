################################################
# Mask of the maritime distribution for the maps
################################################

# Maritime pine distribution from EUFORGEN (shapefile)
PinpinDistri  <- shapefile('../../Pinpin_Clonapin/maps/pinus_pinaster_distribution/Pinus_pinaster_EUFORGEN.shp')

# As I haven't downloaded SRTM data for north-eastern Morocco, I'm going to remove the polygons in this region.
b <- as(extent(5, 10, 33,37.5), 'SpatialPolygons')
crs(b) <- crs(PinpinDistri)
PinpinDistri <- gDifference(PinpinDistri,b)

#To see the mask:
#plot(PinpinDistri)

# saving the mask:
shapefile(PinpinDistri, filename='data/maps/MaskPinpinDistri/PinpinDistri.shp', overwrite=TRUE)
