# Mapping with leaflet of the climatic variables and the buffer zones around each provenance
MapClimRegHeteroInteractive <- function(clim.var.name,clim.var,clim.var.raster,reverse){
  
  pal=colorNumeric("Spectral",raster::values(clim.var.raster),na.color = "transparent", reverse = reverse)
  
  leaflet() %>% 
    addRasterImage(clim.var.raster , colors = pal, opacity = 0.8) %>%
    addLegend("bottomright",title=clim.var.name,pal=pal,values=raster::values(clim.var.raster)) %>% 
    addPolygons(data=buff5,fillColor = "transparent",color = "#8E0152") %>% 
    addPolygons(data=buff10,fillColor = "transparent",color="#C51B7D") %>% 
    addPolygons(data=buff20,fillColor = "transparent",color = "#DE77AE") %>% 
    addPolygons(data=buff50,fillColor = "transparent",color="#7FBC41") %>% 
    addPolygons(data=buff75,fillColor = "transparent",color = "#4D9221") %>% 
    addPolygons(data=buff100,fillColor = "transparent",color="#276419") %>% 
    #addCircleMarkers(data=df,~longitude, ~latitude,label = ~htmlEscape(prov),radius = 0.01,color="black",fillOpacity = 1)
    addCircleMarkers(data=df,~longitude, ~latitude,radius = 0.01,color="black",fillOpacity = 1,
                     label=labels,labelOptions = labelOptions(style = list("font-weight" = "normal", padding = "3px 8px"),
                                                              textsize = "15px", direction = "auto"))
}