# Plotting the regional climatic heterogeneity by radius of the buffer zones
PlotRegHeteroByradius <- function(data){
  data <- data[with(data, order(data$mainGP)),]
  row.names(data) <- data$prov
  data <- data[,c("5","10","20","50","75","100")]
  val.max <- max(data)
  data <- as.data.frame(t(data))
  data[,"radius"] <- factor(row.names(data), levels=c("5","10","20","50","75","100"))
  
  ymax <- -(val.max/15) +(val.max/20) 
  ymin <- -(val.max/15)
  
  rectNA <- data.frame (xmin=0, xmax=2.5, ymin=ymin, ymax=ymax)
  rectC <- data.frame (xmin=2.5, xmax=4.5, ymin=ymin, ymax=ymax)
  rectCS <- data.frame (xmin=4.5, xmax=16.5, ymin=ymin, ymax=ymax)
  rectFA <- data.frame (xmin=16.5, xmax=23.5, ymin=ymin, ymax=ymax)
  rectIA <- data.frame (xmin=23.5, xmax=32.5, ymin=ymin, ymax=ymax)
  rectSES <- data.frame (xmin=32.5, xmax=35, ymin=ymin, ymax=ymax)
  
  mypalette <- c(brewer.pal(n=11,"PiYG")[1:3],brewer.pal(n=11,"PiYG")[9:11])
  
  ggparcoord(data,
             columns = 1:34, groupColumn = "radius",
             scale="globalminmax",
             showPoints = TRUE, 
             alphaLines = 0.8
  ) +
    geom_rect(data=rectNA, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="orangered3", alpha=0.8, inherit.aes = FALSE) +
    geom_rect(data=rectC, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gold2", alpha=0.8, inherit.aes = FALSE) +
    geom_rect(data=rectCS, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="darkorchid3", alpha=0.8, inherit.aes = FALSE) +
    geom_rect(data=rectFA, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="navyblue", alpha=0.8, inherit.aes = FALSE) +
    geom_rect(data=rectIA, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="turquoise2", alpha=0.8, inherit.aes = FALSE) +
    geom_rect(data=rectSES, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="green3", alpha=0.8, inherit.aes = FALSE) +
    geom_line(size=2,alpha=0.5) + 
    #scale_color_viridis(discrete=TRUE) +
    coord_cartesian(ylim=c(0,val.max)) +
    #scale_color_grey(start = 0.8, end = 0.01) +
    scale_color_manual(values = mypalette,labels = c("5 km", 
                                                     "10 km",
                                                     "20 km",
                                                     "50 km",
                                                     "75 km",
                                                     "100 km"),name="Buffer radius") + 
    theme_bw()+
    labs(y="Regional heterogeneity",x="Provenances") +
    theme(plot.title = element_text(size=10))
  
}