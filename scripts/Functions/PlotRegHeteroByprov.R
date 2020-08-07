# Plotting the regional climatic heterogeneity by provenances (grouped by gene pools)
PlotRegHeteroByprov <- function(data){
  ggparcoord(data,
             columns = (length(colnames(data))-5):length(data), 
             groupColumn = "mainGP",
             scale="globalminmax",
             showPoints = TRUE, 
             alphaLines = 0.3) + 
    scale_color_manual(values=c("orangered3",
                                "gold2",
                                "darkorchid3",
                                "navyblue",
                                "turquoise2",
                                "green3"), labels = c("Northern Africa", 
                                                      "Corsica",
                                                      "Central Spain",
                                                      "French Atlantic",
                                                      "Iberian Atlantic",
                                                      "South-eastern Spain"),name="Gene pools") +
    geom_line(size=1,alpha=0.5) +
    labs(y="Regional heterogeneity",x="Buffer radius (in km)") +
    theme_bw() +
    theme(plot.title = element_text(size=10))
  
}