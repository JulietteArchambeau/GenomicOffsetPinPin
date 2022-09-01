library(tidyverse)
library(tidybayes)
library(sp)
library(raster)
library(sf)
library(rnaturalearth)
library(ggthemes)
library(cowplot)

# Figure of the validation in natural populations
#################################################"


DF.go <- readRDS(file="outputs_paper/ValidationNFI/OutputsModels.rds")

model.labeller <- c("Generalised Dissimilarity Modeling", "Gradient Forests")
names(model.labeller) <- c("GDM", "GF")

p <- DF.go %>% 
    filter(grepl("betaGO",term)) %>% 
    filter(!grepl("585",GenomicOffset)) %>% 
    mutate(Model=str_sub(GenomicOffset,18,-1),
           SNPSets=str_sub(GenomicOffset,14,16),
           SNPSets=factor(SNPSets,levels=c("Com","Mer","Ref")),
           Country=case_when(term=="betaGO_country[1]" ~"Spain",
                             term=="betaGO_country[2]" ~ "France"))  %>%
    
    ggplot(aes(x = Country, y = estimate,ymin = conf.low, ymax = conf.high,colour=SNPSets,shape=SNPSets)) +
    geom_pointinterval(position = position_dodge(width = .6),point_size=2.5,size=3,show.legend = c(size = TRUE)) +
    geom_hline(yintercept = 0,color="gray") +
    facet_grid(.~Model,scales="free", space = "free",labeller = labeller(Model = model.labeller)) + 
    ylab(TeX("Regression coefficients")) + xlab("") +
    scale_colour_manual(values=c("forestgreen","chartreuse2","gold1"),labels=c("Common candidate SNPs",
                                                                               "Merged candidate SNPs",
                                                                               "Reference SNPs"),name="SNP sets") +
    scale_shape_manual(values=c(16,17,8),labels=c("Common candidate SNPs",
                                                  "Merged candidate SNPs",
                                                  "Reference SNPs"),name="SNP sets") +
    #scale_x_discrete(labels=c("SSP370" = "SSP3-7.0", "SSP585" = "SSP5-8.5")) +
    theme_bw() +
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=16),
          legend.title=element_text(size=10), 
          legend.text=element_text(size=9),
          #legend.position = c(.71, .87),
          legend.background = element_rect(colour = "grey"),
          strip.text.x = element_text(size = 14),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    guides(color=guide_legend(ncol=1))
  


ggsave(p, file="figs/PosterSSMPG2022/ValidationNFI.pdf",device="pdf",height=4,width=10)




# Map of the experimental design 
################################"

pops <- read.csv("data_DRYAD/PopulationPopStructureEnvData.csv") %>% 
  dplyr::select(prov,latitude,longitude) %>% 
  unique()

# Create a spatial object of the provenance coordinates 
xyprov <- SpatialPoints(pops[,c("longitude","latitude")], 
                        proj4string=CRS("+proj=longlat +datum=WGS84"))

# for mapping
xyprovdf <- xyprov@coords %>% as_tibble()


sites <- readRDS("data/AllDataPhenoClimSoil.RDS") %>% 
  dplyr::select(site,longitude_site,latitude_site) %>% 
  unique() %>% 
  mutate(site=str_to_title(site),
         site=replace(site,site=="Caceres","Cáceres"))

# Create a spatial object of the site coordinates 
xysite <- SpatialPoints(sites[,c("longitude_site","latitude_site")], 
                        proj4string=CRS("+proj=longlat +datum=WGS84"))

# for mapping
xysitedf <- xysite@coords %>% as_tibble()


# Create a buffer of 10 km around the population and site location
buffer.prov <- buffer(xyprov,width=10000)
buffer.site <- buffer(xysite,width=10000)

PinpinDistri  <- shapefile('../../GenomicOffset/GenomicOffsetPinPin/data/maps/MaskPinpinDistri/PinpinDistriEUforgen_NFIplotsBuffer10km.shp')
PinpinDistri <- spTransform(PinpinDistri,crs(xysite))

# Merge maritime pine distribution and 10-km buffer around the population and site location
PinpinDistri <- rgeos::gUnion(PinpinDistri,buffer.prov)
PinpinDistri <- rgeos::gUnion(PinpinDistri,buffer.site)

PinpinDistri <- st_as_sf(PinpinDistri)

world <- ne_countries(scale = "medium", returnclass = "sf")



p <- ggplot() +
  geom_sf(data = world,color="bisque3",size=0.4,fill=alpha("aliceblue",0.6))  +
  geom_sf(data = PinpinDistri,fill="burlywood",size=0.02,color="burlywood")  + # color=alpha("chocolate4",0.8),
  coord_sf(xlim = c(-10, 15), 
           ylim = c(33,51),
           crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") + 
  geom_point(data=xyprovdf,aes(x=longitude,y=latitude),size=2) +
  geom_point(data=sites[sites$site%in%c("Portugal","Asturias","Bordeaux"),],aes(x=longitude_site,y=latitude_site),col="magenta",size=3,shape=17) +
  geom_point(data=sites[!sites$site%in%c("Portugal","Asturias","Bordeaux"),],aes(x=longitude_site,y=latitude_site),col="orangered",size=3,shape=17) +
  # geom_text(data=sites,aes(x=longitude_site,y=latitude_site,label=site),hjust=1.2, vjust=-1.2, col="red",size=6,fontface="bold") +
  theme_map()


p
ggsave(p,file="maps/map_posterSSMPG2022.svg")




# ff
#####"

labels <- c("Common candidate SNPs",
            "Merged candidate SNPs",
            "Reference SNPs")

p <- readRDS(file=paste0("outputs_paper/ValidationCG/MortalityModels/TableCoeffMortalityModels.rds")) %>% 
  filter(term=="betaX1") %>% 
  filter(!str_detect(variable,"CTD")) %>% 
  mutate(Method.GO=case_when(str_detect(variable,"GDM")==TRUE~"GDM",
                             str_detect(variable,"GF")==TRUE~"GF"),
         SNPset=case_when(str_detect(variable,"GO")==TRUE~str_sub(variable,4,6)),
         Method.GO=factor(Method.GO,levels=c("GDM","GF")),
         SNPset=factor(SNPset,levels=c("Com","Mer","Ref",paste0("bio",c(1,5,6,12,15)))),
         site=case_when(site=="caceres"~"Cáceres",
                        site=="madrid"~"Madrid")) %>% 
  
  ggplot(aes(x = Method.GO, y = estimate,ymin = conf.low, ymax = conf.high,color=SNPset,shape=SNPset)) +#
  geom_hline(yintercept = 0,color="gray") +
  geom_pointinterval(position = position_dodge(width = .4),point_size=3.5,size=3) + # 
  facet_grid(.~site,scales="free_x", space = "free") + 
  ylab("") + xlab("") +
  scale_color_manual(values=c("forestgreen","chartreuse2","gold1"),labels=labels,name="SNP sets") +
  scale_shape_manual(values = c(16,17,8),labels=labels,name="SNP sets") +
  theme_bw() +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=1),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=15),
        legend.background = element_rect(colour = "grey"),
        strip.text.x = element_text(size = 16.5),
        strip.text.y = element_text(size=16,angle = 0),
        strip.background = element_rect(fill = alpha("orangered",0.5)),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) +
  guides(color=guide_legend(ncol=1),
         shape = guide_legend(override.aes = list(size =2 )))

part2.fig.paper <- p

DF <- readRDS(file="outputs_paper/ValidationCG/HeightModels/PosteriorsHeightModels.rds")

pdf <- DF %>% 
  filter(term=="betaX1") %>% 
  filter(!str_detect(variable,"CTD")) %>% 
  filter(height.m %in% c("POR_htmay13","AST_htmar14","BDX_htnov18")) %>% 
  mutate(Method.GO=case_when(str_detect(variable,"GDM")==TRUE~"GDM",
                             str_detect(variable,"GF")==TRUE~"GF"),
         Method.GO=factor(Method.GO,levels=c("GDM","GF","CTD")),
         SNPset=case_when(str_detect(variable,"GO")==TRUE~str_sub(variable,4,6)),
         SNPset=factor(SNPset,levels=c("Com","Mer","Ref",paste0("bio",c(1,5,6,12,15)))),
         Site.Age.labels=case_when(height.m=="BDX_htnov18"~ "Bordeaux",
                                   height.m=="POR_htmay13"~ "Portugal",
                                   height.m=="AST_htmar14"~ "Asturias")) 

p <- pdf %>% 
  ggplot(aes(x = Method.GO, y = estimate,ymin = conf.low, ymax = conf.high,colour=SNPset,shape=SNPset)) +#
  geom_hline(yintercept = 0,color="gray") +
  geom_pointinterval(position = position_dodge(width = .4),point_size=3.5,size=3,show.legend = c(size = TRUE)) +
  facet_grid(.~Site.Age.labels,scales="free_x", space = "free") + 
  ylab("") + xlab("") +
  scale_color_manual(values=c("forestgreen","chartreuse2","gold1"),labels=labels) +
  scale_shape_manual(values = c(16,17,8),labels=labels) +
  theme_bw() +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        legend.title=element_text(size=13), 
        legend.text=element_text(size=12),
        legend.background = element_rect(colour = "grey"),
        legend.position = "none",
        strip.text.x = element_text(size = 15.5),
        strip.text.y = element_text(size=14,angle = 0),
        strip.background = element_rect(fill = alpha("magenta",0.4)),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) +
  guides(color=guide_legend(ncol=2))


grid.graphs <- plot_grid(p,part2.fig.paper,nrow=2)


ggsave(grid.graphs, file="figs/PosterSSMPG2022/ValidationCG.pdf",device="pdf",height=8,width=10)
