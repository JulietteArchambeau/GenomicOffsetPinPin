########################################################################################################"
##################                                                               #######################"
##################                     BLUPs for height                          #######################"
##################                  with the package "brms"                      #######################"
##################                                                               #######################"
########################################################################################################"

library(dplyr)
library(tidyverse)
library(brms)
library(bayesplot)
options(mc.cores = parallel::detectCores())
library(xtable)

# Options for the posteriors:
prob=0.95
prob_outer = 0.99
point_est = "mean"

# load the phenotypic data
data <- readRDS(file="data/ClonapinData/PhenoDataNovember2019_AnnualTraits_UpdatedSept2021_AllSites.rds")

# Select the sites and years of measurements in each site
sites <- c("MAD_htdec11","CAC_htdec11","BDX_htnov13","BDX_htnov18","POR_htjan12","POR_htmay13","AST_htdec11","AST_htmar14")


# Experimental design for table in the Supplementary Information:
# --------------------------------------------------------------
ExpDesignTab <- sapply(sites,function(x){
data %>%
  dplyr::select(prov,all_of(x)) %>% 
  dplyr::rename(height=all_of(x)) %>% 
  drop_na() %>% 
  group_by(prov) %>% 
  dplyr::summarise(count =n(),
                   mean.height=mean(height), .groups = 'drop')

  
  
}, USE.NAMES = TRUE,simplify=FALSE ) %>% 
  bind_rows(.id="variable") %>% 
  pivot_wider(names_from=variable,values_from = c(count,mean.height),names_sep="_") %>% 
  dplyr::select(prov,contains("AST_htdec"),contains("AST_htmar"),
                contains("BDX_htnov13"),contains("BDX_htnov18"),
                contains("CAC"),
                contains("MAD"),
                contains("POR_htjan"),contains("POR_htmay"))
  
# Generate the latex table
print(xtable(ExpDesignTab, type = "latex",digits=2), 
      file = paste0("outputs/ValidationCG/ExperimentalDesignTables/HeightData.tex"), 
      include.rownames=FALSE)



# Run the models:
# ---------------
lapply(sites,function(x){
  
subdata <- data %>%
  dplyr::select(prov,block,all_of(x),c(paste0("Q",rep(1:6)))) %>% 
  dplyr::rename(height=all_of(x)) %>% 
  drop_na() %>% 
  mutate(height=(height-mean(height))/sd(height)) #  ! height is centered !

# Population structure
colnames(subdata)[colnames(subdata) %in% c(paste0("Q",rep(1:6)))] <- c(paste0("prop_Q",rep(1:6)))
head(subdata[,c(paste0("prop_Q",rep(1:6)))])
subdata$Q1 <- "Q1"
subdata$Q2 <- "Q2"
subdata$Q3 <- "Q3"
subdata$Q4 <- "Q4"
subdata$Q5 <- "Q5"
subdata$Q6 <- "Q6"
sum(subdata[,c(paste0("prop_Q",rep(1:6)))]<0)
filter_at(subdata,c(paste0("prop_Q",rep(1:6))),any_vars(. < 0))
subdata$prop_Q6[subdata$prop_Q6<0] <- 0


# Specify the priors
priors <- c(prior(normal(0, 10), class = Intercept),
            prior(cauchy(0, 10), class = sd),
            prior(cauchy(0, 10), class = sigma))



# Model not accounting for population structure
# ---------------------------------------------

mod <- brm(height  ~  (1|prov) +  (1|block) ,
           
           prior = priors,
           
           data = subdata, 
           family = "gaussian",
           chain=4,
           iter=2500)

saveRDS(mod,file=paste0("outputs/ValidationCG/ModelsHeightBLUPs/NotAccountingPopStructure/Models/",x,".rds"))

p <- readRDS(file=paste0("outputs/ValidationCG/ModelsHeightBLUPs/NotAccountingPopStructure/Models/",x,".rds")) %>% 
  mcmc_intervals(regex_pars = "^r_prov\\[",
                 prob=prob,
                 prob_outer=prob_outer,
                 point_est = point_est) +  
  theme_bw() + 
  theme(axis.text = element_text(size=16))
ggsave(p,file=paste0("outputs/ValidationCG/ModelsHeightBLUPs/NotAccountingPopStructure/FigsPosteriors/ProvIntercepts_",x,".pdf"),device="pdf")



# Model accounting for population structure
# -----------------------------------------

mod <- brm(height  ~  (1|prov) + 
                      (1|block) +
                      (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))),
            
           prior = priors,
            
           data = subdata, 
           family = "gaussian",
           #control = list(adapt_delta=0.999,max_treedepth =14),
           chain=4,
           iter=2500)

saveRDS(mod,file=paste0("outputs/ValidationCG/ModelsHeightBLUPs/AccountingPopStructure/Models/",x,".rds"))

p <- readRDS(file=paste0("outputs/ValidationCG/ModelsHeightBLUPs/AccountingPopStructure/Models/",x,".rds")) %>% 
  mcmc_intervals(regex_pars = "^r_prov\\[",
                 prob=prob,
                 prob_outer=prob_outer,
                 point_est = point_est) +  
  theme_bw() + 
  theme(axis.text = element_text(size=16))
ggsave(p,file=paste0("outputs/ValidationCG/ModelsHeightBLUPs/AccountingPopStructure/FigsPosteriors/ProvIntercepts_",x,".pdf"),device="pdf")

p <- readRDS(file=paste0("outputs/ValidationCG/ModelsHeightBLUPs/AccountingPopStructure/Models/",x,".rds")) %>% 
  mcmc_intervals(regex_pars = "^r_mm",
                 prob=prob,
                 prob_outer=prob_outer,
                 point_est = point_est) +  
  theme_bw() + 
  theme(axis.text = element_text(size=16))

ggsave(p,file=paste0("outputs/ValidationCG/ModelsHeightBLUPs/AccountingPopStructure/FigsPosteriors/GenePoolIntercepts_",x,".pdf"),device="pdf")

})



