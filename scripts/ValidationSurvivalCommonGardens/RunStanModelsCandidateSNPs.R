
# Validating the candidate SNPs in dry environments
###################################################"

# Juliette Archambeau
# 22/12/2020

# This script has to be run from bash because it crashes in RStudio

# libraries
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Parameters to select:
site="caceres"
covar <- c("bio12","bio14") # all variables significant in Caceres
#covar <- c("bio2","bio12","bio14")
#covar <- c("bio2","bio5","bio6","bio12","bio13","bio14") # all avariables significant in Madrid
SNPset <- c("MergedCand","ComCand","IntCand")
CovSet <- "AvgWater"



# Load the height data
HierMod <- readRDS(file="../../Pinpin_Clonapin/HeightPinpinClonapin/outputs/models/P1/MOD1.rds")
heights <- HierMod %>% broom::tidyMCMC(estimate.method = "mean",conf.int = T) %>% # we take the mean of the prov random intercepts
  filter(str_detect(term, "^(r_prov\\[)")) %>% 
  dplyr::rename(height=estimate,prov=term) %>% 
  mutate(prov=str_sub(prov,8,-12))

# Load the survival data in the common gardens
data <- readRDS(file="data/AllDataPhenoClimSoil.RDS")

# Removing genotypes for which there is no genomic data
data <- data[!(is.na(data$Q1)),]

# Load the genomic dataset
geno <- readRDS("data/GenomicData/5165snps523genotypesImpNA.rds")

# Dataframe in which we are going to stock the results
DF <- crossing(covar,SNPset) %>% mutate(NbSiginficantBeta1=NA,NbCandidates=NA,Beta1SNPcount=NA,CrossZeroSNPcount=NA)

stancode = stan_model("scripts/StanModels/BinomialValidationCommonGardens/BinSurvivalPropSNPAndHeight.stan")

#### START LOOP FOR EACH COVARIATE
for(i in covar){ 
  
  # Dataframe of the survival counts in the common garden of interest
  subdata <- data[data$site==site,] %>% 
    select(prov,survival) %>% 
    group_by(prov) %>% 
    dplyr::summarise(survival=sum(survival),tot.count=n()) %>% 
    left_join(heights[,1:2],by="prov")
  
  # Add the covariate values to the dataframe
  subdata <- data %>% 
    select(prov,paste0(i,"_prov")) %>% 
    dplyr::rename(covariate=paste0(i,"_prov")) %>% 
    unique() %>% 
    left_join(subdata,by="prov")
  
  
  #### START LOOP FOR SET OF SNPs
  for(j in SNPset){
    
    # Load the candidate SNPs (we are only comparing the distribution of candidate SNPs here)
    cand <- try(readRDS(file=paste0("outputs/CandidateSNPs/SetsCandidateSNPs/CovariateSpecific/",i,"_",CovSet,"_",j,".rds")),silent=T)
    
    if(!inherits(cand,"try-error")){ # some sets are empty and there are no corresponding file, so this line allows the loop not to stop if the file is not found
      
      # Number of selected candidate SNPs
      #length(cand$snp)
      
      # Keep only the candidate SNPs
      genocand <- geno[,colnames(geno) %in% cand$snp] %>% 
        rownames_to_column(var="genotype") %>% 
        dplyr::rename(prov=genotype) %>% 
        mutate(prov=str_sub(prov,1,3)) 
      
      candprop <- genocand %>% 
        group_by(prov) %>% 
        summarise_all(~sum(., na.rm = TRUE)/(n()*2))
      
      # alleles with a 
      cand$corr <- NA
      
      for(c in cand$snp){cand[cand$snp==c,"corr"] <- cor(candprop[,c],subdata$covariate)[,1]}
      
      # select the alleles that have to be inverted
      if(i %in% c("bio12","bio14")){
        invert.alleles <- cand$snp[cand$corr>0] # for precipitation variables, the beneficial alleles in Madrid and Caceres have a negative relationship with the precipitation-related covariate, and so the alleles that have to be inverted are those with a positive relationship with the covariate.
      } else {
        invert.alleles <- cand$snp[cand$corr<0] # for temperature covariates, the beneficial alleles in Madrid and Caceres have a positive relationship with the temperature-related covariate.
      }
      
      # Invert the alleles:
      genocand[,invert.alleles] <- genocand %>% select(contains(invert.alleles)) %>% mutate_all(~recode(., `0` = 2,`2`=0))
      
      
      # SNP-by-SNP ANALYSIS
      # ===================#
      
      
      # Recalculate the proportion and join with the dataset with survivors and height
      candprop <- genocand %>% 
        group_by(prov) %>% 
        summarise_all(~sum(., na.rm = TRUE)/(n()*2)) %>%  
        left_join(subdata,by="prov")
      
      output <- tibble(snp=colnames(genocand[,-1]),
                       beta1=NA,
                       LowCI1=NA,
                       HighCI1=NA,
                       CrossZero1=NA,
                       beta2=NA,
                       LowCI2=NA,
                       HighCI2=NA,
                       CrossZero2=NA)
      
      #### START LOOP FOR SNP
      for(c in colnames(genocand[,-1])){
        datalist <- list(N=length(candprop$prov),
                         NbTot=candprop %>% pull(tot.count),
                         NbSurv=candprop %>% pull(survival),
                         heightSC=scale(candprop %>% pull(height))[,1],
                         propsnpSC= scale(candprop %>% pull(all_of(c)))[,1])
        
        # Model 1
        mstan <- sampling(stancode, data = datalist, iter = 2000, chains = 4, cores = 4,init=0)
        
        df <- broom::tidyMCMC(mstan,droppars = NULL, estimate.method = "median", conf.int = T,conf.level = 0.95) 
        
        # SNP
        output[output$snp==c,"beta1"] <- df[df$term=="betaPropSNP","estimate"]$estimate
        output[output$snp==c,"LowCI1"] <- df[df$term=="betaPropSNP","conf.low"]$conf.low
        output[output$snp==c,"HighCI1"] <- df[df$term=="betaPropSNP","conf.high"]$conf.high
        output[output$snp==c,"CrossZero1"] <- ifelse((output$LowCI1[output$snp==c]>0 & output$HighCI1[output$snp==c]<0) | (output$LowCI1[output$snp==c]<0 & output$HighCI1[output$snp==c]>0), "-" ,"*")
        
        # Height
        output[output$snp==c,"beta2"] <- df[df$term=="betaHeight","estimate"]$estimate
        output[output$snp==c,"LowCI2"] <- df[df$term=="betaHeight","conf.low"]$conf.low
        output[output$snp==c,"HighCI2"] <- df[df$term=="betaHeight","conf.high"]$conf.high
        output[output$snp==c,"CrossZero2"] <- ifelse((output$LowCI2[output$snp==c]>0 & output$HighCI2[output$snp==c]<0) | (output$LowCI2[output$snp==c]<0 & output$HighCI2[output$snp==c]>0), "-" ,"*")
        
      }
      
      saveRDS(output,file=paste0("outputs/ValidationSurvivalCommonGardens/CandidateSNPs/SNPbySNPanalyses_",i,"_",CovSet,"_",j,"_",site,".rds"))
      DF[DF$covar==i&DF$SNPset==j,"NbCandidates"] <- length(cand$snp)
      DF[DF$covar==i&DF$SNPset==j,"NbSiginficantBeta1"] <- length(output$snp[output$beta1>0&output$CrossZero1=="*"])
      
    
     # ANALYSIS based on counts of SNPs
     # ===============================#
      
      genocand <-genocand %>%
      group_by(prov) %>% 
      summarise_all(~sum(., na.rm = TRUE)/(n()*2)) %>% 
      mutate(mProp= rowMeans(select(.,starts_with("SNP")), na.rm = TRUE)) %>% 
      select(prov,mProp) %>%  
      left_join(subdata,by="prov")
    
    output <- tibble(beta1=NA,
                     LowCI1=NA,
                     HighCI1=NA,
                     CrossZero1=NA,
                     beta2=NA,
                     LowCI2=NA,
                     HighCI2=NA,
                     CrossZero2=NA)
    
    datalist <- list(N=length(genocand$prov),
                     NbTot=genocand %>% pull(tot.count),
                     NbSurv=genocand %>% pull(survival),
                     heightSC=scale(genocand %>% pull(height))[,1],
                     propsnpSC= scale(genocand %>% pull(mProp))[,1])
    
    # Model 1
    mstan <- sampling(stancode, data = datalist, iter = 2000, chains = 4, cores = 4,init=0)
    
    df <- broom::tidyMCMC(mstan,droppars = NULL, estimate.method = "median", conf.int = T,conf.level = 0.95) 
    
    # SNP
    DF[DF$covar==i&DF$SNPset==j,"Beta1SNPcount"] <- output[1,"beta1"] <- df[df$term=="betaPropSNP","estimate"]$estimate
    output[1,"LowCI1"] <- df[df$term=="betaPropSNP","conf.low"]$conf.low
    output[1,"HighCI1"] <- df[df$term=="betaPropSNP","conf.high"]$conf.high
    DF[DF$covar==i&DF$SNPset==j,"CrossZeroSNPcount"] <-  output[1,"CrossZero1"] <- ifelse((output$LowCI1[1]>0 & output$HighCI1[1]<0) | (output$LowCI1[1]<0 & output$HighCI1[1]>0), "-" ,"*")
    
    # Height
    output[1,"beta2"] <- df[df$term=="betaHeight","estimate"]$estimate
    output[1,"LowCI2"] <- df[df$term=="betaHeight","conf.low"]$conf.low
    output[1,"HighCI2"] <- df[df$term=="betaHeight","conf.high"]$conf.high
    output[1,"CrossZero2"] <- ifelse((output$LowCI2[1]>0 & output$HighCI2[1]<0) | (output$LowCI2[1]<0 & output$HighCI2[1]>0), "-" ,"*")
    
    saveRDS(output,file=paste0("outputs/ValidationSurvivalCommonGardens/CandidateSNPs/SNPCountanalysis_",i,"_",CovSet,"_",j,"_",site,".rds"))
    
    }
    
    
    }}

saveRDS(DF,file=paste0("outputs/ValidationSurvivalCommonGardens/CandidateSNPs/DF_",site,".rds"))

