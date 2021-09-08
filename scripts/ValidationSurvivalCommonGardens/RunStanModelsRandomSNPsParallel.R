
# Validating of the candidate SNPs in dry environments
######################################################"

# Part with the random SNPs
# ----------------------------#

# Juliette Archambeau
# 22/12/2020

# This script has to be run from bash because it crashes in RStudio

setwd('~/Documents/GenomicOffset/GenomicOffsetPinPin/')

start_time <- Sys.time()

# libraries
library(tidyverse)
library(rstan)
library(reshape)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(foreach)
library(doParallel)



# Parameters to select:
site="caceres"
#covar <- c("bio12","bio14") # all variables significant in Caceres
covar <- c("bio14")
#covar <- c("bio2","bio12","bio14")
#covar <- c("bio2","bio5","bio6","bio12","bio13","bio14") # all avariables significant in Madrid
#SNPset <- c("MergedCand","ComCand","IntCand")
SNPset <- c("IntCand")
CovSet <- "AvgWater"
nsample <- 1000



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

stancode = stan_model("scripts/StanModels/BinomialValidationCommonGardens/BinSurvivalPropSNPAndHeight.stan")


# Create a list with the frequency of the random SNPs
cl <- parallel::makeCluster(8) # Options for parallelization
doParallel::registerDoParallel(cl) # starts parallelization
listRand <- foreach(i = covar, .packages=c("reshape","tidyverse")) %dopar% { 

  listRand <- list(NA,NA,NA)
  names(listRand) <- SNPset
  
  for(j in SNPset){
    
    # Load the candidate SNPs (we are only comparing the distribution of candidate SNPs here)
    cand <- try(readRDS(file=paste0("outputs/CandidateSNPs/SetsCandidateSNPs/CovariateSpecific/",i,"_",CovSet,"_",j,".rds")),silent=T)
    
    if(!inherits(cand,"try-error")){ # some sets are empty and there are no corresponding file, so this line allows the loop not to stop if the file is not found
      
      listRand[[j]] <- geno[,!(colnames(geno) %in% cand$snp)]%>% 
        summarise_all(~sum(., na.rm = TRUE)/(n()*2)) %>%  
        melt() %>% 
        dplyr::rename(snp.rand=variable,freq.rand=value)
      
    }}
    
    listRand}

parallel::stopCluster(cl) 
names(listRand) <- covar  


# Create a list with the frequency of the candidate SNPs
cl <- parallel::makeCluster(8) # Options for parallelization
doParallel::registerDoParallel(cl)
listCand <- foreach(i = covar, .packages=c("reshape","tidyverse")) %dopar% { 
  
  listCand <- list(NA,NA,NA)
  names(listCand) <- SNPset
  
  for(j in SNPset){
    
    # Load the candidate SNPs (we are only comparing the distribution of candidate SNPs here)
    cand <- try(readRDS(file=paste0("outputs/CandidateSNPs/SetsCandidateSNPs/CovariateSpecific/",i,"_",CovSet,"_",j,".rds")),silent=T)
    
    if(!inherits(cand,"try-error")){ # some sets are empty and there are no corresponding file, so this line allows the loop not to stop if the file is not found
      
      # Keep only the covariate-associated candidate SNPs in the genomic dataset
      listCand[[j]] <- geno[,colnames(geno) %in% cand$snp]%>%
        summarise_all(~sum(., na.rm = TRUE)/(n()*2))  %>%
        melt() %>%
        dplyr::rename(snp.cand=variable,freq.cand=value)
      
    }}
  
  listCand}

parallel::stopCluster(cl) 
names(listCand) <- covar  


########################################################################################################################################################"
cl <- parallel::makeCluster(8) # Options for parallelization
doParallel::registerDoParallel(cl) # starts parallelization


#### START LOOP FOR EACH COVARIATE
DF <- foreach(i = covar,.packages = c("rstan","tidyverse","reshape")) %dopar% { 
  
  # Dataframe of the survival counts in the common garden of interest
  subdata <- data[data$site==site,] %>% 
    select(prov,survival) %>% 
    group_by(prov) %>% 
    dplyr::summarise(survival=sum(survival),tot.count=n()) %>% 
    left_join(heights[,1:2],by="prov")
  
  # Add the covariate values to the dataframe
  subdata <- data %>% 
    dplyr::select(prov,paste0(i,"_prov")) %>% 
    dplyr::rename(covariate=paste0(i,"_prov")) %>% 
    unique() %>% 
    left_join(subdata,by="prov")
  
  
  #### START LOOP FOR SET OF SNPs (merged candidates, intermediate candidates, common candidates)
  for(j in SNPset){
    
    # Load the candidate SNPs (we are only comparing the distribution of candidate SNPs here)
    cand <- try(readRDS(file=paste0("outputs/CandidateSNPs/SetsCandidateSNPs/CovariateSpecific/",i,"_",CovSet,"_",j,".rds")),silent=T)
    
    if(!inherits(cand,"try-error")){ # some sets are empty and there are no corresponding file, so this line allows the loop not to stop if the file is not found
      
      # Keep only the covariate-associated candidate SNPs in the genomic dataset
      freqcand <- listCand[[i]][[j]]
      
      freqgeno <- listRand[[i]][[j]]
      
      DF <- data.frame(sample.id=1:nsample,
                       NbDup1=NA,
                       NbDup099=NA,
                       NbRand=NA,
                       NbSigniRand=NA,
                       Beta1SNPcount=NA,
                       LowCI1SNPcount=NA,
                       HighCI1SNPcount=NA,
                       CrossZero1SNPcount=NA)
      
      
      #### START LOOP TO DRAW n samples of random SNPs
      for(n in 1:nsample){
        
        freqgenosub <- freqgeno
        
        ## Start here the selection of one group of random SNPs
        for(k in cand$snp){
          freqk <- freqcand[freqcand$snp.cand==k,"freq.cand"]
          sub <- freqgenosub[freqgenosub$freq.rand>(freqk-0.01)&freqgenosub$freq.rand<(freqk+0.01),]
          samp <- sample_n(sub,1) 
          freqcand[freqcand$snp.cand==k,"snp.rand"] <- samp[[1]]
          freqcand[freqcand$snp.cand==k,"freq.rand"] <- samp[[2]]
          freqgenosub <- freqgenosub[!(freqgenosub$snp.rand==samp[[1]]),]
        }
        

        # Keep only the selected SNPs in the genomic dataset
        genorand <- geno[,colnames(geno) %in% freqcand$snp.rand] %>%
          rownames_to_column(var="genotype") %>%
          dplyr::rename(prov=genotype) %>%
          mutate(prov=str_sub(prov,1,3))
        
        randprop <- genorand %>% 
          group_by(prov) %>% 
          summarise_all(~sum(., na.rm = TRUE)/(n()*2))
        
        # Invert the alleles
        freqcand$corr <- NA
        
        for(c in freqcand$snp.rand){freqcand[freqcand$snp.rand==c,"corr"] <- cor(randprop[,c],subdata$covariate)[,1]}
        
        if(i %in% c("bio12","bio14")){ # select the alleles that have to be inverted
          invert.alleles <- freqcand$snp.rand[freqcand$corr>0] # for precipitation variables, the beneficial alleles in Madrid and Caceres have a negative relationship with the precipitation-related covariate, and so the alleles that have to be inverted are those with a positive relationship with the covariate.
        } else {
          invert.alleles <- freqcand$snp.rand[freqcand$corr<0] # for temperature covariates, the beneficial alleles in Madrid and Caceres have a positive relationship with the temperature-related covariate.
        }
        genorand[,as.vector(invert.alleles)] <- genorand %>% dplyr::select(contains(as.vector(invert.alleles))) %>% mutate_all(~recode(., `0` = 2,`2`=0))
        
        
        ## Duplicates
        matcor <- cor(genorand[,-1], use = "pairwise.complete.obs")
        matcor[lower.tri(matcor,diag=T)] <- NA # Keep only the upper triangle of the matrix and remove values in the diagnonal
        DF[DF$sample.id==n,"NbDup1"] <- length(which(abs(matcor)==1)) # Counting the number of pairs of SNPs that have a correlation coefficient of 1
        DF[DF$sample.id==n,"NbDup099"] <-length(which(abs(matcor)>0.99)) # Counting the number of pairs of SNPs that have a correlation coefficient higher than 0.99
        
        
        # SNP-by-SNP ANALYSIS
        # ===================#
        
        # Recalculate the proportion and join with the dataset with survivors and height
        randprop <- genorand %>% 
          group_by(prov) %>% 
          summarise_all(~sum(., na.rm = TRUE)/(n()*2)) %>%  
          left_join(subdata,by="prov")

        output <- tibble(snp.rand=colnames(genorand[,-1]),
                         beta1=NA,
                         LowCI1=NA,
                         HighCI1=NA,
                         CrossZero1=NA,
                         beta2=NA,
                         LowCI2=NA,
                         HighCI2=NA,
                         CrossZero2=NA)
        
        
        for(c in colnames(genorand[,-1])){
          datalist <- list(N=length(randprop$prov),
                           NbTot=randprop %>% pull(tot.count),
                           NbSurv=randprop %>% pull(survival),
                           heightSC=scale(randprop %>% pull(height))[,1],
                           propsnpSC= scale(randprop %>% pull(all_of(c)))[,1])
          # Model 1
          mstan <- sampling(stancode, data = datalist, iter = 2000, chains = 4, cores = 4,init=0)
          
          df <- broom::tidyMCMC(mstan,droppars = NULL, estimate.method = "median", conf.int = T,conf.level = 0.95) 
          
          # SNP
          output[output$snp.rand==c,"beta1"] <- df[df$term=="betaPropSNP","estimate"]$estimate
          output[output$snp.rand==c,"LowCI1"] <- df[df$term=="betaPropSNP","conf.low"]$conf.low
          output[output$snp.rand==c,"HighCI1"] <- df[df$term=="betaPropSNP","conf.high"]$conf.high
          output[output$snp.rand==c,"CrossZero1"] <- ifelse((output$LowCI1[output$snp.rand==c]>0 & output$HighCI1[output$snp.rand==c]<0) | 
                                                              (output$LowCI1[output$snp.rand==c]<0 & output$HighCI1[output$snp.rand==c]>0), "-" ,"*")
          
          # Height
          output[output$snp.rand==c,"beta2"] <- df[df$term=="betaHeight","estimate"]$estimate
          output[output$snp.rand==c,"LowCI2"] <- df[df$term=="betaHeight","conf.low"]$conf.low
          output[output$snp.rand==c,"HighCI2"] <- df[df$term=="betaHeight","conf.high"]$conf.high
          output[output$snp.rand==c,"CrossZero2"] <- ifelse((output$LowCI2[output$snp.rand==c]>0 & output$HighCI2[output$snp.rand==c]<0) | 
                                                              (output$LowCI2[output$snp.rand==c]<0 & output$HighCI2[output$snp.rand==c]>0), "-" ,"*")
          
        }
        
        output <- left_join(freqcand,output,by="snp.rand") %>%
          as_tibble()
        
        # Complete DF 
        DF[DF$sample.id==n,"NbRand"] <- length(output$snp.rand) # total number of candidates (and random) SNPs
        DF[DF$sample.id==n,"NbSigniRand"] <- length(output$snp.rand[output$beta1>0&output$CrossZero1=="*"])  # Add the number coefficient beta1 that did not cross zero
        
        # saveRDS(output,file=paste0("outputs/ValidationInCommonGarden/ValidateAgainstRandomSNPs/SNPbySNPanalyses/sample",n,"_",site,".rds"))
        
        
        # ANALYSIS based on counts of SNPs
        # ===============================#
        
        genorand <-genorand %>%
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
        
        datalist <- list(N=length(genorand$prov),
                         NbTot=genorand %>% pull(tot.count),
                         NbSurv=genorand %>% pull(survival),
                         heightSC=scale(genorand %>% pull(height))[,1],
                         propsnpSC= scale(genorand %>% pull(mProp))[,1])
        
        mstan <- sampling(stancode, data = datalist, iter = 2000, chains = 4, cores = 4,init=0)
        
        df <- broom::tidyMCMC(mstan,droppars = NULL, estimate.method = "median", conf.int = T,conf.level = 0.95) 
        
        # SNP
        DF[DF$sample.id==n,"Beta1SNPcount"] <- output[1,"beta1"] <- df[df$term=="betaPropSNP","estimate"]$estimate
        DF[DF$sample.id==n,"LowCI1SNPcount"] <-output[1,"LowCI1"] <- df[df$term=="betaPropSNP","conf.low"]$conf.low
        DF[DF$sample.id==n,"HighCI1SNPcount"] <-output[1,"HighCI1"] <- df[df$term=="betaPropSNP","conf.high"]$conf.high
        DF[DF$sample.id==n,"CrossZero1SNPcount"] <-  output[1,"CrossZero1"] <- ifelse((output$LowCI1[1]>0 & output$HighCI1[1]<0) | (output$LowCI1[1]<0 & output$HighCI1[1]>0), "-" ,"*")
        
        # Height
        output[1,"beta2"] <- df[df$term=="betaHeight","estimate"]$estimate
        output[1,"LowCI2"] <- df[df$term=="betaHeight","conf.low"]$conf.low
        output[1,"HighCI2"] <- df[df$term=="betaHeight","conf.high"]$conf.high
        output[1,"CrossZero2"] <- ifelse((output$LowCI2[1]>0 & output$HighCI2[1]<0) | (output$LowCI2[1]<0 & output$HighCI2[1]>0), "-" ,"*")
        
        #saveRDS(output,file=paste0("outputs/ValidationSurvivalCommonGardens/CandidateSNPs/SNPCountanalysis_",i,"_",CovSet,"_",j,"_",site,".rds"))
        
      }
      
      saveRDS(DF,file=paste0("outputs/ValidationSurvivalCommonGardens/RandomSNPs/DF_",i,"_",CovSet,"_",j,"_",site,".rds"))
      
      }}}

parallel::stopCluster(cl)

end_time <- Sys.time()
end_time - start_time
