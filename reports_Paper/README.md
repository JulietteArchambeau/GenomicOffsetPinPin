# Data and code for the paper: 'XXXXX'

**Authors:** Juliette Archambeau<sup>1</sup>, Marta Benito Garzón<sup>1</sup>, Marina de Miguel<sup>1,2</sup>, Alexandre Changenet<sup>1</sup>, Camilla Avanzi<sup>3</sup>, Francesca Bagnoli<sup>3</sup>, Frédéric Barraquand<sup>4</sup>, Giovanni G. Vendramin<sup>3</sup> and Santiago C. González-Martínez<sup>1</sup>

**1** INRAE, Univ. Bordeaux, BIOGECO, F-33610 Cestas, France

**2** EGFV, Univ. Bordeaux, Bordeaux Sciences Agro, INRAE, ISVV, F-33882, Villenave d’Ornon, France

**3** Institute of Biosciences and BioResources, National Research Council, 50019 Sesto Fiorentino, Italy

**4** CNRS, Institute of Mathematics of Bordeaux, F-33400 Talence, France

**Corresponding author:** Juliette Archambeau, juli.archambeau@orange.fr


## Data in the DRYAD repository


### Phenotypic, environmental and population structure data

**Dataset `PopulationPopStructureEnvData.csv`**

This file contains information related to each population and clone, i.e. their geographical location (population-level), their environment of origin (population-level) and their proportition of assignement to each gene pool (clone-level), that is the neutral population genetic structure).

Meaning of the columns:

* `prov`: ID of the population.
* `longitude`:  longitude of the geographic location of the population.
* `latitude`: latitude of the geographic location of the population.
* `clon`: ID of the clone (i.e. genotype).
* `altitude`: altitude of the geographic location of the population.
* `bio1` to `bio9`: [Bioclimatic variables](https://www.worldclim.org/data/bioclim.html) from WordlClim (30 seconds resolution, $\sim$ 1 km2 resolution; average for the years 1970-2000) extracted at the geographic location of the population;
    * `bio1`: annual temperature (in C°).
    * `bio10`: mean temperature of the warmest quarter (in C°).
    * `bio11`: mean temperature of the coldest quarter (in C°).
    * `bio12`: annual precipitation (in mm).
    * `bio13`: precipitation of the wettest month (in mm).
    * `bio14`: precipitation of the driest month (in mm).
    * `bio15`:  precipitation seasonality (coefficient of variation).
    * `bio16`: precipitation of the wettest quarter (in mm).
    * `bio17`: precipitation of driest quarter (in mm).
    * `bio18`: precipitation of the warmest quarter (in mm).
    * `bio19`: precipitation of the coldest quarter (in mm).
    * `bio2`: mean diurnal range (mean of monthly (maximum temperature - minimum temperature)).
    * `bio3`: isothermality (bio2/bio7) (×100).
    * `bio4`: temperature seasonality (standard deviation ×100).
    * `bio5`: maximum temperature of the warmest month (in C°).
    * `bio6`: minimum temperature of the coldest month (in C°).
    * `bio7`: temperature annual range (bio5-bio6)
    * `bio8`: mean temperature of the wettest quarter (in C°).
    * `bio9`: mean temperature of the driest quarter (in C°).
  
* `clay_top` to `depth_roots`: soil-related variables from the [European Soil Database Derived data](https://esdac.jrc.ec.europa.eu/content/european-soil-database-derived-data) at 1-km resolution extracted at the geographic location of the population:
    * `clay_top`: clay content in the topsoil (0-30 cm) in %.
    * `silt_top`: silt content in the topsoil (0-30 cm) in %.
    * `sand_top`: sand content in the topsoil (0-30 cm) in %.
    * `water_top`: total available water content (in mm).
    * `depth_roots`: depth available to roots (in cm).
    
* `Q1`: proportion of assignment to the northern African (NA) gene pool for each clone.
* `Q2`: proportion of assignment to the Corsican (C) gene pool for each clone.
* `Q3`: proportion of assignment to the central Spain (CS) gene pool for each clone.
* `Q4`: proportion of assignment to the French Atlantic (FA) gene pool for each clone.
* `Q5`: proportion of assignment to the Iberian Atlantic (IA) gene pool for each clone.
* `Q6`: proportion of assignment to the south-eastern Spain (SES) gene pool for each clone.
* `max.Q`: main gene pool for each clone (i.e. the gene pool constituting the higher proportion of population ancestry)

* `TRI`: topographic ruggedness index (in m) calculated using SAGA v2.3.1 (Conrad *et al*. 2015) based on topographic data generated from NASA’s Shuttle Radar Topography Mission (SRTM) at 90 m resolution. TRI quantifies the terrain heterogeneity, i.e. differences in elevation between adjacent cells (Riley *et al*. 1999)
* `BurnedArea`: average of the monthly burned area from June 1995 to December 2014 (in hectares) extracted from the GFED4 database at 0.25 degrees resolution (∼28 km resolution) (Giglio
*et al*. 2013).


### Raw genomic data

**Dataset `RawGenomicData.csv`**

This file contains the genotype (noted with letters, e.g. A/A) of each clone.

14,016 SNPs, 529 clones.

Missing data are indicated with `---`.

Meaning of the columns:

* `clone`: clone ID.
* `assay`: Assay in which the clone was genotyped, either the Infinium assay (`only_Inf`), the Axiom assay (`only_Affx`) or both assays (`both_Inf_Affx`).
* `snp_1` --> `snp_14016`:  genotype for each of the 14,016 SNPs.

### SNP codes

**Dataset `SnpCodesMatching.csv`**

This file contains the correspondence among the SNP codes of the different assays (Axiom and Illumina Infinium), the codes used in the present study and the original SNP codes.

Meaning of the columns:

* `original_ID`: original SNP ID
* `affx_ID`: SNP ID from the Axiom assay.
* `infinium_ID`: SNP ID from the Illumina Infinium assay.
* `snp_ID`: SNP ID used in the present study.
  
  
  
### Formatted and filtered genomic data

**Dataset `FormattedFilteredGenomicData_454clones_9817snps.csv`**

This file contains the genotype (noted with numbers, i.e. 0, 1 or 2) of each clone after formatting and filtering.

Genomic data were filtered for clones with more than 18% missing data, MAF lower than 1% and SNPs with missing data higher than 20% (more details in the script `1_FormattingGenomicData.Rmd`).

9,817 SNPs (in rows), 454 clones (in columns).

Missing data are indicated with `NA`. 


### Imputed genomic data

**Dataset `ImputedGenomicData_454clones_9817snps.csv`**

Same dataset as `FormattedFilteredGenomicData_454clones_9817snps.csv` but with imputed genomic data. The missing SNP values were imputed based on the most common allele within the main gene pool of the clone (more details in the script `1_FormattingGenomicData.Rmd`).

9,817 SNPs (in rows), 454 clones (in columns).

### Data for BayPass analysis

**Folder `BayPassGEAAnalysis` / Associated script `3_GEAanalyses_BayPass.Rmd`**

Input files required for the BayPass analysis are **allele count data**. The file `AlleleCounts_9817snps454clones.csv` contains the counts of the minor allele in each population.

In BayPass, the input genotyping data file has to be organized as a matrix with *nsnp* rows and *2 ∗ npop* columns (the first two columns correspond to one population, the next two columns to the second population, etc.), with space as row field separator. More details in the [BayPass manual](http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.2.pdf). The file `PreFileBayPass_9817snps454clones` contains the allele count data in the required format. 

The input files with the values of the **environmental covariates** for each population are in the subfolder `EnvironmentalVariables`.

We first estimated the popupation covariance matrix and the output files are located in the folder `CovarianceMatrixOmega`.

For each environmental covariate, we performed five independant runs of the **standard covariate model** using the **Important Sampling (IS) algorithm**. The outfiles of all independant runs can be found in the covariate-specific folders: `ISruns_bio5`, `ISruns_bio6`, `ISruns_bio12`, `ISruns_bio15`, `ISruns_water_top`, `ISruns_depth_roots`, `ISruns_BurnedArea` and `ISruns_TRI`. 

The **list of candidate SNPs** identified in the BayPass analysis is provided in the file `CandSNPsBayPassIS.csv`. Here is the meaning of the columns:

* `snp`: SNP ID used in the present study. 
* `medianBF`: **median estimate of the Bayes Factor** in dB units across the five independant runs (measuring the support of the association of each SNP with each environmental covariate).
* `medianBeta`: **median estimate of the regression coefficients** ($\beta_i$ in the standard covariate model) across the five independant runs (measuring the strength of the association of each SNP with each environmental covariate).
* `medianEBP`: **median estimate of the empirical Bayesian P-values** in the $\log_{10}$ scale (measuring the support in favor of a non-null regression coefficient).
* `COVARIABLE`: **environmental covariate** associated with the candidate SNP.

### Data for RDA analysis

**Folder `RDAGEAAnalysis` / Associated script: `4_GEAanalyses_RDA.Rmd`**

The **list of candidate SNPs** identified in the RDA analysis is provided in the file `CandSNPsRDA.csv`. Here is the meaning of the columns:

* `snp`: SNP ID used in the present study. 
* `bio5` to `BurnedArea`: correlation coefficients between each candidate SNP and each environmental covariate.
* `predictor`: environmental covariate with the highest association with the candidate SNP.
* `correlation`: correlation coefficient between the candidate SNP and the environmental covariate with the highest association with the candidate SNP (i.e. the `predictor` of the previsou column).


### LD of the three SNP sets

**Folder `LD` / Associated script: `5_CreatingFormattingCalculatingLDSubsetsCandidateSNPs.Rmd`**

Outputs of the `LD` function from the *genetics* R package. This function calculates **pairwise linkage disequilibrium between genetic markers** and returns a list of 5 elements, among hich is the LD estimate. 

LD was calculate for the three sets of SNPs: the reference SNPs (file `outputsLDGenetics_Ref.rds`), the merged candidate SNPs (file `outputsLDGenetics_Mer.rds`) and the common candidate SNPs (file `outputsLDGenetics_Com.rds`).

### Genomic input files for the GF and GDM analyses

**Associated script `5_CreatingFormattingCalculatingLDSubsetsCandidateSNPs.Rmd`**

The GF analysis requires **population allele frequencies** as input files. The population allele frequencies of each set of SNPs can be found in the file `ListAlleleFrequencies.rds` (in a list).

The GDM analysis requires **pairwise $F_{ST}$ matrices** as input files. These matrices were calculated for each set of SNPs with the function `pairwise.WCfst` of the *hierfstat* R package and were then scaled, so that the $F_{ST}$ values are between 0 and 1. The pairwise  $F_{ST}$ matrices **before scaling** are stored in a list in the file `ListPairwiseFst.rds`. The pairwise  $F_{ST}$ matrices **after scaling** are stored in the file `ListPairwiseFstScaled.rds`.

### Shapefile of the species distribution

**Folder `Mapping`**

Shapefile corresponding to the maritime pine distribution used in the maps. This distribution was build based on the EUFORGEN distribution (http://www.euforgen.org/) and 10-km radius areas around the National Forest Inventory plots with maritime pines. However, this remains a rough approximation of the actual distribution of maritime pine.

### Survival data in common gardens

**Dataset `SurvivalDataCommonGarden.csv`**

Survival data in the common gardens of Cáceres and Madrid (from the CLONAPIN network) in which a severe summer drought exacerbated by clay soils killed 92% and 72% of the trees.

Meaning of the columns:

* `site`: common garden ID: Cáceres or Madrid.
* `block`: block ID.
* `prov`: population ID.
* `clon`: clone (i.e. genotype) ID.
* `tree`: tree ID. 
* `survival`: binary variable indicating whether the tree was recorded as dead (0) or alive (1).


### Height data in common gardens

**Dataset `HeightDataCommonGarden.csv`**

Height data from the five common gardens of the CLONAPIN network, respectively located in Bordeaux, Portugal, Madrid, Asturias and Cáceres. 

Meaning of the columns:

* `site`: common garden ID.
* `block`: block ID.
* `prov`: population ID.
* `clon`: clone (i.e. genotype) ID.
* `tree`: tree ID. 
* `MAD_htdec11`: height measurements in Madrid in December 2011 when the trees were 13-month old.
* `CAC_htdec11`: height measurements in Cáceres in December 2011 when the trees were 8-month old.
* `BDX_htnov13`: height measurements in Bordeaux in November 2013 when the trees were 25-month old.
* `BDX_htnov18`: height measurements in Bordeaux in November 2018 when the trees were 85-month old.
* `POR_htjan12`: height measurements in Portugal in January 2012 when the trees were 11-month old.
* `POR_htmay13`: height measurements in Portugal in May 2013 when the trees were 27-month old.
* `AST_htdec11`: height measurements in Asturias in December 2011 when trees were 10-month old.
* `AST_htmar14`: height measurements in Asturias in March 2014 when the tree were 37-month old.


### Environmental data in common gardens

**Dataset `EnvDataCommonGarden.csv`**

Values of the environmental variables used in the present study at the location of each common garden. Details about the extraction of each variable is the description of the dataset `PopulationPopStructureEnvData.csv` (see above).

Meaning of the columns:

* `bio1`: annual temperature (in C°).
* `bio5`: maximum temperature of the warmest month (in C°).
* `bio6`: minimum temperature of the coldest month (in C°).
* `bio12`: annual precipitation (in mm).  
* `bio15`: precipitation seasonality (coefficient of variation).
* `water_top`: total available water content (in mm).
* `depth_roots`: depth available to roots (in cm).
* `TRI`: topographic ruggedness index (unitless).
* `BurnedArea`: average of the monthly burned area from June 1995 to December 2014 (hectares).

### Mortality data from National Forest Inventories

**Dataset `NFIdata.csv`**


## Scripts in Zenodo

The code included in the present *zenodo* repository was run on *R version 3.6.3* and *RStudio version 1.1.463* and constitutes the code necessary to replicate the analyses of the present study.

Here is what the different scripts are for: