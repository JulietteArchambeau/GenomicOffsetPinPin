# Data and code for the paper: 'XXXXX'

**Authors:** Juliette Archambeau<sup>1</sup>, Marta Benito Garzón<sup>1</sup>, Marina de Miguel<sup>1,2</sup>, Alexandre Changenet<sup>1</sup>, Camilla Avanzi<sup>3</sup>, Francesca Bagnoli<sup>3</sup>, Frédéric Barraquand<sup>4</sup>, Giovanni G. Vendramin<sup>3</sup> and Santiago C. González-Martínez<sup>1</sup>

**1** INRAE, Univ. Bordeaux, BIOGECO, F-33610 Cestas, France

**2** EGFV, Univ. Bordeaux, Bordeaux Sciences Agro, INRAE, ISVV, F-33882, Villenave d’Ornon, France

**3** Institute of Biosciences and BioResources, National Research Council, 50019 Sesto Fiorentino, Italy

**4** CNRS, Institute of Mathematics of Bordeaux, F-33400 Talence, France

**Corresponding author:** Juliette Archambeau, juli.archambeau@orange.fr


## Data in the DRYAD repository


### Phenotypic, environmental and population structure data

**In the dataset `PhenoEnvPopStructureDataset.csv`**

This file contains information related to the geographic location and the environment of origin of the populations, their proportition of assignement to each gene pool (i.e. the neutral population genetic structure) and the phenotypes (i.e. morality and height) in the different common gardens (used in the validation analysis).

Meaning of the columns:

* `prov`: ID of the population.
* `longitude_prov`:  longitude of the geographic location of the population.
* `latitude_prov`: latitude of the geographic location of the population.
* `site`: common garden in which the phenotype (height and survival) was measured.
* `clon`: ID of the clone (i.e. genotype).
* `block`: ID of the block in which the phenotype was measured.
* `tree`: ID of the tree on which the phenotype was measured.
* `obs`: ID of the observation (i.e. measurement).
* `age`: Tree age (in months) when the phenotype was measured.
* `survival`: binary variable indicating wheter the tree was alive (1) or dead (0).
* `height`: height measurement in cm. 
* `altitude_prov`: altitude of the geographic location of the population.
* `latitude_site`: latitude of the geographic location of the common garden in which the phenotype was measured.
* `longitude_site`: longitude of the geographic location of the common garden in which the phenotype was measured.

* `bio1WC_prov` to `bio9WC_prov`: [Bioclimatic variables](https://www.worldclim.org/data/bioclim.html) from WordlClim (30 seconds resolution, $\sim$ 1 km2 resolution; average for the years 1970-2000) extracted at the geographic location of the population;
    * `bio1WC_prov`: mean annual temperature (in C°).
    * `bio10WC_prov`: mean temperature of the warmest quarter (in C°).
    * `bio11WC_prov`: mean temperature of the coldest quarter (in C°).
    * `bio12WC_prov`: annual precipitation (in mm).
    * `bio13WC_prov`: precipitation of the wettest month (in mm).
    * `bio14WC_prov`: precipitation of the driest month (in mm).
    * `bio15WC_prov`:  precipitation seasonality (coefficient of variation).
    * `bio16WC_prov`: precipitation of the wettest quarter (in mm).
    * `bio17WC_prov`: precipitation of driest quarter (in mm).
    * `bio18WC_prov`: precipitation of the warmest quarter (in mm).
    * `bio19WC_prov`: precipitation of the coldest quarter (in mm).
    * `bio2WC_prov`: mean diurnal range (mean of monthly (maximum temperature - minimum temperature)).
    * `bio3WC_prov`: isothermality (bio2/bio7) (×100).
    * `bio4WC_prov`: temperature seasonality (standard deviation ×100).
    * `bio5WC_prov`:  maximum temperature of the warmest month (in C°).
    * `bio6WC_prov`: minimum temperature of the coldest month (in C°).
    * `bio7WC_prov`: temperature annual range (bio5-bio6)
    * `bio8WC_prov`: mean temperature of the wettest quarter (in C°).
    * `bio9WC_prov`: mean temperature of the driest quarter (in C°).
  
* `clay_top_prov` to `depth_roots_prov`: soil-related variables from the [European Soil Database Derived data](https://esdac.jrc.ec.europa.eu/content/european-soil-database-derived-data) at 1-km resolution extracted at the geographic location of the population:
    * `clay_top_prov`: clay content in the topsoil (0-30 cm) in %.
    * `clay_top_prov`: silt content in the topsoil (0-30 cm) in %.
    * `clay_top_prov`: sand content in the topsoil (0-30 cm) in %.
    * `water_top_prov`: total available water content in mm.
    * `depth_roots_prov`: depth available to roots (in cm).
    
* `Q1`: proportion of assignment to the northern African (NA) gene pool for each clone.
* `Q2`: proportion of assignment to the Corsican (C) gene pool for each clone.
* `Q3`: proportion of assignment to the central Spain (CS) gene pool for each clone.
* `Q4`: proportion of assignment to the French Atlantic (FA) gene pool for each clone.
* `Q5`: proportion of assignment to the Iberian Atlantic (IA) gene pool for each clone.
* `Q6`: proportion of assignment to the south-eastern Spain (SES) gene pool for each clone.
* `max.Q`: main gene pool for each clone (i.e. the gene pool constituting the higher proportion of population ancestry)

* `TRI`: topographic ruggedness index (unitless).
* `BurnedArea`: average of the monthly burned area from June 1995 to December 2014 (hectares)
  
  
### Raw genomic data

**In the dataset `RawGenomicData.csv`**

This file contains the genotype (noted with letters, e.g. A/A) of each clone.

14,016 SNPs, 529 clones.

Missing data are indicated with `---`.

Meaning of the columns:

* `clone`: clone ID.
* `assay`: Assay in which the clone was genotyped, either the Infinium assay (`only_Inf`), the Axiom assay (`only_Affx`) or both assays (`both_Inf_Affx`).
* `snp_1` --> `snp_14016`:  genotype for each of the 14,016 SNPs.

### SNP codes

**In the dataset `SnpCodesMatching.csv`**

This file contains the correspondence among the SNP codes of the different assays (Axiom and Illumina Infinium), the codes used in the present study and the original SNP codes.

Meaning of the columns:

* `original_ID`: original SNP ID
* `affx_ID`: SNP ID from the Axiom assay.
* `infinium_ID`: SNP ID from the Illumina Infinium assay.
* `snp_ID`: SNP ID used in the present study.
  
  
  
### Formatted and filtered genomic data

**In the dataset `FormattedFilteredGenomicData_454clones_9817snps.csv`**

This file contains the genotype (noted with numbers, i.e. 0, 1 or 2) of each clone after formatting and filtering.

Genomic data were filtered for clones with more than 18% missing data, MAF lower than 1% and SNPs with missing data higher than 20% (more details in the script `1_FormattingGenomicData.Rmd`).

9,817 SNPs (in rows), 454 clones (in columns).

Missing data are indicated with `NA`. 


### Imputed genomic data

**In the dataset `ImputedGenomicData_454clones_9817snps.csv`**

Same dataset as `FormattedFilteredGenomicData_454clones_9817snps.csv` but with imputed genomic data. The missing SNP values were imputed based on the most common allele within the main gene pool of the clone (more details in the script `1_FormattingGenomicData.Rmd`).

9,817 SNPs (in rows), 454 clones (in columns).

### Data for BayPass analysis

**In the folder `BayPassGEAAnalysis`**

Associated script: `3_GEAanalyses_BayPass.Rmd`.

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

**In the folder `RDAGEAAnalysis`**

Associated script: `4_GEAanalyses_RDA.Rmd`.

The **list of candidate SNPs** identified in the RDA analysis is provided in the file `CandSNPsRDA.csv`. Here is the meaning of the columns:

* `snp`: SNP ID used in the present study. 
* `bio5` to `BurnedArea`: correlation coefficients between each candidate SNP and each environmental covariate.
* `predictor`: environmental covariate with the highest association with the candidate SNP.
* `correlation`: correlation coefficient between the candidate SNP and the environmental covariate with the highest association with the candidate SNP (i.e. the `predictor` of the previsou column).


### LD of the three SNP sets

**In the folder `LD`**

Associated script: `5_CreatingFormattingCalculatingLDSubsetsCandidateSNPs.Rmd`.

Outputs of the `LD` function from the *genetics* R package. This function calculates **pairwise linkage disequilibrium between genetic markers** and returns a list of 5 elements, among hich is the LD estimate. 

LD was calculate for the three sets of SNPs: the reference SNPs (file `outputsLDGenetics_Ref.rds`), the merged candidate SNPs (file `outputsLDGenetics_Mer.rds`) and the common candidate SNPs (file `outputsLDGenetics_Com.rds`).

### Genomic input files for the GF and GDM analyses

Associated script: `5_CreatingFormattingCalculatingLDSubsetsCandidateSNPs.Rmd`.

The GF analysis requires **population allele frequencies** as input files. The population allele frequencies of each set of SNPs can be found in the file `ListAlleleFrequencies.rds` (in a list).

The GDM analysis requires **pairwise $F_{ST}$ matrices** as input files. These matrices were calculated for each set of SNPs with the function `pairwise.WCfst` of the *hierfstat* R package and were then scaled, so that the $F_{ST}$ values are between 0 and 1. The pairwise  $F_{ST}$ matrices **before scaling** are stored in a list in the file `ListPairwiseFst.rds`. The pairwise  $F_{ST}$ matrices **after scaling** are stored in the file `ListPairwiseFstScaled.rds`.
