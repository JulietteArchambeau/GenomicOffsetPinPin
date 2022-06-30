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

  - `prov`: ID of the population.
  
  - `longitude_prov`:  longitude of the geographic location of the population.
  
  - `latitude_prov`: latitude of the geographic location of the population.

  - `site`: common garden in which the phenotype (height and survival) was measured.
  
  - `clon`: ID of the clone (i.e. genotype).
  
  - `block`: ID of the block in which the phenotype was measured.
  
  - `tree`: ID of the tree on which the phenotype was measured.
  
  - `obs`: ID of the observation (i.e. measurement).
  
  - `age`: Tree age (in months) when the phenotype was measured.
  
  - `survival`: binary variable indicating wheter the tree was alive (1) or dead (0).
  
  - `height`: height measurement in cm. 
  
  - `altitude_prov`: altitude of the geographic location of the population.
  
  - `latitude_site`: latitude of the geographic location of the common garden in which the phenotype was measured.
  
  - `longitude_site`: longitude of the geographic location of the common garden in which the phenotype was measured.
  
  - `bio1WC_prov` to `bio9WC_prov`: [Bioclimatic variables](https://www.worldclim.org/data/bioclim.html) from WordlClim (30 seconds resolution, $\sim$ 1 km2 resolution; average for the years 1970-2000) extracted at the geographic location of the population;
  
    - `bio1WC_prov`: mean annual temperature (in C°).
    
    - `bio10WC_prov`: mean temperature of the warmest quarter (in C°).
    
    - `bio11WC_prov`: mean temperature of the coldest quarter (in C°).
    
    - `bio12WC_prov`: annual precipitation (in mm).
    
    - `bio13WC_prov`: precipitation of the wettest month (in mm).
    
    - `bio14WC_prov`: precipitation of the driest month (in mm).
    
    - `bio15WC_prov`:  precipitation seasonality (coefficient of variation).
    
    - `bio16WC_prov`: precipitation of the wettest quarter (in mm).
    
    - `bio17WC_prov`: precipitation of driest quarter (in mm).
    
    - `bio18WC_prov`: precipitation of the warmest quarter (in mm).
    
    - `bio19WC_prov`: precipitation of the coldest quarter (in mm).
    
    - `bio2WC_prov`: mean diurnal range (mean of monthly (maximum temperature - minimum temperature)).
    
    - `bio3WC_prov`: isothermality (bio2/bio7) (×100).
    
    - `bio4WC_prov`: temperature seasonality (standard deviation ×100).
    
    - `bio5WC_prov`:  maximum temperature of the warmest month (in C°).
    
    - `bio6WC_prov`: minimum temperature of the coldest month (in C°).
    
    - `bio7WC_prov`: temperature annual range (bio5-bio6)
    
    - `bio8WC_prov`: mean temperature of the wettest quarter (in C°).
    
    - `bio9WC_prov`: mean temperature of the driest quarter (in C°).
    
  - `clay_top_prov` to `depth_roots_prov`: soil-related variables from the [European Soil Database Derived data](https://esdac.jrc.ec.europa.eu/content/european-soil-database-derived-data) at 1-km resolution extracted at the geographic location of the population:
  
    - `clay_top_prov`: clay content in the topsoil (0-30 cm) in %.
    
    - `clay_top_prov`: silt content in the topsoil (0-30 cm) in %.
    
    - `clay_top_prov`: sand content in the topsoil (0-30 cm) in %.
    
    - `water_top_prov`: total available water content in mm.
    
    - `depth_roots_prov`: depth available to roots (in cm).
  
  - `Q1`: proportion of assignment to the northern African (NA) gene pool for each clone.
  
  - `Q2`: proportion of assignment to the Corsican (C) gene pool for each clone.
  
  - `Q3`: proportion of assignment to the central Spain (CS) gene pool for each clone.
  
  - `Q4`: proportion of assignment to the French Atlantic (FA) gene pool for each clone.
  
  - `Q5`: proportion of assignment to the Iberian Atlantic (IA) gene pool for each clone.
  
  - `Q6`: proportion of assignment to the south-eastern Spain (SES) gene pool for each clone.
  
  - `max.Q`: main gene pool for each clone (i.e. the gene pool constituting the higher proportion of population ancestry)
    
  - `TRI`: topographic ruggedness index (unitless).
  
  - `BurnedArea`: average of the monthly burned area from June 1995 to December 2014 (hectares)
  
  
### Raw genomic data

**In the dataset `RawGenomicData.csv`**

This file contains the genotype (noted with letters, e.g. A/A) of each clone.

14,016 SNPs, 529 clones.

Missing data are indicated with `---`.

Meaning of the columns:

  - `clone`: clone ID.
  
  - `assay`: Assay in which the clone was genotyped, either the Infinium assay (`only_Inf`), the Axiom assay (`only_Affx`) or both assays (`both_Inf_Affx`).
  
  - `snp_1` --> `snp_14016`:  genotype for each of the 14,016 SNPs.
  
  
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



### SNP codes

**In the dataset `SnpCodesMatching.csv`**

This file contains the correspondence among the SNP codes of the different assays (Axiom and Illumina Infinium), the codes used in the present study and the original SNP codes.

Meaning of the columns:

  - `original_ID`: original SNP ID
  
  - `affx_ID`: SNP ID from the Axiom assay.
  
  - `infinium_ID`: SNP ID from the Illumina Infinium assay.
  
  - `snp_ID`: SNP ID used in the present study.
  
  