# Data and code for the paper: 'XXXXX'

**Authors:** Juliette Archambeau<sup>1</sup>, Marta Benito Garzón<sup>1</sup>, Marina de Miguel<sup>1,2</sup>, Alexandre Changenet<sup>1</sup>, Camilla Avanzi<sup>3</sup>, Francesca Bagnoli<sup>3</sup>, Frédéric Barraquand<sup>4</sup>, Giovanni G. Vendramin<sup>3</sup> and Santiago C. González-Martínez<sup>1</sup>

**1** INRAE, Univ. Bordeaux, BIOGECO, F-33610 Cestas, France

**2** EGFV, Univ. Bordeaux, Bordeaux Sciences Agro, INRAE, ISVV, F-33882, Villenave d’Ornon, France

**3** Institute of Biosciences and BioResources, National Research Council, 50019 Sesto Fiorentino, Italy

**4** CNRS, Institute of Mathematics of Bordeaux, F-33400 Talence, France

**Corresponding author:** Juliette Archambeau, juli.archambeau@orange.fr


## Data in the DRYAD repository

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
  
  