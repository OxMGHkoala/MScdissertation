# EDAM Cost-Effectiveness Model
This repository contains the R script and associated data for the cost-effectiveness analysis of the Electronic Clinical Decision Support for Acute Fever Management (EDAM) application.
S5_EDAM_all_Datasets_21May2025.dta is the clinical data and the only Supplementary Information file which must be imported into R to run the code (see directions below). 
All other files in this repository are complementary and for reading purposes, and do not need to be imported into R.

## How to run
1. Clone this repository.
2. Set working directory in R to the cloned folder.
3. Run S2_1090618_EDAMCEA_MGHDissertation.R script
4. uncomment this line: clinicaldata <- read.dta()
5. In the above line, insert your new file path to S5_EDAM_all_Datasets_21May2025.dta

## Supplementary Information Materials 
S1: Appendix to Main Report 
S2: Model (R Script)*
S3: Presentation of cost inputs used in the model
S4: Presentation of disability weights used in the model
S5: EDAM clinical data*
S6: 'Data Dictionary' including explanations of clinical data variables 

## Authors
YK, CP, CC, EW
