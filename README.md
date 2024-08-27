##  Ensemble Species Distribution Modelling: Assessing Climate Change Impacts

This repository provides a comprehensive guide to ensemble species distribution modeling (SDM) under climate change

### Species of interest




### Cimate and environmetal covaraites



### How to run scripts

Order of running scripts:

1. Run `climate_data.R` to get download climate variables from CHELSA
2. Run `climate_pca.R` to compute the PCA of bioclimatic variables
3. Run `species_data.R` to download the species data
4. Run either `modelling_plant.R` or `modelling_cat.R` to predict the two speecies.
