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


### References

The ensemble model is from:
* Valavi, R., Guillera‐Arroita, G, Lahoz‐Monfort, J.J. & Elith, J. (2021) Predictive performance of presence-only species distribution models: a benchmark study with reproducible code. Ecological Monographs. [DOI:10.1002/ecm.1486](DOI:10.1002/ecm.1486)

* 