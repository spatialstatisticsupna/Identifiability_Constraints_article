# Identifiability constraints in spatio-temporal areal models
This repository contains the R code to fit with INLA the spatio-temporal models described in _"In spatio-temporal disease mapping models, identifiability constraints affect PQL and INLA results"_ [(Goicoa et al., 2018)](https://doi.org/10.1007/s00477-017-1405-0).


## Table of contents

- [Data](#Data)
- [R code](#R-code)
- [References](#References)


# Data
Female breast cancer mortality data (ICD-10 code 50) in Spanish provinces during the period 1990-2010.

- [**BreastCancer_ESP.Rdata**](https://github.com/spatialstatisticsupna/Identifiability_Constraints_article/blob/master/data/BreastCancer_ESP.Rdata)
  
  This .Rdata contains the following objects
  - **_Data_**: `data.frame` object with the number of observed and expected cases (_'Counts'_ and _'Expected'_ variables, respectively) for each province (_'Area'_) and time period (_'Year'_) for female breast cancer mortality data.
  - **_Carto_ESP_**: `sf` object containing the spatial polygons of the Spanish provinces. The data contains a `data.frame` with 50 rows and  _'Area'_ (character vector of geographic identifiers), _'Name'_ (character vector of province names), _'Longitude'_ (numeric vector of longitude values), _'Latitude'_ (numeric vector of latitude values) and _'geometry'_ (sfc_MULTIPOLYGON) variables.
	
  
- [**Esp_prov_nb.graph**](https://github.com/spatialstatisticsupna/Identifiability_Constraints_article/blob/master/data/Esp_prov_nb.graph)
  
  An inla.graph object with the spatial neighbourhood structure of the 50 provinces of Spain.
  

# R code
R code to fit with INLA (http://www.r-inla.org/) the spatio-temporal CAR models described in Goicoa et al. (2018). All the R files are written by the authors of the paper.

- [**CARmodels_INLA.R**](https://github.com/spatialstatisticsupna/Identifiability_Constraints_article/blob/master/R/CARmodels_INLA.R)

  Main script including the required functions to fit in INLA the spatio-temporal CAR models with different types of interactions. Slight modifications of the original code described in the Appendix section of Goicoa et al. (2018) have been introduced in order to be compatible with the 20.12.10 testing version of INLA. More precisely, we eliminate the redundant constraints in the `extraconstr` argument of the `INLA::inla()` function. 
  
- [**posterior_lincombs.R**](https://github.com/spatialstatisticsupna/Identifiability_Constraints_article/blob/master/R/posterior_lincombs.R)

  Auxiliary script including the code to compute the posterior marginal distributions of the spatial, temporal and spatio-temporal patterns (Adin et al., 2017) defined as linear combinations of the log-risks using the `INLA::inla.make.lincombs()` function.


- [**Results.R**](https://github.com/spatialstatisticsupna/Identifiability_Constraints_article/blob/master/R/Results.R)

  Script to analyse the results of previously fitted INLA models.
  

# References
[Adin, A., Martínez-Beneito, M.A., Botella-Rocamora, P., Goicoa, T., and Ugarte, M.D. (2017). Smoothing and high risk areas detection in space-time disease mapping: a comparison of P-splines, autoregressive, and moving average models. _Stochastic Environmental Research and Risk Assessment_, __31(2)__, 403-415.](https://doi.org/10.1007/s00477-016-1269-8)

[Goicoa, T., Adin, A., Ugarte, M.D., and Hodges, J.S. (2018). In spatio-temporal disease mapping models, identifiability constraints affect PQL and INLA results. _Stochastic Environmental Research and Risk Assessment_, __32(3)__, 749-770.](https://doi.org/10.1007/s00477-017-1405-0)
