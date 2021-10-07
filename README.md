## U5MR Estimation Pipeline

This Github Repo documents the implementation of U5MR estimation pipeline described in 
the DHS spatial report (available at https://dhsprogram.com/publications/publication-SAR21-Spatial-Analysis-Reports.cfm)

Scripts are all written in R and are built upon the R package SUMMER.  Required data include DHS survey data (BR recode) and worldpop population density surface.

The pipeline mainly features the implementation of stratified subnational beta-binomial model, accompanied by procedure of acquiring the data and post-estimation visualization. We recommend the readers to check out the vignettes as they contain step-by-step instructions. We also encourage the readers to reference the original report for technical details.  
