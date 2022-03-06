################################################################
#########   load libraries
################################################################
rm(list = ls())
options(warn = -1)
#### Libraries ####
library(SUMMER)
library(classInt)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(rgdal)
library(scales)
library(INLA)
library(survey)
library(ggplot2)
library(raster)
library(maptools)
library(gridExtra)
library(mgcv)
library(caret)
library(geosphere)
library(rgeos)
library(haven)
library(labelled)
library(data.table)
options(gsubfn.engine = "R")
library(sqldf)
library(sp)
library(gstat)

# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home_dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")
countries <- countries <- scan(paste0(home_dir, "/countries_implemented.txt"), character(), quote = "")
country <- countries[length(countries)] # retrieve the country being analyzed
info.name <- paste0(country, "_general_info.Rdata")
load(file = paste0(home_dir,'/Info/',country,"/", info.name, sep='')) # load the country info

################################################################
#########   set directories
################################################################

data_dir <- paste0(home_dir,'/Data/', country) # set the directory to store the data
res_dir <- paste0(home_dir,'/Results/', country) # set the directory to store the results (e.g. fitted R objects, figures, tables in .csv etc.)

################################################################
#########   load files
################################################################

setwd(data_dir)

poly.path <- paste0("shapeFiles_gadm") # specify the folder of the country shape files
poly.layer.adm0 <- paste('gadm36', gadm.abbrev,
                         '0', sep = "_") # specify the name of the national shape file
poly.layer.adm1 <- paste('gadm36', gadm.abbrev,
                         '1', sep = "_") # specify the name of the admin1 shape file
poly.layer.adm2 <- paste('gadm36', gadm.abbrev,
                         '2', sep = "_") # specify the name of the admin2 shape file

poly.adm0 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm0)) # load the national shape file
poly.adm1 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm1)) # load the shape file of admin-1 regions
poly.adm2 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm2)) # load the shape file of admin-2 regions

proj4string(poly.adm0) <- proj4string(poly.adm1) <- proj4string(poly.adm2)
load(paste0('shapeFiles_gadm/', country, '_Amat.rda'))  # load the adjacency matrix
load(paste0('shapeFiles_gadm/', country, '_Amat_Names.rda'))  # load names of admin1 and admin2 regions


################################################################
#########   Final preprocessing
################################################################  

# load the DHS survey data and adjacency matrix
load(paste0(country,'_cluster_dat.rda'),
     envir = .GlobalEnv)
load( paste0(poly.path,'/', country, '_Amat.rda'))
load( paste0(poly.path, '/', country, '_Amat_Names.rda'))


mod.dat<-mod.dat[mod.dat$survey==as.character(survey_year),] # filter the data of the recent survey
mod.dat$years <- as.numeric(as.character(mod.dat$years))
mod.dat<-mod.dat[as.numeric(mod.dat$years)>=beg.year,]

mod.dat$strata <- mod.dat$urban
mod.dat$country <- as.character(country)


################################################################
#########  National stratified BB8 
################################################################

setwd(paste0(res_dir))

# fit the national BB8 stratified model
fit.natl.strat <- smoothCluster(data = mod.dat, family = "betabinomial",
                                Amat = NULL,  strata.time.effect = TRUE,
                                year_label = beg.year:end.year,
                                time.model = "rw2",
                                overdisp.mean = -7.5,
                                overdisp.prec = 0.39,
                                age.groups = levels(mod.dat$age),
                                age.n = c(1,11,12,12,12,12),
                                age.rw.group = c(1,2,3,3,3,3),
                                age.strata.fixed.group = c(1,2,3,3,3,3),
                                survey.effect = FALSE)

saveRDS(fit.natl.strat, paste0('Betabinomial/',
                               country,'_fit.strat.natl.3UR.rds'))

# load the national urban population fraction, this is needed to weight urban/rural estimators
natl.urb.weights <- readRDS(paste0('UR/U5_fraction/','natl_urban_weights.rds'))
natl.urb.weights$rural <- 1- natl.urb.weights$urban

# sample for national U5MR estimators and compute the estimates
res.natl.strat <- getSmoothed(inla_mod = fit.natl.strat, 
                              year_range = beg.year:end.year, 
                              year_label = beg.year:end.year, nsim = 1000, 
                              weight.strata = natl.urb.weights, 
                              weight.frame = NULL,
                              draws = NULL, save.draws = TRUE)

saveRDS(res.natl.strat, paste0('Betabinomial/',
                               country,'_res.strat.natl.3UR.rds'))

################################################################
#########   Admin-1 stratified model
################################################################

mod.dat$region <- mod.dat$admin1.char

# fit the admin1 BB8 stratified model
fit.strat.admin1 <- smoothCluster(data = mod.dat, family = "betabinomial",
                                  Amat = admin1.mat, strata.time.effect = TRUE,
                                  year_label = beg.year:end.year,
                                  time.model = "rw2", st.time.model = "ar1",
                                  pc.st.slope.u = 1, pc.st.slope.alpha = 0.01,
                                  age.strata.fixed.group = c(1,2,3,3,3,3),
                                  type.st = type.st,
                                  bias.adj.by = adj.varnames,
                                  survey.effect = FALSE)

saveRDS(fit.strat.admin1,paste0('Betabinomial/',
                                country,'_fit.strat.admin1.3UR.rds'))

# load the admin1 urban population fraction, this is needed to weight urban/rural estimators
weight.strata.adm1 <- readRDS(paste0('UR/U5_fraction/','admin1_urban_weights.rds'))

# sample for national U5MR estimators and compute the estimates
res.strat.admin1 <- getSmoothed(inla_mod = fit.strat.admin1, 
                                year_range = beg.year:end.year, 
                                year_label = beg.year:end.year, nsim = 1000, 
                                weight.strata = weight.strata.adm1, 
                                weight.frame = NULL,
                                draws = NULL, save.draws = TRUE)

saveRDS(res.strat.admin1,paste0('Betabinomial/',
                                country,'_res.strat.admin1.3UR.rds'))




################################################################
#########   Admin-2 stratified model
################################################################

mod.dat$region <- mod.dat$admin2.char

# fit the admin2 BB8 stratified model
fit.strat.admin2 <- smoothCluster(data = mod.dat, family = "betabinomial",
                                  Amat = admin2.mat, strata.time.effect = TRUE,
                                  year_label = beg.year:end.year,
                                  time.model = "rw2", st.time.model = "ar1",
                                  pc.st.slope.u = 1, pc.st.slope.alpha = 0.01,
                                  type.st = type.st,
                                  age.strata.fixed.group = c(1,2,3,3,3,3),
                                  survey.effect = FALSE)

saveRDS(fit.strat.admin2,paste0('Betabinomial/',
                                country,'_fit.strat.admin2.3UR.rds'))

# load the admin2 urban population fraction, this is needed to weight urban/rural estimators
weight.strata.adm2 <- readRDS(paste0('UR/U5_fraction/','admin2_urban_weights.rds'))

# sample for admin2 U5MR estimators and compute the estimates
res.strat.admin2 <- getSmoothed(inla_mod = fit.strat.admin2, 
                                year_range = beg.year:end.year, 
                                year_label = beg.year:end.year, nsim = 1000, 
                                weight.strata = weight.strata.adm2, 
                                weight.frame = NULL,
                                draws = NULL, save.draws = TRUE)

saveRDS(res.strat.admin2,paste0('Betabinomial/',
                                country,'_res.strat.admin2.3UR.rds'))
options(warn = 0)
