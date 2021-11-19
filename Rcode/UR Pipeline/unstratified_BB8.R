################################################################
#########   load libraries
################################################################
rm(list = ls())

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

code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home_dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")
countries <- countries <- scan(paste0(home_dir, "/countries_implemented.txt"), character(), quote = "")
country <- countries[length(countries)]
info.name <- paste0(country, "_general_info.Rdata")
load(file = paste0(home_dir,'/Info/',country,"/", info.name, sep=''))

################################################################
#########   set directories
################################################################

data_dir <- paste0(home_dir,'/Data/')
res_dir <- paste0(home_dir,'/Results/')

################################################################
#########   load files
################################################################

setwd(paste(data_dir,country,sep=''))

poly.path <- paste0("shapeFiles_gadm")
poly.layer.adm0 <- paste('gadm36', gadm.abbrev,
                         '0', sep = "_")
poly.layer.adm1 <- paste('gadm36', gadm.abbrev,
                         '1', sep = "_")
poly.layer.adm2 <- paste('gadm36', gadm.abbrev,
                         '2', sep = "_")

#poly.path <- paste0(folder.name, poly.file)
poly.adm0 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm0))
poly.adm1 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm1))
poly.adm2 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm2))

proj4string(poly.adm0) <- proj4string(poly.adm1) <- proj4string(poly.adm2)
load(paste0('shapeFiles_gadm/', country, '_Amat.rda'))
load(paste0('shapeFiles_gadm/', country, '_Amat_Names.rda'))


################################################################
#########   Final preprocessing
################################################################  
setwd(paste(data_dir,country,sep=''))

### make sure data are in place ###
load(paste0(country,'_cluster_dat.rda'),
     envir = .GlobalEnv)
load( paste0(poly.path,'/', country, '_Amat.rda'))
load( paste0(poly.path, '/', country, '_Amat_Names.rda'))

### Prepare analysis data set ###

mod.dat<-mod.dat[mod.dat$survey==as.character(survey_year),]
mod.dat$years <- as.numeric(as.character(mod.dat$years))
mod.dat<-mod.dat[as.numeric(mod.dat$years)>=beg.year,]

mod.dat$strata <- NA
mod.dat$country <- as.character(country)

## set benchmark (no benchmark, all 1s)

adj.frame <- expand.grid(years = beg.year:end.year,
                         country = country)
adj.frame$ratio <- 1
adj.varnames <- c("country", "years")

bench.adj <- adj.frame
bench.adj$ratio <- 1.0


################################################################
#########  National level model
################################################################

setwd(paste(res_dir,country,sep=''))

fit.natl <- smoothCluster(data = mod.dat, family = "betabinomial",
                          Amat = NULL, 
                          year_label = c(beg.year:end.year),
                          time.model = "rw2",
                          overdisp.mean = -7.5,
                          overdisp.prec = 0.39,
                          age.groups = levels(mod.dat$age),
                          age.n = c(1,11,12,12,12,12),
                          age.rw.group = c(1,2,3,3,3,3),
                          verbose = FALSE,
                          survey.effect = FALSE)


natl.unweights <- cbind.data.frame(years = c(beg.year:end.year), 
                                    urban = rep(0.5, length(c(beg.year:end.year))))
natl.unweights$rural <- 1- natl.unweights$urban

res.natl <- getSmoothed(inla_mod = fit.natl, 
                        year_range = c(beg.year:end.year),
                        year_label = c(beg.year:end.year),
                        nsim = 1000, 
                        draws = NULL, save.draws = TRUE)

save(fit.natl, file = paste0('Betabinomial/',
                             country, '_rw2_natl.rda'))

save(res.natl, file = paste0('Betabinomial/',
                             country, '_res_rw2_natl.rda'))


################################################################
#########   Admin-1 level model
################################################################

mod.dat$region <- mod.dat$admin1.char
fit.admin1 <- smoothCluster(data = mod.dat, family = "betabinomial",
                            Amat = admin1.mat, 
                            age.strata.fixed.group = c(1,2,3,3,3,3),
                            year_label = c(beg.year:end.year),
                            time.model = "rw2", st.time.model = "ar1",
                            pc.st.slope.u = 1, pc.st.slope.alpha = 0.01,
                            type.st = type.st,
                            bias.adj = bench.adj,
                            bias.adj.by = adj.varnames,
                            survey.effect = FALSE)


## Posterior draws
res.admin1 <- getSmoothed(fit.admin1, 
                          nsim = 1000, 
                          save.draws = TRUE)


## Save results
res.admin1_overall<-res.admin1$overall
admin1_res_merged<-merge(res.admin1_overall, admin1.names, 
                         by.x=c("region"),
                         by.y=c("Internal"))
admin1_res_merged$region<-admin1_res_merged$GADM
class(admin1_res_merged) <- class(res.admin1$overall)
admin1_res_merged_cleaned<-admin1_res_merged[, c("region", "years", "median", "mean", "variance", "lower", "upper")]
write.csv(admin1_res_merged_cleaned,
          paste0(country, '_admin1_U5MR.csv'))

save(fit.admin1, file = paste0("Betabinomial/",
                               country, '_rw2main_randomSlopes_ar1xICAR_admin1.rda'))

save(res.admin1, file = paste0("Betabinomial/",
                               country, '_res_rw2main_randomSlopes_ar1xICAR_admin1.rda'))


################################################################
#########   overall models admin2
################################################################

mod.dat$region <- mod.dat$admin2.char
fit.admin2 <- smoothCluster(data = mod.dat, family = "betabinomial",
                            Amat = admin2.mat, 
                            year_label = c(beg.year:end.year),
                            time.model = "rw2", st.time.model = "ar1",
                            pc.st.slope.u = 1, pc.st.slope.alpha = 0.01,
                            type.st = type.st,
                            bias.adj.by = adj.varnames,
                            survey.effect = FALSE)

## Posterior draws
res.admin2 <- getSmoothed(fit.admin2, 
                          nsim = 1000, 
                          save.draws = TRUE)

## Visualization, save figures: trends
res.admin2_overall<-res.admin2$overall
admin2_res_merged<-merge(res.admin2_overall, admin2.names, 
                         by.x=c("region"),
                         by.y=c("Internal"))
admin2_res_merged$region<-admin2_res_merged$GADM
class(admin2_res_merged) <- class(res.admin2$overall)

## Save results
admin2_res_merged_cleaned<-admin2_res_merged[, c("region", "years", "median", "mean", "variance", "lower", "upper")]
write.csv(admin2_res_merged_cleaned,
          paste0(country, '_admin2_U5MR.csv'))

save(fit.admin2, file = paste0("Betabinomial/",
                               country, '_rw2main_randomSlopes_ar1xICAR_admin2.rda'))

save(res.admin2, file = paste0("Betabinomial/",
                               country, '_res_rw2main_randomSlopes_ar1xICAR_admin2.rda'))

