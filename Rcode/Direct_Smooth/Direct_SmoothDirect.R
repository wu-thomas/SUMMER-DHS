################################################################
#########   Load libraries
################################################################

rm(list = ls())
options(warn=-1)
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
library(gridExtra)
library(parallel)
library(spdep)
library(geosphere)



# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

# retrieve user-specified directory
home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")
countries <- countries <- scan(paste0(home.dir, "/countries_implemented.txt"), character(), quote = "") 
country <- countries[length(countries)] # retrieve the country being analyzed
info.name <- paste0(country, "_general_info.Rdata")
load(file = paste0(home.dir,'/Info/',country,"/", info.name, sep='')) # load the country info

################################################################
#########   set sub-directories
################################################################

data.dir <- paste0(home.dir,'/Data/') # set the directory to store the data
res.dir <- paste0(home.dir,'/Results/') # set the directory to store the results (e.g. fitted R objects, figures, tables in .csv etc.)

################################################################
#########   load polygon files
################################################################

setwd(paste(data.dir,country,sep=''))

poly.path <- paste0("shapeFiles_gadm") # specify the folder of the country shape files

poly.layer.adm0 <- paste('gadm36', gadm.abbrev,
                         '0', sep = "_") # specify the name of the national shape file
poly.layer.adm1 <- paste('gadm36', gadm.abbrev,
                         '1', sep = "_") # specify the name of the admin1 shape file
poly.layer.adm2 <- paste('gadm36', gadm.abbrev,
                         '2', sep = "_") # specify the name of the admin2 shape file


poly.adm0 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm0)) # load the national shape file
# use encoding to read special characters
poly.adm1 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm1)) # load the shape file of admin-1 regions

if(sum(grepl(paste('gadm36', gadm.abbrev,
                   '2', sep = "_"), list.files(poly.path))) != 0){
  poly.adm2 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                       layer = as.character(poly.layer.adm2))} # load the shape file of admin-2 regions

# set coordinate reference system to be equal
if(exists("poly.adm2")){
  proj4string(poly.adm0) <- proj4string(poly.adm1)  <- proj4string(poly.adm2)
}else{
  proj4string(poly.adm0) <- proj4string(poly.adm1)
}

 
load( file = paste0(poly.path,'/', country, '_Amat.rda')) # load the adjacency matrix
load( file = paste0(poly.path, '/', country, '_Amat_Names.rda')) # load names of admin1 and admin2 regions


# load the data
load(paste0(country,'_cluster_dat.rda'),
     envir = .GlobalEnv)


# create repositories to store direct estimates and smoothed direct estimates
setwd(res.dir)
if(!dir.exists(paths = paste0(country))){
  dir.create(path = paste0(country))
}
setwd(paste0(res.dir,country))

if(!dir.exists(paths = paste0('Direct'))){
  dir.create(path = paste0('Direct'))
}

if(!dir.exists(paths = paste0('Smooth_Direct'))){
  dir.create(path = paste0('Smooth_Direct'))
}

# load a helper function from the R script getSmoothed.R in the same folder
source(file=paste0(home.dir, '/Rcode/Direct_Smooth/', 'getSmoothed.R'))

################################################################
#########   direct estimates 
################################################################

doBenchmark <- F # no benchmarking in the direct model
folder.name<-country
beg.years <- seq(beg.year,beg.year+6,3) # since we're computing direct estimate for each 3-year period, we compute the first year for each period.
end.years <- beg.years + 2 # compute the last year of each period
periods <- paste(beg.years, end.years, sep = "-") # convert the periods into string


mod.dat$years <- as.numeric(as.character(mod.dat$years)) # convert the years from string into numbers
mod.dat$period <- as.character(cut(mod.dat$years, breaks = c(beg.years, beg.years[length(beg.years)]+5),
                                   include.lowest = T, right = F, labels = periods)) # generate period label 

# After we process the data, we compute for the direct U5MR estimates at national level

mod.dat$v005 <- mod.dat$v005/1e6

births.dat <- mod.dat %>% as.data.frame()
births.dat$died <- births.dat$Y
births.dat$total <- as.numeric(births.dat$total)

# Now we run the function to compute direct U5MR at national level for each 3-year period. 
direct.natl <-  getDirect(births.dat, periods,
                          regionVar = "admin1.char",
                          timeVar = "period", 
                          clusterVar =  "~cluster",
                          ageVar = "age", Ntrials = "total",
                          weightsVar = "v005",national.only = T)
direct.natl$survey <- 1
direct.natl$surveyYears <- births.dat$survey[1]

# Now we run the function to compute direct U5MR at national level for each year.
direct.natl.yearly <-  getDirect(births.dat, beg.year:end.year,
                                 regionVar = "admin1.char",
                                 timeVar = "years", 
                                 clusterVar =  "~cluster",
                                 ageVar = "age", Ntrials = "total",
                                 weightsVar = "v005",national.only = T)
direct.natl.yearly$survey <- 1
direct.natl.yearly$surveyYears <- births.dat$survey[1]

# Note that the argument "timeVar" in "gerDirect" can control the time period of the model, 


save(direct.natl, file = paste0('Direct/', country, '_direct_natl.rda')) # save the national 3-year direct U5MR
save(direct.natl.yearly, file = paste0('Direct/', country, '_direct_natl_yearly.rda')) # save the national yearly direct U5MR


# compute the direct estimates for U5MR for each 3-year period at admin1 level
direct.admin1 <-  getDirect(births.dat, periods,
                            regionVar = "admin1.char",
                            timeVar = "period", 
                            clusterVar =  "~cluster",
                            ageVar = "age", Ntrials = "total",
                            weightsVar = "v005",national.only = F)
direct.admin1$survey <- 1
direct.admin1$surveyYears <- births.dat$survey[1]


save(direct.admin1, file = paste0('Direct/', country, '_direct_admin1.rda'))  # save the admin1 3-year direct U5MR

# compute the direct estimates for U5MR for each year at admin1 level
direct.admin1.yearly <-  getDirect(births.dat, beg.year:end.year,
                              regionVar = "admin1.char",
                              timeVar = "years", 
                              clusterVar =  "~cluster",
                              ageVar = "age", Ntrials = "total",
                              weightsVar = "v005",national.only = F)
direct.admin1.yearly$survey <- 1
direct.admin1.yearly$surveyYears <- births.dat$survey[1]

save(direct.admin1.yearly, file = paste0('Direct/', country, '_direct_admin1_yearly.rda')) # save the admin1 yearly direct U5MR

# compute 3-year direct U5MR for at admin2 level. But this may fail due to the data sparsity at admin2 level.
direct.admin2 <-  getDirect(births.dat, periods,
                            regionVar = "admin2.char",
                            timeVar = "period", 
                            clusterVar =  "~cluster",
                            ageVar = "age", Ntrials = "total",
                            weightsVar = "v005",national.only = F)
direct.admin2$survey <- 1
direct.admin2$surveyYears <- births.dat$survey[1]


save(direct.admin2, file = paste0('Direct/', country, '_direct_admin2.rda')) # save the admin2 3-year direct U5MR


################################################################
#########   smooth direct estimates 
################################################################
setwd(paste0(res.dir,country))

#### load direct estimates ####
load( file = paste0('Direct/', country, '_direct_natl.rda'))
load( file = paste0('Direct/', country, '_direct_natl_yearly.rda'))
load(file = paste0('Direct/', country, '_direct_admin1.rda'))
load(file = paste0('Direct/', country, '_direct_admin2.rda'))


# compute smoothed direct U5MR estimates for 3-year periods at national level
proj.per <- paste(end.year+1, end.year+3, sep = "-") # add the 3-year period to be projected
fit.natl <- fitINLA(direct.natl, geo = NULL, Amat = NULL, # national level model doesn't need to specify adjacency matrix since it would just be 1.
                    year_label = c(periods, proj.per),
                    year_range = c(beg.year, end.year + 3), is.yearly = F) # fit the smoothed direct model. Changing the year label and year range can change the years the estimators to be computed, 
# even for future years where DHS data is not yeat available. But this would lead to less accurate estimates and larger uncertainty level.

res.natl <- getSmoothed(fit.natl, year_range = c(beg.year, end.year+3),
                        year_label = c(periods, proj.per)) # sample for smoothed direct estimates

res.natl$years.num <- seq(beg.year+1, end.year+3, 3)
res.natl$region.gadm <- country
file.out <- paste0(country, "_res_natl_SmoothedDirect.rda")
save(res.natl, file = paste0('Smooth_Direct/',file.out)) # save the national 3-year smoothed direct U5MR

# compute yearly smoothed direct U5MR estimates for at national level 
fit.natl.yearly <- fitINLA(direct.natl.yearly, geo = NULL, Amat = NULL,
                           year_label = as.character(beg.year:(end.year + 3)),
                           year_range = c(beg.year, end.year + 3), is.yearly = F)
res.natl.yearly <- getSmoothed(fit.natl.yearly, year_range = c(beg.year, end.year + 3),
                               year_label = as.character(beg.year:(end.year + 3)))
res.natl.yearly$years.num <- beg.year:(end.year + 3)
res.natl.yearly$region.gadm <- country
file.out <- paste0(country, "_res_natl_yearly_SmoothedDirect.rda")
save(res.natl.yearly, file = paste0('Smooth_Direct/', file.out)) # save the national yearly smoothed direct U5MR


#### Admin-1 Model ####

# compute smoothed direct U5MR estimates for 3-year periods at admin1 level
data.admin1 <- direct.admin1[direct.admin1$region!='All',] # direct.admin1 is a matrix containing all national and admin1 level estimates. Only admin1 estimates are interested here.
fit.admin1 <- fitINLA(data.admin1, geo = poly.adm1, Amat = admin1.mat,
                      year_label = c(periods, proj.per),
                      year_range = c(beg.year, end.year + 3), is.yearly = F)
res.admin1 <- getSmoothed_sd(fit.admin1, Amat = admin1.mat,
                             year_range = c(beg.year, end.year + 3),
                             year_label = c(periods, proj.per),
                             save.draws = TRUE)

res.admin1$results$years.num <- seq(beg.year+1, end.year+3, 3)[match(res.admin1$results$years, c(periods, proj.per))]
res.admin1$results$region.gadm <- admin1.names$GADM[match(res.admin1$results$region, admin1.names$Internal)]

file.out <- paste0(country, "_res_admin1_SmoothedDirect.rda")
save(res.admin1, file = paste0('Smooth_Direct/', file.out))# save the admin1 3-year smoothed direct U5MR


# compute yearly smoothed direct U5MR estimates at admin1 level
direct.admin1.yearly <- direct.admin1.yearly[direct.admin1.yearly$region!='All',]
fit.admin1.yearly <- fitINLA(direct.admin1.yearly, Amat = admin1.mat,
                             year_label = as.character(beg.year:end.year),
                             year_range = c(beg.year, end.year), is.yearly = F)
sd.admin1.yearly <- getSmoothed_sd(fit.admin1.yearly, Amat = admin1.mat,
                                   year_range = c(beg.year, end.year),
                                   year_label = as.character(beg.year:end.year),
                                   save.draws = TRUE)
sd.admin1.yearly$results$region.gadm <- admin1.names$GADM[match(sd.admin1.yearly$results$region, admin1.names$Internal)]


file.out <- paste0(country, "_res_admin1_SmoothedDirect_yearly.rda")

save(sd.admin1.yearly, file = paste0('Smooth_Direct/', file.out)) # save the admin1 yearly smoothed direct U5MR

options(warn=-1)