# USER INPUT REQUIRED AT LINE 79~83! PLEASE INPUT AS INSTRUCTED BY  THE COMMENTS!

################################################################
#########   Load libraries
################################################################

# Download most recent version of SUMMER from Github
# library(devtools)
# devtools::install_github("bryandmartin/SUMMER",
#                          build_vignettes = F, force = T)
# 
# # Install the stable version of INLA
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

rm(list = ls())
library(utils)
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
library(sf)
library(spdep)
library(readstata13)
library(geosphere)
library(mapproj)



# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

# retrieve user-specified directory
home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-2)], collapse = "/")

################################################################
#########   set sub-directories
################################################################

data.dir <- paste0(home.dir,'/Data/') # set the directory to store the data
res.dir <- paste0(home.dir,'/Results/') # set the directory to store the results (e.g. fitted R objects, figures, tables in .csv etc.)

################################################################
#########   set parameters
################################################################

setwd(data.dir) # set the working directory to the folder of Data

# Files info (For those lines with ### xxx ### above, please fill in as commented)
countries <- countries <- scan(paste0(home.dir, "/countries_implemented.txt"), character(), quote = "")
country <- countries[length(countries)]

# retrieve gadm file names
gadm.abbrev <- strsplit(list.files(path = paste0(data.dir, country, "/shapeFiles_gadm"))[[1]][1], "_")[[1]][2]
country.abbrev <- tolower(gadm.abbrev)           # lower the country gadm abbreviation 
poly.path <- paste0(country,"/shapeFiles_gadm")  # specify the folder of the country shape files

# # retrieve the folder containing the DHS data and the name of the DHS data file inside, separated by "/" ###
# dhsStata.folder.path <- strsplit(list.dirs(path = paste0(data.dir, country, "/dhsStata"), recursive = F), "/")[[1]]
# dhsStata.folder <- dhsStata.folder.path[length(dhsStata.folder.path)]
# dhsStata.file <- paste0(dhsStata.folder, "/", dhsStata.folder, ".dta")

# retrieve the file name containing the DHS GPS data ###
dhsFlat.folder <- strsplit(list.dirs(path = paste0(data.dir, country, "/dhsFlat"), recursive = F), "/")[[1]]
dhsFlat.file <- dhsFlat.folder[length(dhsFlat.folder)]


### please fill in the following information ####


dhsStata.file = "SNBR8BFL/SNBR8BFL.dta" # the path to the file containing DHS data
beg.year = 2010   # the first year of the interest. In the DHS report, we considered 9 years before the recent DHS report year
end.year = 2019   # the last year of interest. In the DHS report, we use the year of the recent DHS report.
survey_year<-2019 # year of the DHS survey
frame.year<-2013  # year of the national census frame (Google country + census might give you the year of the most recent census)

#################################################

type.st =  4      # type of space-time interaction in the models. We use type 4 interaction for the models in the DHS report
info.name <- paste0(country, "_general_info.Rdata")
save.image(file = paste0(home.dir,'/Info/',country,"/", info.name, sep=''))

################################################################
#########   load polygon files
################################################################
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


################################################################
#########   Create Adj Matrix
################################################################ 

# Adjacency matrix is a symmetric matrix with each entry of 1's or 0's indicating if the two administrative regions are adjacent. 
# Each row or column represents an administrative region.
# The codes below generates spatial adjacency matrix based on the spatial polygon file.

if(exists("poly.adm1")){ # create the adjacency matrix for admin1 regions.
  admin1.mat <- poly2nb(SpatialPolygons(poly.adm1@polygons))
  admin1.mat <- nb2mat(admin1.mat, zero.policy = TRUE)
  colnames(admin1.mat) <- rownames(admin1.mat) <- paste0("admin1_", 1:dim(admin1.mat)[1])
  admin1.names <- data.frame(GADM = poly.adm1@data$NAME_1,
                             Internal = rownames(admin1.mat))
}else{
  message("There is no Admin1 polygon file.")
}
if(exists("poly.adm2")){  # create the adjacency matrix for admin2 regions.
  admin2.mat <- poly2nb(SpatialPolygons(poly.adm2@polygons))
  admin2.mat <- nb2mat(admin2.mat, zero.policy = TRUE)
  colnames(admin2.mat) <- rownames(admin2.mat) <- paste0("admin2_", 1:dim(admin2.mat)[1])
  admin2.names <- data.frame(GADM = poly.adm2@data$NAME_2,
                             Internal = rownames(admin2.mat))
}else{
  message("There is no Admin2 polygon file.")
}

save(admin1.mat, admin2.mat, file = paste0(poly.path,'/', country, '_Amat.rda')) # save the admin1 and admin2 adjacency matrix
save(admin1.names, admin2.names, file = paste0(poly.path, '/', country, '_Amat_Names.rda')) # save the admin1 and admin2 names


################################################################
#########   Process DHS data
################################################################ 

# The codes below first loads the raw DHS data, then it assigns the GPS coordinates to each sampling cluster and admin regions where
# the sampling is conducted and assigns the admin regions where the clusters are located.

# read DHS data
dat.tmp <- getBirths(filepath = paste0(country,'/dhsStata/',dhsStata.file),
                     surveyyear = survey_year,
                     year.cut = seq(beg.year, end.year + 1, 1),
                     strata = c("v022"), compact = T)

# retrieve the some columns of the full data
dat.tmp <- dat.tmp[ ,c("v001", "v024", "time", "total",
                       "age", "v005", "v025", "strata", "died")]

# specify the name of DHS GPS file, which contains the GPS coordinates of the sampling cluster where the data is sampled
points.path <- paste0(country, "/dhsFlat/", dhsFlat.file)
points <- readOGR(dsn = path.expand(points.path), # read the GPS file
                  layer = as.character(dhsFlat.file))

# detect points in the DHS GPS file with mis-specified coordinates and remove them if any
wrong.points <- which(points@data$LATNUM == 0.0 & points@data$LONGNUM == 0.0)
if(!is.null(dim(wrong.points))){message("There are wrong GPS points: (Longitude, Latitude) = (0, 0)")}

# remove wrong points in the data if any
dat.tmp <- dat.tmp[!(dat.tmp$v001 %in% points@data$DHSCLUST[wrong.points]),]
points@data$DHSCLUST[wrong.points] %in% unique(dat.tmp$v001)

# add the column for GPS coordinate in the data
dat.tmp$LONGNUM <- dat.tmp$LATNUM <- NA
for(i in 1:dim(points)[1]){
  dat.tmp$LATNUM[dat.tmp$v001 == points@data$DHSCLUST[i]] <- points@data$LATNUM[i] # assign latitude to DHS cluster location
  dat.tmp$LONGNUM[dat.tmp$v001 == points@data$DHSCLUST[i]] <- points@data$LONGNUM[i] # assign longitude to DHS cluster location
}

# remove missing points in the data if any
miss <- which(dat.tmp$LATNUM == 0 & dat.tmp$LONGNUM == 0)
if(length(miss != 0)){
  dat.tmp <- dat.tmp[-miss,]
}

message("\n Assigned LAT & LONG")


# assign admin regions based on coordinates and polygon files
adm1.ind <- exists("poly.adm1")
adm2.ind <- exists("poly.adm2")

points.frame <- as.data.frame(dat.tmp[,c("LONGNUM", "LATNUM")]) # retrieve GPS coordinates where data is sampled.
points.frame <- SpatialPoints(points.frame) # convert the GPS coordinates into "sp" object.
if(adm2.ind){
  poly.over.adm2 <- SpatialPolygons(poly.adm2@polygons)
  proj4string(points.frame) <- proj4string(poly.over.adm2) <- 
    proj4string(poly.adm2)  <- 
    proj4string(poly.adm1)  
  admin2.key <- over(points.frame, poly.over.adm2)
  miss.frame.adm2 <- unique(points.frame@coords[which(is.na(admin2.key)),])
  
  if(dim(miss.frame.adm2)[1] != 0){
    miss.poly.adm2 <- dist2Line( miss.frame.adm2, poly.over.adm2)
    
    for(i in 1:dim(miss.poly.adm2)[1]){
      long.ids <- which(points.frame@coords[,c("LONGNUM")] %in% miss.frame.adm2[i,1])
      lat.ids <- which(points.frame@coords[,c("LATNUM")] %in% miss.frame.adm2[i,2])
      ids <- intersect(long.ids, lat.ids)
      admin2.key[ids] <- rep(miss.poly.adm2[i, 'ID'], length(ids))
    }
  }
  
  dat.tmp$admin2 <- admin2.key
  dat.tmp$admin2.char <- paste0("admin2_", admin2.key)
  dat.tmp$admin2.name <- as.character(poly.adm2@data$NAME_2)[admin2.key]
}else{
  dat.tmp$admin2 <- dat.tmp$admin2.name <- NA
  message("There is no Admin2 polygon to assign points to.")
}

if(adm1.ind){
  poly.over.adm1 <- SpatialPolygons(poly.adm1@polygons)
  proj4string(points.frame) <- proj4string(poly.over.adm1) <- 
    proj4string(poly.adm1) 
  admin1.key <- over(points.frame, poly.over.adm1)
  miss.frame.adm1 <- unique(points.frame@coords[which(is.na(admin1.key)),])
  
  if(dim(miss.frame.adm1)[1] != 0){
    miss.poly.adm1 <- dist2Line( miss.frame.adm1, poly.over.adm1)
    
    for(i in 1:dim(miss.poly.adm1)[1]){
      long.ids <- which(points.frame@coords[,c("LONGNUM")] %in% miss.frame.adm1[i,1])
      lat.ids <- which(points.frame@coords[,c("LATNUM")] %in% miss.frame.adm1[i,2])
      ids <- intersect(long.ids, lat.ids)
      admin1.key[ids] <- rep(miss.poly.adm1[i, 'ID'], length(ids))
    }
  }
  
  dat.tmp$admin1 <- admin1.key
  dat.tmp$admin1.char <- paste0("admin1_", admin1.key)
  dat.tmp$admin1.name <- as.character(poly.adm1@data$NAME_1)[admin1.key]
}else{
  dat.tmp$admin2 <- dat.tmp$admin2.name <- NA
  message("There is no Admin1 polygon to assign points to.")
}  

if(FALSE){
  check <- dat.tmp$strata
  check <- gsub(" - rural", "", check)
  check <- gsub(" - urban", "", check)
  inconsist <- which(check != tolower(dat.tmp$admin1.name))
  table(check[inconsist], dat.tmp$admin1.name[inconsist])
  unique(dat.tmp[which(check == "neno" & dat.tmp$admin1.name == "Balaka"), "v001"])
}

# prepare and save the raw data ###
dat.tmp <- dat.tmp[,c("v001", "age", "time", "total", "died", "v005", 
                      "strata", "v025", "LONGNUM", "LATNUM",
                      "admin1", "admin2", "admin1.char", "admin2.char", "admin1.name", "admin2.name")]
colnames(dat.tmp) <- c("cluster", "age", "years", "total",
                       "Y", "v005", "strata", "urban", "LONGNUM", "LATNUM",
                       "admin1", "admin2", "admin1.char", "admin2.char", "admin1.name", "admin2.name")
dat.tmp$survey<-survey_year
dat.tmp$survey.id<-1

mod.dat<-dat.tmp
save(mod.dat, file = paste0(country, '/', country, '_cluster_dat.rda'))


################################################################
#########   Final preprocessing
################################################################  

# try to load the data and make sure all the data processed above have been correctly stored
load(paste0(country,'/',country,'_cluster_dat.rda'),
     envir = .GlobalEnv)
load( paste0(poly.path,'/', country, '_Amat.rda'))
load( paste0(poly.path, '/', country, '_Amat_Names.rda'))

### Prepare analysis data set ###

mod.dat$years <- as.numeric(as.character(mod.dat$years))
mod.dat$strata <- NA # setting strata to be NA will directs to unstratified models in the following pipelines
mod.dat$country <- as.character(country)


## set benchmark (no benchmark, all 1s)

adj.frame <- expand.grid(years = beg.year:end.year,
                         country = country)
adj.frame$ratio <- 1
adj.varnames <- c("country", "years")

bench.adj <- adj.frame
bench.adj$ratio <- 1.0


### create repositories for results
setwd(res.dir)
if(!dir.exists(paths = paste0(country))){
  dir.create(path = paste0(country))
}

setwd(paste0(res.dir,'/',country))

if(!dir.exists(paths = paste0('Betabinomial'))){ # this folder will store the 
  dir.create(path = paste0('Betabinomial'))
}

if(!dir.exists(paths = paste0('Figures'))){
  dir.create(path = paste0('Figures'))
}

