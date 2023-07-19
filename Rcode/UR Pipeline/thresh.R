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
library(stringdist)
library(openxlsx)

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
pop_dir <- paste0(data_dir,'/Population') # set the directory to store the population surface files

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
#########   load worldpop
################################################################


setwd(pop_dir)

# UNadjusted population counts
worldpop <- raster(paste0(country.abbrev,
                          '_ppp_',frame.year,
                          '_1km_Aggregated_UNadj.tif',sep=''))

################################################################
#########   sampling frame urban proportion at admin1 
################################################################


setwd(data_dir)

# read the .txt or .xlsx file containing urban population fraction at admin1 level.
if (file.exists(paste(country.abbrev, "frame_urb_prop.txt", sep = "_"))){
  frame <- read.delim(paste(country.abbrev, "frame_urb_prop.txt", sep = "_"), 
                      header = FALSE,sep=' ')
}
if (file.exists(paste(country.abbrev, "frame_urb_prop.xlsx", sep = "_"))){
  frame <- read.xlsx(paste(country.abbrev, "frame_urb_prop.xlsx", sep = "_"))
}


# # identify column for fraction (need additional processing in general)
frame[,c(2,4)] <- lapply(frame[,c(2,4)],   ## function to remove comma in numbers
                         function(x){as.numeric(gsub(",", "", x))})
frame$frac <- frame[, 2]/frame[, 4]

# greedy algorithm to match admin names 
adm1.ref <- expand.grid(tolower(frame[, 1]),
                        tolower(admin1.names$GADM)) # Distance matrix in long form
names(adm1.ref) <- c("frame_name","gadm_name")
### string distance,  jw=jaro winkler distance, try 'dl' if not working
adm1.ref$dist <- stringdist(adm1.ref$frame_name,
                            adm1.ref$gadm_name, method="jw") 

greedyAssign <- function(a,b,d){
  x <- numeric(length(a)) # assgn variable: 0 for unassigned but assignable, 
  # 1 for already assigned, -1 for unassigned and unassignable
  while(any(x==0)){
    min_d <- min(d[x==0]) # identify closest pair, arbitrarily selecting 1st if multiple pairs
    a_sel <- a[d==min_d & x==0][1] 
    b_sel <- b[d==min_d & a == a_sel & x==0][1] 
    x[a==a_sel & b == b_sel] <- 1
    x[x==0 & (a==a_sel|b==b_sel)] <- -1
  }
  cbind(a=a[x==1],b=b[x==1],d=d[x==1])
}

match_order<-data.frame(greedyAssign(adm1.ref$frame_name,
                                     adm1.ref$gadm_name,
                                     adm1.ref$dist))

# create reference table 
ref.tab <- admin1.names
ref.tab$matched_name <- frame$V1[match_order$a] ### check!!!
ref.tab$urb_frac <- frame$frac[match_order$a] 




################################################################
#########   admin 1 threshold 
################################################################

## load grid
setwd(data_dir)
load(file='prepared_dat/natl_grid.rda')

# index the grid
urb_dat$index <- c(1:nrow(urb_dat))
adm1_dat <- split( urb_dat , f = urb_dat$admin1 )


# This function computes the urban population threshold for a given admin1 area.
# This is done by keep counting the urban locations until the urban population fraction in the reference table is reached.
thresh_urb<-function(adm_grid,ref_tab){
  
  # sort grid population
  vals <- adm_grid$pop_den
  vals[is.na(vals)] <- 0
  sort.obj <- sort.int(vals, decreasing = TRUE, index.return = TRUE, method = 'shell')
  svals <- sort.obj$x
  svals.int <- sort.obj$ix
  
  # extract cutoff proportion based on admin1
  adm.idx <- adm_grid$admin1.char[1]
  cutoff <- ref_tab[ref_tab$Internal==adm.idx,]$urb_frac
  
  # determine population threshold and urban rural
  csvals <- cumsum(svals)/sum(svals)
  is.urb <- csvals <= cutoff
  org.isurb <- is.urb[invPerm(svals.int)]
  threshold <- min(vals[org.isurb == 1]) #cutoff
  
  # prepare return object (grid with urban/rural)
  adm_grid$threshold <- threshold
  adm_grid$urban <- as.numeric(org.isurb)
  #adm_grid[is.na(adm_grid$pop_den),]$urban<-NA
  
  return(adm_grid)
  
}

urb_list<-lapply(adm1_dat, FUN=thresh_urb,ref_tab=ref.tab)

urb_class <- do.call("rbind", urb_list)


urb_grid <- urb_dat
urb_grid$urb_ind <-NA
urb_grid[urb_class$index,]$urb_ind <- urb_class$urban


urb_surf<-worldpop
values(urb_surf)<-urb_grid$urb_ind


## save reference table along with calculated threshold 
thresh_ref <- urb_class[!duplicated(urb_class[,c('admin1')]),]
ref.tab$threshold <- thresh_ref$threshold # check whether the thresholds are sensible (shouldn't be NA or all 0)
write.xlsx(ref.tab, file='prepared_dat/reference_table.xlsx',
           row.names = FALSE)

################################################################
#########   check classification accuracy based on clusters
################################################################

setwd(data_dir)

load('prepared_dat/crc_dat.rda')
load('prepared_dat/uncrc_dat.rda')

### remove rows with missing covariates, could also build model with missing data
crc_dat<-crc_dat_final[complete.cases(crc_dat_final), ]
uncrc_dat<-uncrc_dat_final[complete.cases(uncrc_dat_final), ]


xy_crc <- as.matrix(crc_dat[c('x','y')])
xy_uncrc <- as.matrix(uncrc_dat[c('x','y')])

# extract the urban/rural prediction
crc_dat$urb_pred<-raster::extract(urb_surf,xy_crc)
uncrc_dat$urb_pred<-raster::extract(urb_surf,xy_uncrc)
pred_crc <- factor( ifelse(crc_dat$urb_pred ==1 ,"urban","rural" ))
pred_crc  <- relevel(pred_crc, "urban") # make sure levels are same 
pred_uncrc <- factor( ifelse(uncrc_dat$urb_pred ==1 ,"urban","rural" ))
pred_uncrc  <- relevel(pred_uncrc, "urban") # make sure levels are same 



### create directory for results
setwd(res_dir)

if(!dir.exists(paths = paste0('UR/'))){
  dir.create(path = paste0('UR/'))
}

if(!dir.exists(paths = paste0('UR/Threshold/'))){
  dir.create(path = paste0('UR/Threshold/'))
}

if(!dir.exists(paths = paste0('UR/U5_fraction/'))){
  dir.create(path = paste0('UR/U5_fraction/'))
}

# compute the confusion to evaluate the accuracy
confmatrix_crc<-caret::confusionMatrix(
  data = pred_crc,
  reference = crc_dat$urban
)

confmatrix_crc
save(confmatrix_crc,file='UR/Threshold/confmatrix_crc.rda')

confmatrix_uncrc<-caret::confusionMatrix(
  data = pred_uncrc,
  reference = uncrc_dat$urban
)

confmatrix_uncrc
save(confmatrix_uncrc,file='UR/Threshold/confmatrix_uncrc.rda')


################################################################
#########  load function for urban fraction
################################################################


get_subnatl_frac<-function(adm.names,adm.idx,wp,poly_file,wp_adm=NULL,
                           urb_vec){
  
  poly_file <- spTransform(poly_file, crs(wp)) # comment out the line if throws error
  
  
  
  if(is.null(wp_adm))
    wp_adm <- lapply(1:nrow(poly_file), function(x) {
      list(state_id = x, state_raster = mask(crop(wp,poly_file[x,]), poly_file[x,]))
    })
  
  pred_surf <- wp
  values(pred_surf)<-urb_vec
  
  urb_adm <- lapply(1:nrow(poly_file), function(x) {
    list(state_id = x, state_raster = mask(crop(pred_surf,poly_file[x,]), poly_file[x,]))
  })
  
  frac_vec<-vector()
  
  for(j in 1:length(adm.names)){
    urb_j<-urb_adm[[j]]
    wp_j<-wp_adm[[j]]
    
    val_urb_j<-values(urb_j$state_raster)
    val_wp_j<-values(wp_j$state_raster)
    
    frac_vec[j]<-sum(val_urb_j*val_wp_j,na.rm=TRUE)/
      sum(val_wp_j,na.rm=TRUE)
  }
  
  subnatl_frac<-data.frame(adm_name=adm.names,adm_idx=adm.idx,urb_frac=frac_vec)
  return(subnatl_frac)
}


################################################################
#########  national U5 urban proportion 
################################################################

years <- c(beg.year:end.year)
natl.u5.urb <- vector()

for ( t in 1:length(years)){
  
  print(t)
  year <- years[t]
  
  # load U5 population at year t
  setwd(pop_dir)
  u5_pop<-raster(paste0(country.abbrev,'_u5_',year,'_1km.tif'))
  
  # national urban fraction for U5 population at year t
  u5_natl <- sum(urb_grid$urb_ind*values(u5_pop),na.rm=TRUE)/
    sum(values(u5_pop),na.rm=TRUE)
  natl.u5.urb[t] <- u5_natl
  
}

natl.urb.weights <- data.frame(years= years, urban=natl.u5.urb)


setwd(res_dir)
saveRDS(natl.urb.weights,paste0('UR/U5_fraction/','natl_urban_weights.rds'))

################################################################
#########  subnatl U5 urban proportion 
################################################################


years <- c(beg.year:end.year)
adm1.weight.frame <- data.frame()
adm2.weight.frame <- data.frame()

for ( t in 1:length(years)){
  
  print(t)
  year <- years[t]
  
  # load U5 population at year t
  setwd(pop_dir)
  u5_pop<-raster(paste0(country.abbrev,'_u5_',year,'_1km.tif'))
  
  # admin1 urban fraction for U5 population at year t
  u5_urb_admin1<-get_subnatl_frac(adm.names = admin1.names$GADM,
                                  adm.idx = admin1.names$Internal,
                                  wp=u5_pop,
                                  poly_file = poly.adm1,
                                  wp_adm = NULL,
                                  urb_vec = urb_grid$urb_ind)
  
  u5_urb_admin1$years <- year
  adm1.weight.frame <- rbind(adm1.weight.frame,u5_urb_admin1)
  
  # admin2 urban fraction for U5 population at year t
  u5_urb_admin2<-get_subnatl_frac(adm.names = admin2.names$GADM,
                                  adm.idx = admin2.names$Internal,
                                  wp=u5_pop,
                                  poly_file = poly.adm2,
                                  wp_adm = NULL,
                                  urb_vec = urb_grid$urb_ind)
  
  u5_urb_admin2$years <- year
  adm2.weight.frame <- rbind(adm2.weight.frame,u5_urb_admin2)
  
  setwd(res_dir)
  
  # save calculated urban fractions
  saveRDS(u5_urb_admin1,file=paste0('UR/U5_fraction/','admin1_',
                                    year, '_urban_frac.rds'))
  saveRDS(u5_urb_admin2,file=paste0('UR/U5_fraction/','admin2_',
                                    year, '_urban_frac.rds'))
  
}


# process admin 1 urban rural weights data frame
adm1.weight.frame <- adm1.weight.frame[,c('adm_idx','years','urb_frac')]
colnames(adm1.weight.frame) <- c('region','years','urban')
adm1.weight.frame$rural <- 1 - adm1.weight.frame$urban

# process admin 2 urban rural weights data frame
adm2.weight.frame <- adm2.weight.frame[,c('adm_idx','years','urb_frac')]
colnames(adm2.weight.frame) <- c('region','years','urban')
adm2.weight.frame$rural <- 1 - adm2.weight.frame$urban

# save weights frames
saveRDS(adm1.weight.frame,paste0('UR/U5_fraction/','admin1_urban_weights.rds'))
saveRDS(adm2.weight.frame,paste0('UR/U5_fraction/','admin2_urban_weights.rds'))
options(warn = 0)
