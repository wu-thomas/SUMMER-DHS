################################################################
#########   Load libraries
################################################################

## Download most recent version of SUMMER from Github
#library(devtools)
#devtools::install_github("bryandmartin/SUMMER",
#                         build_vignettes = F, force = T)


rm(list = ls())
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
library(Rfast)
library(parallel)
library(spdep)
library(geosphere)
library(SimDesign)

# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home_dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")
countries <- countries <- scan(paste0(home_dir, "/countries_implemented.txt"), character(), quote = "")
country <- countries[length(countries)] # retrieve the country being analyzed
info.name <- paste0(country, "_general_info.Rdata")
load(file = paste0(home_dir,'/Info/',country,"/", info.name, sep='')) # load the country info

compare.year <- beg.year: end.year
sd.year <- seq(beg.year + 1, end.year - 1, 3) # year range for smooth direct

################################################################
#########   set directories
################################################################

data.dir <- paste0(home_dir,'/Data/')
res.dir <- paste0(home_dir,'/Results/')

################################################################
#########   load polygon files
################################################################


setwd(paste(data.dir,country,sep=''))

poly.path <- paste0("shapeFiles_gadm")

poly.layer.adm0 <- paste('gadm36', gadm.abbrev,
                         '0', sep = "_")  # specify the name of the national shape file
poly.layer.adm1 <- paste('gadm36', gadm.abbrev,
                         '1', sep = "_") # specify the name of the admin1 shape file
poly.layer.adm2 <- paste('gadm36', gadm.abbrev,
                         '2', sep = "_") # specify the name of the admin2 shape file


poly.adm0 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm0))  # load the national shape file
# use encoding to read special characters
poly.adm1 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm1)) # load the shape file of admin-1 regions

if(sum(grepl(paste('gadm36', gadm.abbrev,
                   '2', sep = "_"), list.files(poly.path))) != 0){ # load the shape file of admin-2 regions
  poly.adm2 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                       layer = as.character(poly.layer.adm2))}

# set coordinate reference system to be equal
if(exists("poly.adm2")){
  proj4string(poly.adm0) <- proj4string(poly.adm1)  <- proj4string(poly.adm2)
}else{
  proj4string(poly.adm0) <- proj4string(poly.adm1)
}

### admin name dictionary
load( paste0(poly.path,'/', country, '_Amat.rda'))
load( paste0(poly.path, '/', country, '_Amat_Names.rda'))


################################################################
#########   create directories
################################################################  

setwd(paste0(res.dir,'/',country))

if(!dir.exists(paths = paste0('Figures'))){
  dir.create(path = paste0('Figures'))
}

if(!dir.exists(paths = paste0('Figures/Comparison'))){
  dir.create(path = paste0('Figures/Comparison'))
}

if(!dir.exists(paths = paste0('Figures/Comparison/National'))){
  dir.create(path = paste0('Figures/Comparison/National'))
}

###################################################################################
#########  Load national level models
###################################################################################

### national yearly direct (direct.natl.yearly)
load(file = paste0('Direct/', country, '_direct_natl_yearly.rda'))  


### national yearly smooth direct (res.natl.yearly)
file.out <- paste0(country, "_res_natl_yearly_SmoothedDirect.rda")
load( file = paste0('Smooth_Direct/', file.out))


###################################################################################
#########  Load admin-1 level models
###################################################################################

### Direct estimate admin1 3-year window (direct.admin1)
load(paste0('Direct/',country,'_direct_admin1.rda'))
direct.admin1<-direct.admin1[direct.admin1$region!='All',]

### smooth direct admin1 3-year window (admin1_sd)
file.out <- paste0(country, "_res_admin1_SmoothedDirect.rda")
load( file = paste0('Smooth_Direct/', file.out))
admin1.sd<-res.admin1$results
admin1.sd.draws<-res.admin1$draws.est

### smooth direct admin1 yearly (sd.admin1.yearly)
file.out <- paste0(country, "_res_admin1_SmoothedDirect_yearly.rda")
load( file = paste0('Smooth_Direct/', file.out))
admin1.sd.yearly<-sd.admin1.yearly$results
admin1.sd.yearly.draws<-sd.admin1.yearly$draws.est


### BB8 admin1 stratified (admin1.strat.BB8)
res.strat.admin1 <- readRDS(paste0('Betabinomial/',
                                country,'_res.strat.admin1.3UR.rds'))
admin1.strat.BB8<-res.strat.admin1$overall



###################################################################################
#########  Load admin-2 level models
###################################################################################

### BB8 admin2 stratified (admin2.strat.BB8)
res.strat.admin2 <- readRDS(paste0('Betabinomial/',
                                   country,'_res.strat.admin2.3UR.rds'))
admin2.strat.BB8<-res.strat.admin2$overall

################################################################
#########   prepare function and dictionary (each row is one admin2 region)
################################################################

# link accross admin1 and admin2
adm_link <- poly.adm2@data
adm_link <- adm_link[,c('NAME_1','NAME_2')]
colnames(adm_link)<-c('admin1_name','admin2_name')

# get unique ID for admin2
adm_link$admin2_idx<- admin2.names$Internal

# get admin1 index
adm1_match<-match(adm_link$admin1_name, admin1.names$GADM)
adm_link$admin1_idx<- admin1.names$Internal[adm1_match]

# not using this because of repeated admin2 name
# adm2_match<-match(adm_link$admin2_name, admin2.names$GADM)#
# adm_link$admin2_idx<- admin2.names$Internal[adm2_match]

# function that calculates population in each admin2 area
pop_adm2<-function(adm2.shp, wp,admin_pop_dat){
  
  # make sure polygons have same crs as population raster
  adm2.shp <- spTransform(adm2.shp, wp@crs)
  
  # list of admin2 regions
  adm2.names <- adm2.shp$NAME_2
  
  # admin 2 level population
  wp.adm2.list <- lapply(1:nrow(adm2.shp), function(x) {
    list(state_id = x, state_raster = mask(crop(wp,adm2.shp[x,]),
                                           adm2.shp[x,]))
  })
  
  # store total population at admin 2
  pop.adm2<-vector()
  for ( j in 1:nrow(adm2.shp)){
    pop_j<-wp.adm2.list[[j]]
    pop.adm2[j]<-sum(values(pop_j$state_raster),na.rm=TRUE)
    
  }
  
  # add admin2 population 
  admin_pop_dat$admin2_pop<-pop.adm2
  # not using matching because repeated admin2 names
  # assume order in admin.link is the same as polygon
  
  return (admin_pop_dat)
}

################################################################
#########  function to get BB8 posterior 
################################################################  

## function to get posterior draws from BB8 
draw_1y_adm<-function(admin_draws, year_num,admin_vec, nsim=1000){
  
  # year_num: year of comparison
  # nsim: number of posterior draws
  # admin_vec: vector of admin index
  # admin_draws: posterior draws (as a list from SUMMER output)
  
  # prepare reference frame for draws 
  # ID corresponds to specific year, region
  draw_ID<-c(1:length(admin_draws))
  draw_year<-vector()
  draw_region<-vector()
  
  for( i in draw_ID){
    tmp_d<-admin_draws[[i]]
    draw_year[i]<-tmp_d$years
    draw_region[i]<-tmp_d$region
  }
  
  draw_ref<-data.frame(id=draw_ID,year=draw_year,
                       region=draw_region)
  
  draw_frame<-matrix( nrow = nsim, ncol = length(admin_vec))
  
  for(i in 1:length(admin_vec)){
    admin_i<-admin_vec[i]
    id_draw_set<-draw_ref[draw_ref$year==year_num&
                            draw_ref$region==admin_i,]$id 
    
    draw_set<-admin_draws[[id_draw_set]]$draws
    
    draw_frame[,i]<-draw_set
    #print(mean(r_frame[,c(admin_i)]))
  }
  
  colnames(draw_frame)<-admin_vec
  
  return(draw_frame)
}



################################################################
#########  prepare admin fractions 
################################################################ 
setwd(paste0(data.dir,'/',country))

adm1.frac.list <- list()
adm2.frac.list <- list()

for(year in compare.year){
  
  print(year)

  
  pop_u5<- raster (paste0('Population/',country.abbrev,'_u5_',year,'_100m','.tif'))

  # admin 2 population fraction 
  adm2_pop<-pop_adm2(adm2.shp=poly.adm2,
                     wp=pop_u5,
                     admin_pop_dat=adm_link)
  
  adm2_pop<-adm2_pop %>% 
    group_by(admin1_idx) %>% 
    mutate(admin1_pop = sum(admin2_pop))
  
  
  # fraction of admin2 w.r.t. admin1
  adm2_pop$admin2_frac<-adm2_pop$admin2_pop/
    adm2_pop$admin1_pop
  
  
  # sanity check, fraction for admin2 in each admin1 sum up to 1
  aggregate(admin2_frac~admin1_idx, data = adm2_pop, FUN = sum)
  
  
  save(adm2_pop, file = paste('worldpop/', 'admin2_pop_frac_', year, '.rda', sep = ''))
  
  # delete the world pop tif files
  # unlink(f_0_name)
  # unlink(f_1_name)
  # unlink(m_0_name)
  # unlink(m_1_name)
}


for (i in 1: length(compare.year)){
  #i = 1
  year = compare.year[i]
  load(paste0( 'worldpop/', 'admin2_pop_frac_', year, '.rda', sep = ''))  
  
  ### admin-2 level population fraction w.r.t. natl and admin-1
  
  adm2.pop <- adm2_pop
  adm2.pop$adm2.adm1.frac <- adm2.pop$admin2_frac
  adm2.pop$adm2.natl.frac <- adm2.pop$admin2_pop/sum(adm2.pop$admin2_pop)
  
  adm2.frac.list [[i]] <- adm2.pop
  ### admin-1 level population fraction w.r.t. natl population
  
  adm1.pop<-adm2_pop[!duplicated(adm2_pop[,c('admin1_idx')]),]
  # create an ordered admin1 list
  match.order = match(paste("admin1", 1: nrow(adm1.pop), 
                            sep = "_"), adm1.pop$admin1_idx)
  adm1.pop = adm1.pop[match.order, ]
  
  adm1.pop<-adm1.pop[,c('admin1_name','admin1_idx','admin1_pop')]
  adm1.pop$admin1_frac<-adm1.pop$admin1_pop/sum(adm1.pop$admin1_pop)
  adm1.frac.list [[i]] <- adm1.pop
  
}

names(adm1.frac.list) <- compare.year
names(adm2.frac.list) <- compare.year

################################################################
#########  prepare national level estimates 
################################################################ 

setwd(paste0(data.dir,'/',country))


### direct estimate draws
natl.direct.yearly.draw = expit(Rfast::rmvnorm(10000, mu = direct.natl.yearly$logit.est,
                                               sigma = direct.natl.yearly$var.est*diag(nrow(direct.natl.yearly))))
natl.direct.yearly = t(apply(natl.direct.yearly.draw, 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(natl.direct.yearly) = c("lower", "median", "upper")

natl.direct.yearly.frame <- as.data.frame(natl.direct.yearly)
natl.direct.yearly.frame$method <- 'natl.direct'
natl.direct.yearly.frame$years <- c(beg.year: end.year)



# Load smoothed direct estimate at national level
natl.sd.frame<-data.frame()

natl.sd.est <- res.natl.yearly[res.natl.yearly$years %in% compare.year, "median"]
natl.sd.lower <- res.natl.yearly[res.natl.yearly$years %in% compare.year, "lower"]
natl.sd.upper <- res.natl.yearly[res.natl.yearly$years %in% compare.year, "upper"]
natl.sd.year <- res.natl.yearly[res.natl.yearly$years %in% compare.year, "years"]


natl.sd.frame<-data.frame(median=natl.sd.est,
                          lower=natl.sd.lower,
                          upper=natl.sd.upper,
                          method='natl.sd',
                          years=natl.sd.year)



# aggregated 3-year smoothed direct estimates of admin-1 level

sd.adm1.to.natl.frame = matrix(NA, nrow = length(sd.year),
                              ncol =  ncol(natl.direct.yearly))

for (i in 1: length(sd.year)){
  # i = 1
  year = sd.year[i]
  adm1.pop <- adm1.frac.list[[as.character(year)]]
  
  tmp.res<-admin1.sd[order(admin1.sd$region.gadm),]
  sd.idx<-which(tmp.res$years.num==year) #assuming admin1_idx has the correct order
  t.sd.draw <- t(admin1.sd.draws)
  t.sd.draw<-t.sd.draw[,sd.idx]
  colnames(t.sd.draw)<-adm1.pop$admin1_idx
  
  adm1.sd.natl.draw <- t.sd.draw %*% adm1.pop$admin1_frac
  sd.adm1.to.natl.frame[i, ] = quantile(adm1.sd.natl.draw, probs = c(0.025, 0.5, 0.975))
}


sd.adm1.to.natl.frame<-as.data.frame(sd.adm1.to.natl.frame)
colnames(sd.adm1.to.natl.frame) = c("lower", "median", "upper")
sd.adm1.to.natl.frame$method <- "aggre.sd.adm1"
sd.adm1.to.natl.frame$years = sd.year


# aggregated yearly smoothed direct estimates of admin-1 level

sd.adm1.yl.to.natl.frame = matrix(NA, nrow = length(compare.year),
                               ncol =  ncol(natl.direct.yearly))

for (i in 1: length(compare.year)){
  # i = 1
  year = compare.year[i]
  adm1.pop <- adm1.frac.list[[as.character(year)]]
  
  tmp.res<-admin1.sd.yearly[order(admin1.sd.yearly$region.gadm),]
  sd.idx<-which(tmp.res$years.num==year) #assuming admin1_idx has the correct order
  t.sd.draw <- t(admin1.sd.yearly.draws)
  t.sd.draw<-t.sd.draw[,sd.idx]
  colnames(t.sd.draw)<-adm1.pop$admin1_idx
  
  adm1.sd.natl.draw <- t.sd.draw %*% adm1.pop$admin1_frac
  sd.adm1.yl.to.natl.frame[i, ] = quantile(adm1.sd.natl.draw, probs = c(0.025, 0.5, 0.975))
}


sd.adm1.yl.to.natl.frame<-as.data.frame(sd.adm1.yl.to.natl.frame)
colnames(sd.adm1.yl.to.natl.frame) = c("lower", "median", "upper")
sd.adm1.yl.to.natl.frame$method <- "aggre.sd.yearly.adm1"
sd.adm1.yl.to.natl.frame$years = compare.year


# aagregated stratified BB8 estimates of admin-1 level

BB8.adm1.to.natl.frame <- matrix(NA, nrow = length(compare.year),
                                 ncol =  ncol(natl.direct.yearly)+1)

for (i in 1: length(compare.year)){
  # i = 1
  year = compare.year[i]
  
  adm1.pop <- adm1.frac.list[[as.character(year)]]
  
  adm1.strat.bb8.draw<-draw_1y_adm(admin_draws=res.strat.admin1$draws.est.overall,
                                   year_num=year,
                                   admin_vec=admin1.names$Internal)
  
  natl.tmp <- adm1.strat.bb8.draw %*% adm1.pop$admin1_frac
  BB8.adm1.to.natl.frame[i, ] = c(year, quantile(natl.tmp, probs = c(0.025, 0.5, 0.975)))
  
}

BB8.adm1.to.natl.frame<-as.data.frame(BB8.adm1.to.natl.frame)
colnames(BB8.adm1.to.natl.frame) = c("years","lower", "median", "upper")
BB8.adm1.to.natl.frame$method <- "aggre.adm1.strat.BB8"


# aggregated stratified BB8 estimates of admin-2 level

BB8.adm2.to.natl.frame <- matrix(NA, nrow = length(compare.year),
                                 ncol =  ncol(natl.direct.yearly)+1)

for (i in 1: length(compare.year)){
  # i = 1
  year = compare.year[i]
  
  adm2.pop <- adm2.frac.list[[as.character(year)]]
  
  adm2.strat.bb8.draw<-draw_1y_adm(admin_draws=res.strat.admin2$draws.est.overall,
                             year_num=year,
                             admin_vec=admin2.names$Internal)
  
  natl.tmp <- adm2.strat.bb8.draw %*% adm2.pop$adm2.natl.frac
  BB8.adm2.to.natl.frame[i, ] = c(year, quantile(natl.tmp, probs = c(0.025, 0.5, 0.975)))
  
}

BB8.adm2.to.natl.frame<-as.data.frame(BB8.adm2.to.natl.frame)
colnames(BB8.adm2.to.natl.frame) = c("years","lower", "median", "upper")
BB8.adm2.to.natl.frame$method <- "aggre.adm2.strat.BB8"


### combine all methods

natl.all<-rbind(natl.direct.yearly.frame,
                sd.adm1.yl.to.natl.frame,
                sd.adm1.to.natl.frame,
                natl.sd.frame,
                BB8.adm1.to.natl.frame,
                BB8.adm2.to.natl.frame)


### save data frame for plotting 
setwd(paste0(res.dir,'/',country))

#saveRDS(natl.all, file = paste0("Comparison/", country, "-natl-compare.rds") )



### plot 

natl.to.plot <- natl.all[natl.all$method %in% c('natl.direct','aggre.adm2.strat.BB8'), ]
natl.to.plot$region <- natl.to.plot$method
natl.to.plot$years.num <- as.numeric(natl.to.plot$years)
natl.to.plot$is.yearly <- FALSE
class(natl.to.plot)<-class(res.strat.admin1$overall)


g1 <- plot(natl.to.plot, plot.CI = TRUE, dodge.width = 0.25, proj_year = end.year + 1, per1000=TRUE) +
  theme(legend.position = 'bottom') + scale_linetype(guide='none')+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = c(beg.year:end.year))+
  ylab('Deaths per 1000 live births')+
  labs(color='Method')+
  scale_colour_discrete(labels = c( "Aggregated BB8 admin-2", "National direct yearly"))


setwd(paste0(res.dir,'/',country))
ggsave(g1, width=8, height = 6, file = paste0("Figures/Comparison/National/", country, "-natl-compare-BB8-direct.pdf"))

ggsave(g1, width=8, height = 6, file = paste0("Figures/Report/", country, "-natl-compare-BB8-direct.pdf"))
ggsave(g1, width=8, height = 6, file = paste0("Figures/Report/", country, "-natl-compare-BB8-direct.tiff"))

