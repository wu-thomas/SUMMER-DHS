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
library(sqldf)
library(sp)
library(gstat)
library(ggridges)

# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

# retrieve user-specified directory
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

data_dir <- paste0(home_dir,'/Data/')
res_dir <- paste0(home_dir,'/Results/')



################################################################
#########   load shapefiles
################################################################

setwd(paste(data_dir,country,sep=''))

poly.path <- paste0("shapeFiles_gadm") # specify the folder of the country shape files

poly.layer.adm0 <- paste('gadm36', gadm.abbrev,
                         '0', sep = "_") # specify the name of the national shape file
poly.layer.adm1 <- paste('gadm36', gadm.abbrev,
                         '1', sep = "_") # specify the name of the admin1 shape file
poly.layer.adm2 <- paste('gadm36', gadm.abbrev,
                         '2', sep = "_") # specify the name of the admin2 shape file

#poly.path <- paste0(folder.name, poly.file)
poly.adm0 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm0)) # load the national shape file
poly.adm1 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm1)) # load the shape file of admin-1 regions
poly.adm2 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm2)) # load the shape file of admin-2 regions

# set coordinate reference system to be equal
proj4string(poly.adm0) <- proj4string(poly.adm1) <- proj4string(poly.adm2)

load(paste0('shapeFiles_gadm/', country, '_Amat.rda'))
load(paste0('shapeFiles_gadm/', country, '_Amat_Names.rda'))

################################################################
#########   load fitted BB8 models and results
################################################################

setwd(paste0(res_dir,'/',country))


### load fitted objects 
if (FALSE){ # don't need here 
  fit.strat.natl <- readRDS( paste0('Betabinomial/',
                                    country,'_fit.strat.natl.3UR.rds'))
  
  fit.strat.admin1 <- readRDS(paste0('Betabinomial/',
                                     country,'_fit.strat.admin1.3UR.rds'))
  
  fit.strat.admin2 <- readRDS(paste0('Betabinomial/',
                                     country,'_fit.strat.admin2.3UR.rds'))
}



### load results 
res.strat.natl <- readRDS( paste0('Betabinomial/',
                                  country,'_res.strat.natl.3UR.rds'))

res.strat.admin1 <- readRDS(paste0('Betabinomial/',
                                   country,'_res.strat.admin1.3UR.rds'))

res.strat.admin2 <- readRDS(paste0('Betabinomial/',
                                   country,'_res.strat.admin2.3UR.rds'))



### load direct national estimates  

load(paste0('Direct/', country, '_direct_natl.rda'))

################################################################
#########   create directories 
################################################################

setwd(paste0(res_dir,'/',country))

if(!dir.exists(paths = paste0('Betabinomial/Postsamp'))){
  dir.create(path = paste0('Betabinomial/Postsamp'))
}


################################################################
#########   posterior sampling
################################################################
# The codes below retrieves the posterior sample of BB8 estimates. The sample will
# be used to produce some plots. (e.g. ridge plot)
setwd(paste0(res_dir,'/',country))
# posterior sampling parameters
years_vt <- c(beg.year:end.year)
n_years <- length(years_vt)
n_samp <- 1000
admin_vt <- c("admin1","admin2") ## admin level
level_vt <- c(1,2) ## polygon level
strat <- 'strat' ## using stratified model

for(admin in admin_vt){
  
  res_obj <- get(paste0("res.strat.", admin ))
  
  draws_list <- res_obj$draws.est.overall
  n_admin <- length(draws_list)/n_years 
  
  postsamp_mt_list <- vector(mode = "list", length = n_years)
  
  for (i in 1:n_years){
    # i <- 1
    postsamp_mt_list[[i]]$years <- years_vt[i]
    postsamp_mt <- matrix(0, nrow = n_admin, ncol = n_samp)
    
    for(j in 1:n_admin){
      # j <- 1
      postsamp_mt[j, ] <- draws_list[[n_years*(j-1)+i]]$draws
    }
    
    postsamp_mt_list[[i]]$postsamp_mt <- postsamp_mt
  }
  
  save(postsamp_mt_list, 
       file = paste0("Betabinomial/Postsamp/", country, "_", admin,
                     "_", strat, "_postsamp.RData"))
  
}

################################################################
######### load posterior draws 
################################################################


# posterior sampling parameters
years_vt <- c(beg.year:end.year)
n_years <- length(years_vt)
n_samp <- 1000
admin_vt <- c("admin1","admin2") ## admin level
level_vt <- c(1,2) ## polygon level
strat <- 'strat' ## using stratified model

setwd(paste0(res_dir,'/',country))

admin = "admin2"
map_shp <- poly.adm2
admin_name_dt <- as.data.table(get(paste0(admin, ".names")))
admin_name_dt$toPlot <- paste0(map_shp$NAME_2)
admin_name_dt$admin1.name <-  map_shp$NAME_1


### load draws

load( file = paste0("Betabinomial/Postsamp/", country, "_", admin,
                    "_", strat, "_postsamp.RData"))





################################################################
#########   load functions
################################################################

setwd(home_dir)
source(paste0("Rcode", "/Visualization/Calculate_exceed_prob.R"))




################################################################
#########   create directories 
################################################################

### create repositories for results
setwd(res_dir)
if(!dir.exists(paths = paste0(country))){
  dir.create(path = paste0(country))
}

setwd(paste0(res_dir,'/',country))

if(!dir.exists(paths = paste0('Figures'))){
  dir.create(path = paste0('Figures'))
}

if(!dir.exists(paths = paste0('Figures/Report'))){
  dir.create(path = paste0('Figures/Report'))
}



################################################################
#########   Figure 1a : admin-2 U5MR
################################################################

# We plot the U5MR estimates of the latest year fitted by beta-binomial model at admin2 level on the map of the given country.

# prepare results
res.admin2_overall<-res.strat.admin2$overall
admin2_res_merged<-merge(res.admin2_overall, admin2.names, 
                         by.x=c("region"),
                         by.y=c("Internal"))
admin2_res_merged$region<-admin2_res_merged$GADM
admin2_res_merged$width <- admin2_res_merged$upper- admin2_res_merged$lower
class(admin2_res_merged) <- class(res.strat.admin2$overall)


# admin2 level maps, last year
last_year <- subset(admin2_res_merged, years == end.year)
g1a <- mapPlot(last_year, geo = poly.adm2, by.data = "region", by.geo = "NAME_2", variable = "median", 
              is.long=FALSE, per1000=TRUE, removetab=TRUE, legend.label = "U5MR", direction = -1, 
              size= 0.1)

g1a <- g1a+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=16),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='U5MR (deaths per 1000 live births)',
                                label.position = "bottom"))

ggsave(g1a, width=8, height = 10, file = paste0("Figures/Report/", country, "-", end.year, "-admin2-map.tiff"))
ggsave(g1a, width=8, height = 10, file = paste0("Figures/Report/", country, "-", end.year, "-admin2-map.pdf"))



################################################################
#########   Figure 1b : admin-2 U5MR CI width
################################################################

# We plot the width of credible of our U5MR estimates of the latest year fitted by beta-binomial model at admin2 level on the map of the given country.

g1b <- mapPlot(last_year, geo = poly.adm2, by.data = "region", by.geo = "NAME_2", variable = "width",
              is.long=FALSE, per1000=TRUE,  removetab=TRUE,legend.label = "Width of 95% CI", direction = -1,
              size= 0.1)

g1b <- g1b + scale_fill_viridis_c("Width of 95% CI",option="B",direction=-1)+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=16),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))


ggsave(g1b, width=8, height = 10, file = paste0("Figures/Report/", country,"-", end.year, "-wid-95CI-admin2-map.tiff"))
ggsave(g1b, width=8, height = 10, file = paste0("Figures/Report/", country,"-", end.year, "-wid-95CI-admin2-map.pdf"))


################################################################
#########   Figure 1c : admin-2 U5MR over year
################################################################

# We plot the width of credible of our U5MR estimates for all years fitted by beta-binomial model at admin2 level on the map of the given country.

g1c <- mapPlot(admin2_res_merged, geo = poly.adm2, by.data = "region", by.geo = "NAME_2", variable = "years",
              values = "median", is.long=TRUE, per1000=TRUE, legend.label = "U5MR", direction = -1, ncol = 3, 
              size= 0.1) +
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=15),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=15),
         strip.text = element_text(size = 14))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='U5MR (deaths per 1000 live births)',
                                label.position = "bottom"))

ggsave(g1c, width=10, height = 10, file = paste0("Figures/Report/", country, "-all-years-admin2-map.pdf"))
ggsave(g1c, width=10, height = 10, file = paste0("Figures/Report/", country, "-all-years-admin2-map.tiff"))


################################################################
#########   Figure 1d : admin-2 U5MR CI width over year
################################################################

# We plot the width of credible of our U5MR estimates for all years fitted by beta-binomial model at admin2 level on the map of the given country.

g1d <- mapPlot(admin2_res_merged, geo = poly.adm2, by.data = "region", by.geo = "NAME_2", variable = "years",
              values = "width", is.long=TRUE, per1000=TRUE, legend.label = "Width of 95% CI", 
              direction = -1, ncol = 3)+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=15),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=15),
         strip.text = element_text(size = 14))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))



ggsave(g1d, width=10, height = 10, file = paste0("Figures/Report/", country, "-all-years-wid-95CI-admin2-map.pdf"))
ggsave(g1d, width=10, height = 10, file = paste0("Figures/Report/", country, "-all-years-wid-95CI-admin2-map.tiff"))



################################################################
#########   Figure 1e : admin-1 U5MR over year
################################################################

# We plot our U5MR estimates for all years fitted by beta-binomial model at admin1 level on the map of the given country.

# prepare results
res.admin1_overall<-res.strat.admin1$overall
admin1_res_merged<-merge(res.admin1_overall, admin1.names, 
                         by.x=c("region"),
                         by.y=c("Internal"))
admin1_res_merged$region<-admin1_res_merged$GADM
admin1_res_merged$width <- admin1_res_merged$upper- admin1_res_merged$lower
class(admin1_res_merged) <- class(res.strat.admin1$overall)


g1e <- mapPlot(admin1_res_merged, geo = poly.adm1, by.data = "region", by.geo = "NAME_1", variable = "years",
               values = "median", is.long=TRUE, per1000=TRUE, legend.label = "U5MR", direction = -1, ncol = 3, 
               size= 0.1) +
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=15),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=15),
         strip.text = element_text(size = 14))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='U5MR (deaths per 1000 live births)',
                                label.position = "bottom"))

ggsave(g1e, width=10, height = 10, file = paste0("Figures/Report/", country, "-all-years-admin1-map.pdf"))
ggsave(g1e, width=10, height = 10, file = paste0("Figures/Report/", country, "-all-years-admin1-map.tiff"))


################################################################
#########   Figure 1f : admin-1 U5MR CI width over year
################################################################

# We plot the width of credible of our U5MR estimates for all years fitted by beta-binomial model at admin1 level on the map of the given country.

g1f <- mapPlot(admin1_res_merged, geo = poly.adm1, by.data = "region", by.geo = "NAME_1", variable = "years",
               values = "width", is.long=TRUE, per1000=TRUE, legend.label = "Width of 95% CI", 
               direction = -1, ncol = 3)+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=15),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=15),
         strip.text = element_text(size = 14))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))

ggsave(g1f, width=10, height = 10, file = paste0("Figures/Report/", country, "-all-years-wid-95CI-admin1-map.pdf"))
ggsave(g1f, width=10, height = 10, file = paste0("Figures/Report/", country, "-all-years-wid-95CI-admin1-map.tiff"))

################################################################
#########   Figure 2 : admin-1 U5MR trend plot
################################################################

# We plot our U5MR estimates versus years fitted by beta-binomial model at admin1 level.

# prepare results
res.admin1_overall<-res.strat.admin1$overall
admin1_res_merged<-merge(res.admin1_overall, admin1.names, 
                         by.x=c("region"),
                         by.y=c("Internal"))
admin1_res_merged$region<-admin1_res_merged$GADM
admin1_res_merged$width <- admin1_res_merged$upper- admin1_res_merged$lower
class(admin1_res_merged) <- class(res.strat.admin1$overall)

### Admin 1 level trend plot
g2 <- plot(admin1_res_merged, plot.CI = TRUE, dodge.width = 0.5, proj_year = end.year + 1, per1000=TRUE) +
  theme(legend.position = 'bottom') + scale_linetype(guide='none')+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = c(beg.year:end.year))+
  ylab('Deaths per 1000 live births')+
  theme(text = element_text(size=16))

ggsave(g2, width=8, height = 10, file = paste0("Figures/Report/", country, "-trends-admin1.tiff"))
ggsave(g2, width=8, height = 10, file = paste0("Figures/Report/", country, "-trends-admin1.pdf"))



################################################################
#########   Figure 3 : admin-2 U5MR trend plot, no legends
################################################################

# We plot our U5MR estimates versus years fitted by beta-binomial model at admin2 level.

range_adm2 <- range(c(admin2_res_merged$median, admin2_res_merged$median)) * 1000
admin2_trend_plot <- admin2_res_merged
admin2_trend_plot$region <- res.strat.admin2$overall$region
g3 <- plot(admin2_trend_plot, plot.CI = FALSE, dodge.width = 0.5, proj_year = end.year + 1, per1000=TRUE) +
  theme(legend.position = 'none') + scale_linetype(guide='none')+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = c(beg.year:end.year))+
  ylim(range_adm2)+
  ylab('Deaths per 1000 live births')+
  theme(text = element_text(size=16))

ggsave(g3, width=10, height = 8, file = paste0("Figures/Report/", country, "-trends-admin2.pdf"))
ggsave(g3, width=10, height = 8, file = paste0("Figures/Report/", country, "-trends-admin2.tiff"))


################################################################
#########   Figure 4 : Exceedance probability 
################################################################

# We compute the probability our U5MR estimates fitted by beta-binoimial model exceed the national direct estimate for the latest year

# set parameters
plot.year <- sd.year[length(sd.year)] + 1 # last year of the last 3-year window
first.year <- sd.year[1] + 1             # last year of the first 3-year windows


# prepare grouping 
grouping <- expand.grid(years.1 = c(first.year, plot.year), 
                        region.1 = unique(res.strat.admin2$overall$region))
grouping$years.2[grouping$years.1 == first.year] <- paste0(first.year-2,'-',first.year)
grouping$years.2[grouping$years.1 == plot.year] <- paste0(plot.year-2,'-',plot.year)

grouping$strata.1 = "strata_all"
grouping$region.2 <- "All"
grouping$strata.2 <- "strata_all"


### calculate exceedance probability
prob2 <- exceedProb(est1 = res.strat.admin2, est2 = direct.natl, 
                    grouping = grouping)


### last year 
prob.last.year <- subset(prob2, years == plot.year)
prob.last.year <- merge(prob.last.year, admin2.names, 
                        by.x = c("region"), 
                        by.y = c("Internal"))

g4 <- mapPlot(prob.last.year, geo = poly.adm2, by.data = "GADM", by.geo = "NAME_2", variable = "prob",
              is.long=FALSE, removetab=TRUE,legend.label = "Exceedance Probability", direction = -1,
              size= 0.1)

g4 <- g4+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=16),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='Exceedance Probability',
                                label.position = "bottom"))

ggsave(g4, width=10, height = 8, file = paste0("Figures/Report/", country, "-", plot.year, "-exceed-map.pdf"))
ggsave(g4, width=10, height = 8, file = paste0("Figures/Report/", country, "-", plot.year, "-exceed-map.tiff"))



################################################################
#########   posterior sampling for ridge, rank, TCP
################################################################

# some parameters about the posterior sample of our U5MR estimator
years_vt <- c(beg.year:end.year)
n_years <- length(years_vt)
n_samp <- 1000
admin_vt <- c("admin1","admin2") ## admin level
level_vt <- c(1,2) ## polygon level
strat <- 'strat' ## using stratified model


setwd(paste0(res_dir,'/',country))

for(admin in admin_vt){
  
  res_obj <- get(paste0("res.strat.", admin ))
  
  draws_list <- res_obj$draws.est.overall
  n_admin <- length(draws_list)/n_years 
  
  postsamp_mt_list <- vector(mode = "list", length = n_years)
  
  for (i in 1:n_years){
    # i <- 1
    postsamp_mt_list[[i]]$years <- years_vt[i]
    postsamp_mt <- matrix(0, nrow = n_admin, ncol = n_samp)
    
    for(j in 1:n_admin){
      # j <- 1
      postsamp_mt[j, ] <- draws_list[[n_years*(j-1)+i]]$draws
    }
    
    postsamp_mt_list[[i]]$postsamp_mt <- postsamp_mt
  }
  
  save(postsamp_mt_list, 
       file = paste0("Betabinomial/Postsamp/", country, "_", admin,
                     "_", strat, "_postsamp.RData"))
  
  
}


################################################################
#########   Figure 5 : Rank Plot 
################################################################

# We produce the rank plot here, which gives the probability of the rank of the U5MR estimates among all the admin areas.
setwd(paste0(res_dir,'/',country))

admin = "admin2"
map_shp <- poly.adm2
admin_name_dt <- as.data.table(get(paste0(admin, ".names")))
admin_name_dt$toPlot <- paste0(map_shp$NAME_2)
admin_name_dt$admin1.name <-  map_shp$NAME_1


load( file = paste0("Betabinomial/Postsamp/", country, "_", admin,
                    "_", strat, "_postsamp.RData"))

year <- end.year

cond <- lapply(postsamp_mt_list, function(x){x$years == year})
postsamp_mt <- postsamp_mt_list[unlist(cond)][[1]]$postsamp_mt

# prepare the posterior sample for rank plot

rank_mt <- apply(postsamp_mt, 2, rank)

pred_dt <- admin_name_dt
pred_dt[, "ID"] <- 1:nrow(pred_dt)
pred_dt[, "avg_rank"] <- apply(rank_mt, 1, mean)
pred_dt[, "low_rank"] <- apply(rank_mt, 1, min)
pred_dt[, "up_rank"] <- apply(rank_mt, 1, max)

pred_dt_order <- pred_dt[order(avg_rank)]


### plotting 

## pdf format 


# We make the rank plot of 5 admin2 areas having the highest possible U5MR estimates.

pdf(paste0("Figures/Report/", country,'-',year, "-rank-top5.pdf"),
     width=6, height = 10)

rank_cut <- round (nrow(pred_dt_order)/2)

par(mar = c(3, 1, 2.5, 1), mfrow = c(5, 1))

pred_dt_order <- pred_dt[order(avg_rank)]

for (i in 1:5){
  
  id <- pred_dt_order[i, ID]
  name <- pred_dt_order[i, toPlot]
  
  rank_vt <- rank_mt[id, ]
  rank_vt <- ifelse(rank_vt <= rank_cut, rank_vt, rank_cut)
  
  avg_rank <- pred_dt_order[i, avg_rank]
  
  ranktable <- as.data.table(table(rank_vt))
  ranktable <- merge(data.table(rank =  as.character(1:rank_cut)), ranktable, 
                     by.x = "rank", by.y = "rank_vt", all.x = T)
  ranktable[, "rank" := as.integer(rank)]
  ranktable <- ranktable[order(rank)]
  ranktable[is.na(N), "N"] <- 0
  
  barplot(ranktable$N, width = 0.825, xlim = c(rank_cut, 0), 
          xlab = "", ylab = "",
          main = paste0(name, "\nER = ", format(round(avg_rank, 1), nsmall = 1)),
          xaxt = "n", yaxt = "n", col = "#31a354", border = F,
          cex.main = 1.5)
    axis(1, at = rank_cut:1-0.5, labels = as.character(c(paste0(rank_cut,'+'),(rank_cut-1):1)), tick = F,
         cex.axis=1.25)
}
dev.off()


### tiff format

tiff(paste0("Figures/Report/", country,'-',year, "-rank-top5.tiff"),
    width=6, height = 10,res = 300,units = 'in')


par(mar = c(3, 1, 2.5, 1), mfrow = c(5, 1))

pred_dt_order <- pred_dt[order(avg_rank)]

for (i in 1:5){
  
  id <- pred_dt_order[i, ID]
  name <- pred_dt_order[i, toPlot]
  
  rank_vt <- rank_mt[id, ]
  rank_vt <- ifelse(rank_vt <= rank_cut, rank_vt, rank_cut)
  
  avg_rank <- pred_dt_order[i, avg_rank]
  
  ranktable <- as.data.table(table(rank_vt))
  ranktable <- merge(data.table(rank =  as.character(1:rank_cut)), ranktable, 
                     by.x = "rank", by.y = "rank_vt", all.x = T)
  ranktable[, "rank" := as.integer(rank)]
  ranktable <- ranktable[order(rank)]
  ranktable[is.na(N), "N"] <- 0
  
  barplot(ranktable$N, width = 0.825, xlim = c(rank_cut, 0), 
          xlab = "", ylab = "",
          main = paste0(name, "\nER = ", format(round(avg_rank, 1), nsmall = 1)),
          xaxt = "n", yaxt = "n", col = "#31a354", border = F,
          cex.main = 1.5)
  axis(1, at = rank_cut:1-0.5, labels = as.character(c(paste0(rank_cut,'+'),(rank_cut-1):1)), tick = F,
       cex.axis=1.25)
}
dev.off()


# We make the rank plot of 5 admin2 areas having the lowest possible U5MR estimates.

rank_cut <- nrow(pred_dt_order)-15
  
  
pdf(paste0("Figures/Report/", country,'-',year, "-rank-bottom5.pdf"),
    width=6, height = 10)


par(mar = c(3, 1, 2.5, 1), mfrow = c(5, 1))

pred_dt_order <- pred_dt[order(-avg_rank)]

for (i in 1:5){
  
  id <- pred_dt_order[i, ID]
  name <- pred_dt_order[i, toPlot]
  
  rank_vt <- rank_mt[id, ]
  rank_vt <- ifelse(rank_vt >= rank_cut, rank_vt, rank_cut)
  
  avg_rank <- pred_dt_order[i, avg_rank]
  
  ranktable <- as.data.table(table(rank_vt))
  ranktable <- merge(data.table(rank = as.character(rank_cut:nrow(pred_dt_order))), ranktable, 
                     by.x = "rank", by.y = "rank_vt", all.x = T)
  ranktable[, "rank" := as.integer(rank)]
  ranktable <- ranktable[order(rank)]
  ranktable[is.na(N), "N"] <- 0
  
  barplot(ranktable$N, width = 0.825, 
          xlim = c((nrow(pred_dt_order)-rank_cut+1), 0), xlab = "", ylab = "",
          main = paste0(name, "\nER = ", format(round(avg_rank, 1), nsmall = 1)),
          xaxt = "n", yaxt = "n", col = "#31a354", border = F,
          cex.main = 1.5)
  axis(1, at = (nrow(pred_dt_order)-rank_cut+1):1-0.5, labels = 
         as.character(c(nrow(pred_dt_order):(rank_cut+1),paste0(rank_cut,'-'))), tick = F,
       cex.axis=1.25)
}
dev.off()


### tiff format

tiff(paste0("Figures/Report/", country,'-',year, "-rank-bottom5.tiff"),
     width=6, height = 10,res = 300,units = 'in')


par(mar = c(3, 1, 2.5, 1), mfrow = c(5, 1))

pred_dt_order <- pred_dt[order(-avg_rank)]

for (i in 1:5){
  
  id <- pred_dt_order[i, ID]
  name <- pred_dt_order[i, toPlot]
  
  rank_vt <- rank_mt[id, ]
  rank_vt <- ifelse(rank_vt >= rank_cut, rank_vt, rank_cut)
  
  avg_rank <- pred_dt_order[i, avg_rank]
  
  ranktable <- as.data.table(table(rank_vt))
  ranktable <- merge(data.table(rank = as.character(rank_cut:nrow(pred_dt_order))), ranktable, 
                     by.x = "rank", by.y = "rank_vt", all.x = T)
  ranktable[, "rank" := as.integer(rank)]
  ranktable <- ranktable[order(rank)]
  ranktable[is.na(N), "N"] <- 0
  
  barplot(ranktable$N, width = 0.825, 
          xlim = c((nrow(pred_dt_order)-rank_cut+1), 0), xlab = "", ylab = "",
          main = paste0(name, "\nER = ", format(round(avg_rank, 1), nsmall = 1)),
          xaxt = "n", yaxt = "n", col = "#31a354", border = F,
          cex.main = 1.5)
  axis(1, at = (nrow(pred_dt_order)-rank_cut+1):1-0.5, labels = 
         as.character(c(nrow(pred_dt_order):(rank_cut+1),paste0(rank_cut,'-'))), tick = F,
       cex.axis=1.25)
}
dev.off()






################################################################
#########   Figure 6 : Ridge Plot
################################################################

# We make ridge plots here, which shows the posterior distribution of the U5MR estimator.

# prepare posterior draws 

data_plot_dt <- NULL

year <- end.year

cond <- lapply(postsamp_mt_list, function(x){x$years == year})
postsamp_mt <- postsamp_mt_list[unlist(cond)][[1]]$postsamp_mt
  
data_plot_dt_year <- data.table(Year = year, 
                                  Internal = rep(admin_name_dt[, Internal], 1000),
                                  GADM = rep(admin_name_dt[, GADM], 1000),
                                  toPlot = rep(admin_name_dt[, toPlot], 1000),
                                  admin1.name = rep(admin_name_dt[, admin1.name], 1000),
                                  U5MR = as.numeric(postsamp_mt))
  
data_plot_dt <- rbind(data_plot_dt, data_plot_dt_year)

manual.col <- colorRampPalette(viridis_pal(direction = -1)(10))(10)[c(1,3,6,9)]# 1000) ## color palette

# prepare results
res.admin2_overall<-res.strat.admin2$overall
admin2_res_merged<-merge(res.admin2_overall, admin2.names, 
                         by.x=c("region"),
                         by.y=c("Internal"))
admin2_res_merged$region<-admin2_res_merged$GADM
admin2_res_merged$width <- admin2_res_merged$upper- admin2_res_merged$lower
class(admin2_res_merged) <- class(res.strat.admin2$overall)

### prepare plotting data frame for top 10 and bottom 10 (if # of admin-2 areas>=20)
### if not, then plot all together
num_admin2 <- dim(admin2.names)[1]

# rank by U5MR for most recent year

### >20 admin-2 regions
if(num_admin2>20){

  # order admin-2 by posterior median  
  last_year <- subset(admin2_res_merged, years == end.year)
  last_year <- data.table(last_year)
  
  # get names for top10 and bottom 10 
  ridge_top10_name <- last_year[order(median)]$GADM[1:10]
  ridge_bottom10_name <- last_year[order(-median)]$GADM[1:10]
  
  ridge.max <- max(data_plot_dt$U5MR)*1000
  ridge.min <- min(data_plot_dt$U5MR)*1000
  
  # plot top 10
  top10_plot_dt <- data_plot_dt[GADM %in% ridge_top10_name, ]
  
  top10_plot_dt_1 <- top10_plot_dt
  top10_plot_dt_1[, "U5MR_med" := median(U5MR), by = c("Year", "Internal")]
  top10_plot_dt_1_order <- top10_plot_dt_1[order(U5MR_med, U5MR)]
  top10_plot_order <- unique(top10_plot_dt_1_order$toPlot)
  
  top10_plot_dt[, Area := factor(toPlot, levels = rev(top10_plot_order))]
  top10_plot_dt[, U5MRperc := as.numeric(U5MR*1000)]
  top10_plot_dt$rank <- 'Top 10'
  
  g6.top <- ggplot(top10_plot_dt, aes(x = U5MRperc, y = Area)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = max(15/length(top10_plot_order),3), size = 0.3) +
    scale_fill_gradientn(colours = manual.col,limits=c(0,ridge.max),
                         name = "U5MR") +
    theme(axis.text.y = element_text(size = 16)) +
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    #theme(strip.background = element_blank(), strip.text = element_blank())+
    ggtitle(paste0('Top 10 admin-2', ' regions within ',country))+
    xlab("")+
    theme(text = element_text(size=16))+
    guides(fill = guide_colourbar(title.position = "top",
                                  title.hjust = .5,
                                  title='U5MR (deaths per 1000 live births)',
                                  label.position = "bottom"))+
    theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=16),
           legend.key.width = unit(2,'cm'),legend.title = element_text(size=16))
  
  ggsave(g6.top, width = 12, height = 8 , file = paste0("Figures/Report/", country,'-',year, "-ridge-top10.tiff"))
  ggsave(g6.top,  width = 12, height = 8 ,file = paste0("Figures/Report/", country,'-',year, "-ridge-top10.pdf"))
  
  
  # plot bottom 10
  bottom10_plot_dt <- data_plot_dt[GADM %in% ridge_bottom10_name, ]
  
  bottom10_plot_dt_1 <- bottom10_plot_dt
  bottom10_plot_dt_1[, "U5MR_med" := median(U5MR), by = c("Year", "Internal")]
  bottom10_plot_dt_1_order <- bottom10_plot_dt_1[order(U5MR_med, U5MR)]
  bottom10_plot_order <- unique(bottom10_plot_dt_1_order$toPlot)
  
  bottom10_plot_dt[, Area := factor(toPlot, levels = rev(bottom10_plot_order))]
  bottom10_plot_dt[, U5MRperc := as.numeric(U5MR*1000)]
  bottom10_plot_dt$rank <- 'Bottom 10'
  
  g6.bottom <- ggplot(bottom10_plot_dt, aes(x = U5MRperc, y = Area)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = max(15/length(bottom10_plot_order),3), size = 0.3) +
    scale_fill_gradientn(colours = manual.col,limits=c(0,ridge.max),
                         name = "U5MR") +
    theme(axis.text.y = element_text(size = 16)) +
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    #theme(strip.background = element_blank(), strip.text = element_blank())+
    ggtitle(paste0('Bottom 10 admin-2', ' regions within ',country))+
    xlab("")+
    theme(text = element_text(size=16))+
    guides(fill = guide_colourbar(title.position = "top",
                                  title.hjust = .5,
                                  title='U5MR (deaths per 1000 live births)',
                                  label.position = "bottom"))+
    theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=16),
           legend.key.width = unit(2,'cm'),legend.title = element_text(size=16))
  
  
  ggsave(g6.bottom, width = 12, height = 8 , file = paste0("Figures/Report/", country,'-',year, "-ridge-bottom10.tiff"))
  ggsave(g6.bottom,  width = 12, height = 8 ,file = paste0("Figures/Report/", country,'-',year, "-ridge-bottom10.pdf"))
  
  
  # altogether, 
  bottom_top <- rbind(top10_plot_dt,bottom10_plot_dt)
  bottom_top$rank <- factor(bottom_top$rank, levels = c("Top 10","Bottom 10"))
  
  
  g6.together.long <- ggplot(bottom_top, aes(x = U5MRperc, y = Area)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = max(15/length(bottom10_plot_order)/2,3), size = 0.3) +
    scale_fill_gradientn(colours = manual.col,
                         name = "U5MR") +
    theme(axis.text.y = element_text(size = 16)) +
    facet_grid(rank ~ . ,scales = "free") + 
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    #theme(strip.background = element_blank(), strip.text = element_blank())+
    ggtitle(paste0('Top and bottom 10 admin-2', ' regions within ',country))+
    xlab("")+
    theme(text = element_text(size=16))+
    guides(fill = guide_colourbar(title.position = "top",
                                  title.hjust = .5,
                                  title='U5MR (deaths per 1000 live births)',
                                  label.position = "bottom"))+
    theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=16),
           legend.key.width = unit(2,'cm'),legend.title = element_text(size=16),
           legend.box.margin=margin(-15,0,0,0),
           legend.margin=margin(0,0,0,0))+
    theme(panel.spacing = unit(1.5, "lines"))
  
  ggsave(g6.together.long, width = 10, height = 10 , file = paste0("Figures/Report/", country,'-',year, "-ridge-bottom-top10.tiff"))
  ggsave(g6.together.long,  width = 10, height = 10 ,file = paste0("Figures/Report/", country,'-',year, "-ridge-bottom-top10.pdf"))
  
  
  
  g6.together.wide <- ggplot(bottom_top, aes(x = U5MRperc, y = Area)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = max(15/length(bottom10_plot_order)/2,3), size = 0.3) +
    scale_fill_gradientn(colours = manual.col,
                         name = "U5MR") +
    theme(axis.text.y = element_text(size = 16)) +
    facet_wrap( rank~ . ,scales = "free") + 
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    #theme(strip.background = element_blank(), strip.text = element_blank())+
    ggtitle(paste0('Top and bottom 10 admin-2', ' regions within ',country))+
    xlab("")+
    theme(text = element_text(size=16))+
    guides(fill = guide_colourbar(title.position = "top",
                                  title.hjust = .5,
                                  title='U5MR (deaths per 1000 live births)',
                                  label.position = "bottom"))+
    theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=16),
           legend.key.width = unit(2,'cm'),legend.title = element_text(size=16))
  
  ggsave(g6.together.wide, width = 12, height = 8 , file = paste0("Figures/Report/", country,'-',year, "-ridge-bottom-top10-wide.tiff"))
  ggsave(g6.together.wide,  width = 12, height = 8 ,file = paste0("Figures/Report/", country,'-',year, "-ridge-bottom-top10-wide.pdf"))
  
  
  
}


### <=20 admin-2 regions

if(num_admin2<21){

  # range of U5MR
  ridge.max <- max(data_plot_dt$U5MR)*1000
  ridge.min <- min(data_plot_dt$U5MR)*1000
  
  # plot all

  data_plot_dt_1 <- data_plot_dt
  data_plot_dt_1[, "U5MR_med" := median(U5MR), by = c("Year", "Internal")]
  data_plot_dt_1_order <- data_plot_dt_1[order(U5MR_med, U5MR)]
  data_plot_order <- unique(data_plot_dt_1_order$toPlot)
  
  data_plot_dt_1[, Area := factor(toPlot, levels = rev(data_plot_order))]
  data_plot_dt_1[, U5MRperc := as.numeric(U5MR*1000)]

  g6.all <- ggplot(data_plot_dt_1, aes(x = U5MRperc, y = Area)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = max(15/length(data_plot_order),3), size = 0.3) +
    scale_fill_gradientn(colours = manual.col,limits=c(0,ridge.max),
                         name = "U5MR") +
    theme(axis.text.y = element_text(size = 16)) +
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    #theme(strip.background = element_blank(), strip.text = element_blank())+
    ggtitle(paste0('All admin-2', ' regions within ',country))+
    xlab("")+
    theme(text = element_text(size=16))+
    guides(fill = guide_colourbar(title.position = "top",
                                  title.hjust = .5,
                                  title='U5MR (deaths per 1000 live births)',
                                  label.position = "bottom"))+
    theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=16),
           legend.key.width = unit(2,'cm'),legend.title = element_text(size=16))
  
  ggsave(g6.all, width = 12, height = max(5,length(data_plot_order)/15*10) , file = paste0("Figures/Report/", country,'-',year, "-ridge-all.tiff"))
  ggsave(g6.all,  width = 12, height = max(5,length(data_plot_order)/15*10) ,file = paste0("Figures/Report/", country,'-',year, "-ridge-all.pdf"))
  
}
################################################################
#########   Figure 7 : TCP plots
################################################################

# We make true classification probability plot, which colors different admin regions based on the range their U5MR lies in.

setwd(paste0(res_dir,'/',country))

### parameter setup 
K_vt <- c(2, 3, 4)
year_vt <- c(end.year:end.year)
admin_level_dt <- data.table(Admin = admin_vt,
                             level = level_vt)

### load function
get_measure_dt <- function(postsamp_mt, pred_dt, grp_thresh){
  # postsamp_mt <- state_postsamp_mt
  # pred_dt <- pred_dt_state
  # grp_thresh <- c(0, 0.5, 1)
  
  # create group lookup table based on thresholds
  n_grp <- length(grp_thresh)-1
  grp_dt <- data.table(grp = 1:n_grp,
                       lower = grp_thresh[1:(n_grp)],
                       upper = grp_thresh[2:(n_grp+1)])
  
  n_area <- nrow(postsamp_mt)
  n_postsamp <- ncol(postsamp_mt)
  
  
  #### assign group based on model max post prob ####
  
  
  grp_cnt_mt <- matrix(0, nrow = n_area, ncol = n_grp)
  for (i in 1:n_grp){
    grp_cnt_mt[, i] <- apply(postsamp_mt, 1, 
                             function(x){sum(x > grp_dt[i, lower] & x <= grp_dt[i, upper])})
  }
  
  DF <- data.frame(grp_cnt_mt)
  DT <- data.table(value = unlist(DF, use.names=FALSE), 
                   colid = 1:nrow(DF), 
                   rowid = rep(names(DF), each=nrow(DF)))
  setkey(DT, colid, value)
  
  grp_cnt_dt_max <- as.data.table(DF)
  grp_cnt_dt_max[, "grp"] <- DT[J(unique(colid), DT[J(unique(colid)), value, mult="last"]), rowid, mult="first"]
  grp_cnt_dt_max[, "TCP" := 0]
  for (i in 1:n_grp){
    idx_grp <- which(grp_cnt_dt_max[, grp] == paste0("X", i))
    grp_cnt_dt_max[idx_grp, "TCP"] <- grp_cnt_dt_max[idx_grp, paste0("X", i), with = F]/n_postsamp
  }
  
  pred_dt[, "grp"] <- grp_cnt_dt_max[, grp]
  pred_dt[, "TCP"] <- grp_cnt_dt_max[, TCP]
  
  return(pred_dt)
}


### prepare results and plot

admin_name_dt <- as.data.table(get(paste0(admin, ".names")))

load( file = paste0("Betabinomial/Postsamp/", country, "_", admin,
                    "_", strat, "_postsamp.RData"))

year <- end.year

cond <- lapply(postsamp_mt_list, function(x){x$years == year})
postsamp_mt <- postsamp_mt_list[unlist(cond)][[1]]$postsamp_mt
  

for(K in K_vt){
  
  intv <- 1/K
  quant_vt <- seq(0, 1, intv)
  quant_val_vt <- (quant_vt[1:K]+quant_vt[2:(K+1)])/2
  
  post_med_vt <- as.numeric(postsamp_mt) # apply(postsamp_mt, 1, median)
  
  L_vt <- quantile(post_med_vt, quant_vt)
  L_vt[1] <- quantile(post_med_vt, 0.01)
  L_vt[K+1] <- quantile(post_med_vt, 0.99)
  L_val_vt <- quantile(post_med_vt, quant_val_vt)
  L_dt <- data.table(grp = paste0("X", 1:K),
                     grp_low = L_vt[1:K],
                     grp_up = L_vt[2:(K+1)],
                     grp_val = L_val_vt)
  
  # set plotting names (different for admin 1 and 2)
  if (admin == "admin1") {
    namesToPlot <- map_shp$NAME_1
  } else if (admin == "admin2") {
    namesToPlot <- paste0(map_shp$NAME_2,
                          ",\n", map_shp$NAME_1)
  }
  
  measure_dt <- get_measure_dt(postsamp_mt = postsamp_mt,
                               pred_dt = data.table(adm_name = paste0(admin, "_", 1:nrow(postsamp_mt)),
                                                    adm_name_toPlot = namesToPlot),
                               grp_thresh = L_vt)
  measure_dt <- merge(measure_dt, L_dt, by = "grp", all.x = T)
  measure_dt[, "K"] <- K
  
  measure_dt[, "Internal" := adm_name]
  measure_dt <- merge(measure_dt, admin_name_dt, by = "Internal", all.x = T)
  measure_dt[, paste0("NAME_", admin_level_dt[Admin == admin, level]) := GADM]
  
  save(measure_dt, L_dt,
       file = paste0("Betabinomial/Postsamp/", country, "_", admin,
                     "_", 
                     strat, "_Y", year,
                     "_K", K, "_measure.RData"))
  
  # set map names to be unique for merging (add in Internal names)
  map_shp$NAME_INTERNAL <- as.character(admin_name_dt[, Internal])
  
  shp_plot <- merge(map_shp, measure_dt,
                    by.x = "NAME_INTERNAL",
                    by.y = "adm_name")
  
  #manual.col <- rev(colorRampPalette(brewer.pal(8, "RdYlGn"))(length(L_dt$grp_val))) # 1000)
  manual.col <- colorRampPalette(viridis_pal(direction = -1)(10))(10)[c(1,3,6,9)]# 1000)
  manual.col[4] <- "navyblue"
  # color.match <- manual.col[1:length(L_dt$grp_val)] # round(sort(L_dt$grp_val)*1000)]
  if(K == 2){
    color.match <- manual.col[c(1,4)]
  }else if(K == 3){
    color.match <- manual.col[c(1,3,4)]
  }else if(K == 4){
    color.match <- manual.col
  }
  lookup_dt <- data.table(grp_val = sort(L_dt$grp_val),
                          col = color.match)
  shp_plot <- merge(shp_plot, lookup_dt, by = "grp_val", duplicateGeoms = TRUE)
  shp_plot$grp_val <- as.factor(shp_plot$grp_val)
  col_regions <- as.vector(lookup_dt[grp_val %in% shp_plot$grp_val, col])
  
  labelat <- sort(unique(c(L_dt$grp_low, 
                           L_dt$grp_up)))
  labeltext <- format(round(labelat*1000, 2), nsmall = 2)
  
  #### measure map ####
  
  pdf(paste0("Figures/Report/",
             country, "_", admin,
             "_", strat, "_Y", year,
             "_K", K, "_measuremap.pdf"),
      width = 8, height = 8)
  sp_plot <- spplot(shp_plot, zcol = "grp_val",
                    main = list(label=paste0("ATCP = ", format(round(mean(measure_dt$TCP), 2), nsmall = 2)),
                                cex=1.5),
                    xlab = "", ylab = "",
                    # sp.layout = list(scale_bar, text_x, text_y),
                    col.regions = col_regions,
                    colorkey = list(col = color.match,
                                    at = labelat,
                                    labels = list(at = labelat, labels = labeltext,
                                                  cex=1.5)))
  print(sp_plot)
  dev.off()
  
  
  tiff(paste0("Figures/Report/",
             country, "_", admin,
             "_", strat, "_Y", year,
             "_K", K, "_measuremap.tiff"),
      width = 8, height = 8,res = 300,units = 'in')
  sp_plot <-  spplot(shp_plot, zcol = "grp_val",
                     main = list(label=paste0("ATCP = ", format(round(mean(measure_dt$TCP), 2), nsmall = 2)),
                                 cex=1.5),
                     xlab = "", ylab = "",
                     # sp.layout = list(scale_bar, text_x, text_y),
                     col.regions = col_regions,
                     colorkey = list(col = color.match,
                                     at = labelat,
                                     labels = list(at = labelat, labels = labeltext,
                                                   cex=1.5)))
  print(sp_plot)
  dev.off()
  

}



################################################################
#########   Map of direct estimates vs. smoothed direct
################################################################

# We plot the 3-year aggregated estimate on the map at admin1 level.
setwd(paste0(res_dir,'/',country))

last_window = paste0(sd.year[length(sd.year)] - 1,'-',sd.year[length(sd.year)] + 1)
### Direct estimate admin1 most recent 3-year window (direct.admin1)
load(paste0('Direct/',country,'_direct_admin1.rda'))
direct.admin1<-direct.admin1[direct.admin1$region!='All'& direct.admin1$years == last_window,]

direct.admin1.merged<-merge(direct.admin1, admin1.names, 
                         by.x=c("region"),
                         by.y=c("Internal"))
direct.admin1.merged$region<-direct.admin1.merged$GADM

### Smoothed direct estimate admin1 most recent 3-year window (res.admin1)
load(paste0('Smooth_Direct/',country,'_res_admin1_SmoothedDirect.rda'))
sd.admin1 <- res.admin1$results
sd.admin1<-sd.admin1[sd.admin1$region!='All'& sd.admin1$years == last_window,]


### set range 

ymax <- ceiling(  max(c(sd.admin1$median,direct.admin1$mean))*1000 )
ymin <- floor(  min(c(sd.admin1$median,direct.admin1$mean)) *1000)

### Map direct estimates 
g9a <- mapPlot(direct.admin1.merged, geo = poly.adm1, by.data = "region", by.geo = "NAME_1", variable = "mean", 
               is.long=FALSE, per1000=TRUE, removetab=TRUE, legend.label = "U5MR", direction = -1, 
               size= 0.1,ylim= c(ymin,ymax))
g9a <- g9a +
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=16),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='U5MR (deaths per 1000 live births)',
                                label.position = "bottom"))


g9b <- mapPlot(sd.admin1, geo = poly.adm1, by.data = "region.gadm", by.geo = "NAME_1", variable = "median", 
               is.long=FALSE, per1000=TRUE, removetab=TRUE, legend.label = "U5MR", direction = -1, 
               size= 0.1,ylim= c(ymin,ymax))
g9b <- g9b +
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),legend.text=element_text(size=16),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='U5MR (deaths per 1000 live births)',
                                label.position = "bottom"))

ggsave(g9a, width=8, height = 10, file = paste0("Figures/Report/", country, "-",last_window, "-admin1-direct-map.tiff"))
ggsave(g9a, width=8, height = 10, file = paste0("Figures/Report/", country, "-",last_window, "-admin1-direct-map.pdf"))


ggsave(g9b, width=8, height = 10, file = paste0("Figures/Report/", country, "-",last_window, "-admin1-smoothed-direct-map.tiff"))
ggsave(g9b, width=8, height = 10, file = paste0("Figures/Report/", country, "-",last_window, "-admin1-smoothed-direct-map.pdf"))


### scatter plot 
match.order <- match(sd.admin1$region,direct.admin1$region)
direct.admin1 <- direct.admin1[match.order,]
range_min <- min(c(direct.admin1$mean,sd.admin1$median))*1000
range_max <- max(c(direct.admin1$mean,sd.admin1$median))*1000

plot.dat <- data.frame(direct.est=direct.admin1$mean*1000, smoothed.est=sd.admin1$median*1000)

g9c <- ggplot(plot.dat, aes(x=direct.est, y=smoothed.est)) +
  geom_point(size=1.5,col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.7))+xlab('Direct estimates')+
  ylab('Smoothed direct estimates')+
  ggtitle(paste('Smoothed direct vs. direct, admin-1,', country, paste0(sd.year[length(sd.year)] - 1,'-',sd.year[length(sd.year)] + 1))) +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_abline(intercept = 0, slope = 1)+
  xlim(range_min,range_max)+  ylim(range_min,range_max)+coord_fixed()+
  theme(text = element_text(size=12))


ggsave(g9c, width=6, height = 6, file = paste0("Figures/Report/", country, "-",last_window, "-admin1-direct-compare.tiff"))
ggsave(g9c, width=6, height = 6, file = paste0("Figures/Report/", country, "-",last_window, "-admin1-direct-compare.pdf"))



################################################################
#########   save stratified (final) estimates as excel 
################################################################

setwd(paste0(res_dir,'/',country))


### save admin-1 estimates
res.admin1_overall<-res.strat.admin1$overall
admin1_res_merged<-merge(res.admin1_overall, admin1.names, 
                         by.x=c("region"),
                         by.y=c("Internal"))
admin1_res_merged$region<-admin1_res_merged$GADM
admin1_res_merged$width <- admin1_res_merged$upper- admin1_res_merged$lower
class(admin1_res_merged) <- class(res.strat.admin1$overall)

admin1_res_merged_cleaned<-admin1_res_merged[, c("region", "years", "median", "mean", "variance", "lower", "upper")]

write.csv(admin1_res_merged_cleaned,
          paste0(country, '_final_3UR_admin1_U5MR.csv'))

### save admin-2 estimates

### prepare results
res.admin2_overall<-res.strat.admin2$overall
admin2_res_merged<-merge(res.admin2_overall, admin2.names, 
                         by.x=c("region"),
                         by.y=c("Internal"))
admin2_res_merged$region<-admin2_res_merged$GADM
admin2_res_merged$width <- admin2_res_merged$upper- admin2_res_merged$lower
class(admin2_res_merged) <- class(res.strat.admin2$overall)

admin2_res_merged_cleaned<-admin2_res_merged[, c("region", "years", "median", "mean", "variance", "lower", "upper")]
write.csv(admin2_res_merged_cleaned,
          paste0(country, '_final_3UR_admin2_U5MR.csv'))
