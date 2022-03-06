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
res_dir <- paste0(home_dir,'/Results/') # set the directory to store the results (e.g. fitted R objects, figures, tables in .csv etc.)


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
#########   load fitted BB8 models and results
################################################################

setwd(paste0(res_dir,country))


### load fitted objects 
fit.strat.natl <- readRDS( paste0('Betabinomial/',
                                  country,'_fit.strat.natl.3UR.rds'))

fit.strat.admin1 <- readRDS(paste0('Betabinomial/',
                                   country,'_fit.strat.admin1.3UR.rds'))

fit.strat.admin2 <- readRDS(paste0('Betabinomial/',
                                   country,'_fit.strat.admin2.3UR.rds'))



### load results 
res.strat.natl <- readRDS( paste0('Betabinomial/',
                                  country,'_res.strat.natl.3UR.rds'))

res.strat.admin1 <- readRDS(paste0('Betabinomial/',
                                   country,'_res.strat.admin1.3UR.rds'))

res.strat.admin2 <- readRDS(paste0('Betabinomial/',
                                   country,'_res.strat.admin2.3UR.rds'))



################################################################
#########   create directories 
################################################################

setwd(paste0(res_dir,'/',country))

if(!dir.exists(paths = paste0('Figures'))){
  dir.create(path = paste0('Figures'))
}

if(!dir.exists(paths = paste0('Figures/Comparison'))){
  dir.create(path = paste0('Figures/Comparison'))
}

if(!dir.exists(paths = paste0('Figures/Comparison/UR'))){
  dir.create(path = paste0('Figures/Comparison/UR'))
}


if(!dir.exists(paths = paste0('Figures'))){
  dir.create(path = paste0('Figures'))
}

if(!dir.exists(paths = paste0('Figures/Report'))){
  dir.create(path = paste0('Figures/Report'))
}





################################################################
#########  Function for urban/rural trend and hazard/odds ratio
################################################################

getDiag_UR <- function(inla_mod, CI = 0.95,sampAll=NULL){
  
  year_label = inla_mod$year_label
  
  lower <- (1 - CI) / 2
  upper <- 1 - lower
  
  label <- rep(year_label, length(inla_mod$age.rw.group))
  group <- rep(inla_mod$age.groups, each = length(year_label))
  
  
  fixed <- rownames(inla_mod$fit$summary.fixed)
  fixed <- c(fixed, "time.struct")
  select <- list()
  for(i in 1:length(fixed)){
    select[[i]] <- 0
    names(select)[i] <- fixed[i]
  }
  
  # posterior samples 
  if(is.null(sampAll)){
    sampAll <- INLA::inla.posterior.sample(n = 1e3, 
                                           result = inla_mod$fit, intern = TRUE, 
                                           selection = select, verbose = FALSE)        	
  }
  
  #inla_mod$fit$marginals.random$time.struct 
  re <- grep("time.struct", rownames(sampAll[[1]]$latent)) # random effect
  fe <- grep("age", rownames(sampAll[[1]]$latent)) # fixed effect
  
  
  if(!is.null(inla_mod$age.rw.group)) {
    expand <-  max(inla_mod$age.rw.group) / 
      length(inla_mod$age.rw.group) 
  }
  
  struct.all <- matrix(0, length(re)/expand, length(sampAll))
  
  ages <- gsub( ":rural", "", inla_mod$age.groups)
  ages <- gsub( ":urban", "", ages)
  urs <- gsub(".*:", "", inla_mod$age.groups)
  
  # store posterior samples in a matrix
  for(j in 1:length(sampAll)){
    print(j)
    struct.onesamp <- NULL
    
    for(i in 1:length(inla_mod$age.groups)){
      
      where <- (inla_mod$age.rw.group[i] - 1) * length(year_label)+
        c(1:length(year_label))
      
      
      where.fe1 <- grep(paste0("intercept", ages[i]), rownames(sampAll[[1]]$latent)[fe])
      where.fe2 <- intersect(
        grep(paste0("diff"), rownames(sampAll[[1]]$latent)[fe]), 
        intersect(
          grep(paste0(ages[i]), rownames(sampAll[[1]]$latent)[fe]), 
          grep(paste0(urs[i]), rownames(sampAll[[1]]$latent)[fe])
        ))
      where.fe <- c(where.fe1, where.fe2)
      
      struct.onegroup <- as.matrix(sampAll[[j]]$latent[re[where], 1])+
        sum(as.numeric(sampAll[[j]]$latent[fe[where.fe], 1]))
      
      struct.onesamp <- rbind (struct.onesamp, struct.onegroup)
      
    }
    
    struct.all[, j] <- struct.onesamp
    
  }
  
  # calculate quantiles
  quants <- data.frame(t(apply(struct.all, 1, stats::quantile, c(lower, 0.5, upper))))
  colnames(quants) <- c("lower", "median", "upper")
  
  quants$years <- label
  quants$group <- group
  quants$years.num <- suppressWarnings(as.numeric(as.character(quants$years)))
  quants$label <- 'RW'
  
  # calculate urban rural hazard odds ratio
  urban.idx <-  grep("urban",group) # rows contain urban trends
  rural.idx <-  grep("rural",group) # rows contain rural trends
  
  urban.samp <- struct.all[urban.idx,]
  rural.samp <- struct.all[rural.idx,]
  
  # urban rural specific hazard and odds
  
  hazard.urban.samp <- exp(urban.samp)
  hazard.rural.samp <- exp(rural.samp)
  
  odds.urban.samp <- expit(urban.samp)
  odds.rural.samp <- expit(rural.samp)
  
  ### form frame
  hazard.urban.frame <- data.frame(t(apply(hazard.urban.samp, 1, stats::quantile, c(lower, 0.5, upper))))
  hazard.rural.frame <- data.frame(t(apply(hazard.rural.samp, 1, stats::quantile, c(lower, 0.5, upper))))
  
  odds.urban.frame <- data.frame(t(apply(odds.urban.samp, 1, stats::quantile, c(lower, 0.5, upper))))
  odds.rural.frame <- data.frame(t(apply(odds.rural.samp, 1, stats::quantile, c(lower, 0.5, upper))))
  
  ### assign column name
  colnames(hazard.urban.frame) <- c("lower", "median", "upper")
  colnames(hazard.rural.frame) <- c("lower", "median", "upper")
  
  colnames(odds.urban.frame) <- c("lower", "median", "upper")
  colnames(odds.rural.frame) <- c("lower", "median", "upper")
  
  ### assign year
  hazard.urban.frame$years <- label[1:(length(label)/2)]
  hazard.rural.frame$years <- label[1:(length(label)/2)]
  
  odds.urban.frame$years <- label[1:(length(label)/2)]
  odds.rural.frame$years <- label[1:(length(label)/2)]
  
  
  ### assign age group
  hazard.urban.frame$group <- gsub( ":.*$", "", group)[1:(length(group)/2)]
  hazard.rural.frame$group <- gsub( ":.*$", "", group)[1:(length(group)/2)]
  
  odds.urban.frame$group <- gsub( ":.*$", "", group)[1:(length(group)/2)]
  odds.rural.frame$group <- gsub( ":.*$", "", group)[1:(length(group)/2)]
  
  ### assign year number
  hazard.urban.frame$years.num <- suppressWarnings(as.numeric(as.character(hazard.urban.frame$years)))
  hazard.rural.frame$years.num <- suppressWarnings(as.numeric(as.character(hazard.rural.frame$years)))
  
  odds.urban.frame$years.num <- suppressWarnings(as.numeric(as.character(odds.urban.frame$years)))
  odds.rural.frame$years.num <- suppressWarnings(as.numeric(as.character(odds.rural.frame$years)))
  
  ### assign label
  hazard.urban.frame$label <- 'Hazard: Urban'
  hazard.rural.frame$label <- 'Hazard: Rural'
  
  odds.urban.frame$label <- 'Odds: Urban'
  odds.rural.frame$label <- 'Odds: Rural'
  
  # hazard and odds ratio
  hazard.r.samp <- exp(urban.samp)/exp(rural.samp)
  odds.r.samp <- expit(urban.samp)/expit(rural.samp)
  
  hazard.r.frame <- data.frame(t(apply(hazard.r.samp, 1, stats::quantile, c(lower, 0.5, upper))))
  odds.r.frame <- data.frame(t(apply(odds.r.samp, 1, stats::quantile, c(lower, 0.5, upper))))
  
  colnames(hazard.r.frame) <- c("lower", "median", "upper")
  colnames(odds.r.frame) <- c("lower", "median", "upper")
  
  hazard.r.frame$years <- label[1:(length(label)/2)]
  odds.r.frame$years <- label[1:(length(label)/2)]
  
  hazard.r.frame$group <- gsub( ":.*$", "", group)[1:(length(group)/2)]
  odds.r.frame$group <- gsub( ":.*$", "", group)[1:(length(group)/2)]
  
  hazard.r.frame$years.num <- suppressWarnings(as.numeric(as.character(hazard.r.frame$years)))
  odds.r.frame$years.num <- suppressWarnings(as.numeric(as.character(odds.r.frame$years)))
  
  hazard.r.frame$label <- 'Hazard Ratio'
  odds.r.frame$label <- 'Odds Ratio'
  
  # prepare return object
  return.obj <- list (quants, hazard.r.frame, odds.r.frame,hazard.urban.frame,
                      hazard.rural.frame,odds.urban.frame,odds.rural.frame)
  names(return.obj) <- c('temporal.trends','hazard.ratio','odds.ratio',
                         'hazard.urban','hazard.rural','odds.urban','odds.rural')
  return (return.obj)
  
}


################################################################
#########  admin-2 model urban rural odds ratio
################################################################

trend.res <- getDiag_UR (inla_mod= fit.strat.admin2, CI = 0.95,sampAll=res.strat.admin2$draws)

odds.ratio <- trend.res$odds.ratio
hazard.ratio <- trend.res$hazard.ratio

odds.ratio <- odds.ratio[odds.ratio$group %in% c('0','1-11','12-23'),]

odds.ratio[odds.ratio$group=='0',]$group <- '0-1 months old'
odds.ratio[odds.ratio$group=='1-11',]$group <- '1-12 months old'
odds.ratio[odds.ratio$group=='12-23',]$group <- '>12 months old'

odds.ratio <-  odds.ratio %>% mutate(group = factor(group, levels=c('0-1 months old','1-12 months old','>12 months old')))


### visualization odds ratio
g.odds <- ggplot(odds.ratio, aes(x = years, y = median, ymin=lower, ymax=upper,fill=group)) +
  geom_line()  + geom_ribbon(fill="red", alpha = 0.3) +
  facet_wrap(~group) + ggtitle("Urban-rural odds ratio") +
  ylab('Odds ratio')+xlab('Year')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

setwd(paste0(res_dir,'/',country))

ggsave(g.odds, width=8, height = 4, file = paste0("Figures/Comparison/UR/", country, "-3UR-odds-ratio.pdf"))
ggsave(g.odds, width=8, height = 4, file = paste0("Figures/Report/", country, "-3UR-odds-ratio.pdf"))
ggsave(g.odds, width=8, height = 4, file = paste0("Figures/Report/", country, "-3UR-odds-ratio.tiff"))


if(FALSE){
g.hazard <- ggplot(hazard.ratio, aes(x = years, y = median, ymin=lower, ymax=upper,fill=group)) +
  geom_line()  + geom_ribbon(fill="red", alpha = 0.3) +
  facet_wrap(~group) + ggtitle("Urban-rural hazard ratio") +
  ylab('hazard ratio')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
}


################################################################
#########  admin-2 model urban rural specific hazard
################################################################

hazard.urban <- trend.res$hazard.urban
hazard.urban$strata <- 'urban'

hazard.rural <- trend.res$hazard.rural
hazard.rural$strata <- 'rural'

hazard.all <- rbind(hazard.urban,hazard.rural)
hazard.all$median <-hazard.all$median*1000

g.hazard <- ggplot(data=hazard.all, aes(x = years, y = median,color=group,linetype = strata)) +
  geom_line(size=1) +
  ggtitle("Urban rural hazards by age group") +
  ylab('Monthly hazard')+
  xlab('Years')+
  scale_y_continuous(trans = 'log10',limits=c(0.1,100))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=20))+
  theme(text = element_text(size=20))+
  scale_colour_discrete(labels = c( "0-1", "1-12",
                                    '12-24','24-36',
                                    '36-48','48-60'))

ggsave(g.hazard, width=12, height = 8, file = paste0("Figures/Comparison/UR/", country, "-3UR-hazard.pdf"))
ggsave(g.hazard, width=12, height = 8, file = paste0("Figures/Report/", country, "-3UR-hazard.pdf"))
ggsave(g.hazard, width=12, height = 8, file = paste0("Figures/Report/", country, "-3UR-hazard.tiff"))

options(warn = 0)
