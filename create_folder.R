################################################################
#########   User Input 
################################################################
rm(list = ls()) # clear the R environment and prepare for the pipeline

### Please type the name of the country you would like to analyze ### 
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.

country <- 'Senegal'

################################################################
#########   Main Code Body (NO USER INPUT NEEDED)
################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set the directory, which is the user-specified folder


# create a folder called Results to store all the results for each country, including the fitted R models, figures and tables in .csv etc. No actions if the folder is already there.
if(!dir.exists(paths = paste0('Results'))){ 
  dir.create(path = paste0('Results'))
}

# create a folder called Data to store all required data for each country including DHS survey data, country shape files etc. No actions if the folder is already there.
if(!dir.exists(paths = paste0('Data'))){ 
  dir.create(path = paste0('Data'))
}

# create a folder called Info to store R objects for each country containing all required information to run the model, 
# including the name of the country being analyzed, the years of the DHS survey data etc. No actions if the folder is already there.
if(!dir.exists(paths = paste0('Info'))){ 
  dir.create(path = paste0('Info'))
}

# create a text file to store all the countries analyzed for the first run. If created, this file is loaded and the name of the countries analyzed are retrieved.
if(!file.exists("countries_implemented.txt")){ 
  countries <- c()
}else{
  countries <- scan("countries_implemented.txt", character(), quote = "")
}


# create a folder for the user-specified country in the fold of Data, which contains the following sub-folders:
# dhsFlat: This folder contains all the GPS locations where the DHS survey is carried out.
# dhsStata: This folder contains all the survey data.
# Population: This folder contains the under-five population raster objects (i.e. population surface) of each year at the resolution of 100m and 1km.
# shapeFiles_gadm: This folder contains all the shape files of the country
# worldpop: This folder contains all the under-five population fraction of each admin2 region for each year as R objects.
if(!dir.exists(paths = paste0('Data/',country))){
  dir.create(path = paste0('Data/',country))
  if(!dir.exists(paths = paste0('Data/',country, '/dhsFlat'))){
    dir.create(path = paste0('Data/',country, '/dhsFlat'))
  }
  if(!dir.exists(paths = paste0('Data/',country, '/dhsStata'))){
    dir.create(path = paste0('Data/',country, '/dhsStata'))
  }
  if(!dir.exists(paths = paste0('Data/',country, '/Population'))){
    dir.create(path = paste0('Data/',country, '/Population'))
  }
  if(!dir.exists(paths = paste0('Data/',country, '/shapeFiles_gadm'))){
    dir.create(path = paste0('Data/',country, '/shapeFiles_gadm'))
  }
  if(!dir.exists(paths = paste0('Data/',country, '/worldpop'))){
    dir.create(path = paste0('Data/',country, '/worldpop'))
  }
  
}

# create a sub-folder for the country in the folder of Results to store the results.
if(!dir.exists(paths = paste0('Results/',country))){
  dir.create(path = paste0('Results/',country))
}

# create a sub-folder for the country in the folder of Info to store the country info.
if(!dir.exists(paths = paste0('Info/',country))){
  dir.create(path = paste0('Info/',country))
}
countries <- c(countries, country)
write(countries, file = "countries_implemented.txt") # save the newly added country into the text file.


