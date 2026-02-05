#-----------loading_libraries-----------
library(tidyverse)
library(curl)

source("PlatformsList.R")

#-----------functions------------------


#Install the packeges needed by name directly from the website
SimpleBrainArrayInstall = function(packages){
  
  #create temporary directory
  tmpDir = tempdir()
  
  #explicitly name version and source
  version = "25.0.0"
  annotationSource = "entrezg"
  
  for(packageName in packages){
    
    #build the correct url
    packageFileName = paste(packageName, "_", version, ".tar.gz", sep = "")
    packageUrl = paste("http://brainarray.mbni.med.umich.edu/customcdf/", version,
                       "/", annotationSource, ".download/", packageFileName,
                       sep = "")
    
    #build temporary destination path
    tempPackageFilePath = paste(tmpDir, packageFileName, sep = "")
    
    download.file(packageUrl, tempPackageFilePath)
    
    #install packages while ignoring errors
    try(install.packages(tempPackageFilePath, repos = NULL, type = "source"), silent = TRUE)
  }
}


  
#-----------installing_array_packages-----------

#get list of each package needed with no repeats
packages = unlist(platforms_list) %>%
  unique()
packages = packages[!is.na(packages)]

#install from brain array website
SimpleBrainArrayInstall(packages)

