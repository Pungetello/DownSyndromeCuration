#-----------loading_libraries-----------
library(GEOquery)
library(affy)
library(tidyverse)
library(BiocManager)
library(SCAN.UPC)
library(arrayQualityMetrics)
library(osfr)

source("PlatformsList.R")

#-----------functions------------------

BrainInstall2 = function(cel_files){ 
  print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
  # Specifying only one file needed
  celFilePath = cel_files[1]
  
  #Explicitly naming the source, organism, and version
  annotationSource = "entrezg"
  organism = "hs"
  version = "25.0.0"
  platform = cleancdfname(read.celfile.header(celFilePath, 
                                              info = "full")$cdfName)
  platform = sub("cdf", "", platform)
  platform = sub("stv1", "st", platform)
  platform = sub("stv2", "st", platform)
  packageName = paste(platform, organism, annotationSource, 
                      "probe", sep = "")
  
  #creates a temporary directory to hold the packages
  tmpDir = tempdir()
  packageFileName = paste(packageName, "_", version, ".tar.gz",
                          sep = "")
  
  print(packageFileName)
  tempPackageFilePath = paste(tmpDir, packageFileName, sep = "")
  packageUrl = paste("http://mbni.org/customcdf/", version,
                     "/", annotationSource, ".download/", packageFileName,
                     sep = "")
  
  
  download.file(packageUrl, tempPackageFilePath)
  
  #installs packages while ignoring errors
  try(install.packages(tempPackageFilePath, repos = NULL, type = "source"), silent = TRUE)
  print("TEST")
  return(packageName)
}


BrainInstallFromOSF = function(packageName){
  version = "25.0.0"
  
  #creates a temporary directory to hold the packages
  tmpDir = tempdir()
  packageFileName = paste(packageName, "_", version, ".tar.gz",
                          sep = "")
  
  tempPackagePath = file.path(tmpDir)
  
  files = osf_retrieve_node("b7r3g") %>%
    osf_ls_files(n_max = Inf)
  
  packageID = filter(files, name == packageFileName )%>%
    pull(id)
  
  print(tempPackagePath)
  
  osf_retrieve_file(packageID) %>%
    osf_download(tempPackagePath, conflicts = "overwrite")
  
  tempPackageFilePath = paste0(tempPackagePath, "/", libpackageFileName)
  
  install.packages(tempPackageFilePath, repos = NULL, type = "source")
}

  
#-----------installing_array_packages-----------


for (geo_id in create_unique_vector(platforms_list_new)){ #TODO: refactor so it only downloads the files list once
  
  pkgName = platforms_list_new[[geo_id]]
  
  BrainInstallFromOSF(pkgName)
}