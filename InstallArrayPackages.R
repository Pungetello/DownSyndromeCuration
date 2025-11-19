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
  
  print(packageFileName)
  
  tempPackageFilePath = paste(tmpDir, packageFileName, sep = "")
  
  files = osf_retrieve_node("b7r3g") %>%
    osf_ls_files()
  
  print(files)
  
  packageID = filter(files, names == packageFileName )%>%
    pull(id)
  
  osf_retrieve_file(packageID) %>%
    osf_download(path = tempPackageFilePath, conflicts = "overwrite")
}

  
#-----------installing_array_packages-----------

print(platform_list)

for (geo_id in create_unique_vector(platform_list)){
  
  pkgName = platform_list[[geo_id]]
  print(pkgName)
  
  BrainInstallFromOSF(pkgName)
}