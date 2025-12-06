#-----------loading_libraries-----------
library(GEOquery)
library(affy)
library(tidyverse)
library(BiocManager)
library(SCAN.UPC)
library(arrayQualityMetrics)
library(osfr)
library(archive)

source("PlatformsList.R")

#-----------functions------------------

BrainInstallFromOSF = function(project_link, packages){
  version = "25.0.0"
  
  #creates a temporary directory to hold the packages
  tmpDir = tempdir()
  tempPackagePath = file.path(tmpDir)
  
  #gets tibble of all the files in the OSF project
  files = osf_retrieve_node(project_link) %>%
    osf_ls_files(n_max = Inf)%>%
    arrange(name) #debug
  print(files, n=Inf) #debug
  
  for (packageName in packages){
  
    packageFileName = paste(packageName, "_", version, ".tar.gz",
                            sep = "")
    print(packageFileName) #debug
    #Use the file list to find the file ID of the package
    packageID = filter(files, name == packageFileName )%>%
      pull(id)
    
    #Download and install package by its ID
    osf_retrieve_file(packageID) %>%
      osf_download(tempPackagePath, conflicts = "overwrite")
    print(list.files(tempPackagePath))
    
    #tempPackageFilePath = paste0(tempPackagePath, "/", packageFileName)
    tempPackageFilePath = file.path(tempPackagePath, packageFileName)
  }
}

BrainInstallFromZenodo = function(packages){
  version = "25.0.0"
  zipFile = paste0(getwd(), "/BrainArray/BrainArrays.zip")
  
  if(!file.exists(zipFile)){
    options(timeout = Inf)
    download.file("https://zenodo.org/records/17808641/files/BrainArray.zip?download=1", zipFile)
  }
  
  packages = unlist(lapply(packages, function(x) paste0("BrainArray/", x, "_", version, ".tar.gz", method = curl)))
  #pull out packages from zip
  
  unzip(zipFile, files=packages)
  
}

  
#-----------installing_array_packages-----------

if (!dir.exists("BrainArray")){
  dir.create("BrainArray", recursive = TRUE)
}

packages = unlist(platforms_list) %>%
  unique()
packages = packages[!is.na(packages)]

BrainInstallFromZenodo(packages)
#BrainInstallFromOSF("b7r3g", packages)
