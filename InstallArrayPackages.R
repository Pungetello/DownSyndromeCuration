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
    
    #tempPackageFilePath = paste0(tempPackagePath, "/", packageFileName)
    tempPackageFilePath = file.path(tempPackagePath, packageFileName)
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
    
    #tempPackageFilePath = paste0(tempPackagePath, "/", packageFileName)
    tempPackageFilePath = file.path(tempPackagePath, packageFileName)

    install.packages(tempPackageFilePath, repos = NULL, type = "source")
  }
}

  
#-----------installing_array_packages-----------

packages = unlist(platforms_list) %>%
  unique()
packages = packages[!is.na(packages)]

print(packages)


BrainInstallFromOSF("b7r3g", packages)
    install.packages(tempPackageFilePath, repos = NULL, type = "source")
  }
}

  
#-----------installing_array_packages-----------

packages = unlist(platforms_list) %>%
  unique()
packages = packages[!is.na(packages)]

print(packages)


BrainInstallFromOSF("b7r3g", packages)
    install.packages(tempPackageFilePath, repos = NULL, type = "source")
  }
}

  
#-----------installing_array_packages-----------

packages = unlist(platforms_list) %>%
  unique()
packages = packages[!is.na(packages)]

print(packages)


BrainInstallFromOSF("b7r3g", packages)
