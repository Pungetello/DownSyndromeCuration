#-----------loading_libraries-----------
library(GEOquery)
library(affy)
library(tidyverse)
library(BiocManager)
library(SCAN.UPC)
library(arrayQualityMetrics)
library(osfr)
library(archive)
library(curl)

source("PlatformsList.R")

#-----------functions------------------

BrainInstallFromZenodo = function(packages){
  version = "25.0.0"
  zipFile = paste0(getwd(), "/BrainArray/BrainArrays.zip")
  
  if (!file.exists(zipFile)) {
    options(timeout = 99999)
    curl::curl_download(
      "https://zenodo.org/records/17808641/files/BrainArray.zip?download=1",
      destfile = "BrainArray/BrainArrays.zip"
    )
  }
  
  packages = unlist(lapply(packages, function(x) paste0("BrainArray/", x, "_", version, ".tar.gz")))
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