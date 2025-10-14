#-----------loading_libraries-----------
library(GEOquery)
library(affy)
library(tidyverse)
library(BiocManager)
library(SCAN.UPC)
library(arrayQualityMetrics)

source("PlatformsList.R")
  
#-----------installing_array_packages-----------
#platform_to_package_list = list()


for (geo_id in create_unique_vector(platform_list)){
  geo_id_dir = sprintf("%s/BrainArrayExamples/%s", getwd(), geo_id)  #directory where the tar file is
  
  tar_file = sprintf("%s/%s_RAW.tar", geo_id_dir, geo_id)
  untar(tar_file, exdir = geo_id_dir)                           #gets all the cel files from the tar file
  
  platform = platform_list[[geo_id]]
  cel_files = list.files(path = geo_id_dir, pattern="^[^.]*\\.CEL\\.gz$", full.names= TRUE, ignore.case = TRUE) # TODO: find a way to just untar one cel file (I looked it up and it's actually very complicated, so maybe not...)
}
 

 
BrainInstall2 = function(cel_files){ 
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
  tempPackageFilePath = paste(tmpDir, packageFileName, sep = "")
  packageUrl = paste("http://mbni.org/customcdf/", version,
                     "/", annotationSource, ".download/", packageFileName,
                     sep = "")
  download.file(packageUrl, tempPackageFilePath)
  
  #installs packages while ignoring errors
  try(install.packages(tempPackageFilePath, repos = NULL, type = "source"), silent = TRUE)
  print("TEST")
}



BrainInstall2(cel_files)
#delete BrainArrayExamples directory and contents
unlink(sprintf("%s/BrainArrayExamples", getwd()), recursive = TRUE)

  #platform_to_package_list[[platform]] = pkgName
