library(GEOquery)
library(affy)
library(tidyverse)
library(BiocManager)
library(SCAN.UPC) 
#library(lubridate)
library(arrayQualityMetrics)

platform_list = list(
  "GSE110064" = "GPL570",
  "GSE48611"= "GPL570",
  "GSE1789" = "GPL96"
)

#temp_dir = tempdir()
#setwd(temp_dir)
  
platform_to_package_list = list()

#setwd("~/DownSyndromeCuration")
  
for (geo_id in names(platform_list)){
  original_wd = getwd()
  print(original_wd)
  new_wd = sprintf("%s/Platforms/%s", original_wd, geo_id)
  print(new_wd)
  setwd(new_wd)
  
  tar_file = sprintf("%s_RAW.tar", geo_id)
  cel_folder = dir.create(tempfile("cel_files"))
  untar(tar_file, exdir = cel_folder)
  
  platform = platform_list[[geo_id]]
  cel_files = list.files(path = cel_folder, pattern="^[^.]*\\.CEL\\.gz$", full.names= TRUE, ignore.case = TRUE)
  pkgName = InstallBrainArrayPackage(cel_files[1], "25.0.0", "hs", "entrezg")   # in SCAN.UPC package, finds what package needed to read the type of file in arg 1
      platform_to_package_list[[platform]] = pkgName
  
  setwd(original_wd)
}