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
  
platform_to_package_list = list()
  
for (geo_id in names(platform_list)){
  geo_id_dir = sprintf("%s/Platforms/%s", original_wd, geo_id)
  
  tar_file = sprintf("%s/%s_RAW.tar", geo_id_dir, geo_id)
  untar(tar_file, exdir = geo_id_dir)
  
  platform = platform_list[[geo_id]]
  cel_files = list.files(path = geo_id_dir, pattern="^[^.]*\.CEL\.gz$", full.names= TRUE, ignore.case = TRUE)
  print(cell_files[1])
  pkgName = InstallBrainArrayPackage(cel_files[1], "25.0.0", "hs", "entrezg")           # in SCAN.UPC package, finds what package needed to read the type of file in arg 1
      platform_to_package_list[[platform]] = pkgName
}
