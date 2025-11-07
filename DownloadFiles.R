#-----------loading_libraries-----------
library(GEOquery)
source("PlatformsList.R")

#-----------functions-----------
DownloadData = function(accession_id_list, platform_folder){
  # Create the working directory to your desired path
  if (!dir.exists(platform_folder)){
    dir.create(platform_folder, recursive = TRUE)
  }
  # Specify the GSE accession number and then download to a new subdirectory named after accession number
  for (geo_id in accession_id_list) {
    gse_accession = geo_id
    getGEOSuppFiles(gse_accession, makeDirectory = TRUE, baseDir = platform_folder)
    
    geo_id_dir = sprintf("%s/%s/%s", getwd(), platform_folder, geo_id)  #directory where the tar file is
    
    tar_file = sprintf("%s/%s_RAW.tar", geo_id_dir, geo_id)
    untar(tar_file, exdir = geo_id_dir)
  }
  #TODO: make delete tar file
  
  
}
#TODO: refactor to do one at a time for space efficiency?