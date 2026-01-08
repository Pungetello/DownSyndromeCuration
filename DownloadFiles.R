#-----------loading_libraries-----------
library(GEOquery)
source("PlatformsList.R")

#-----------functions-----------
DownloadData = function(platforms_list, platform_folder){
  accession_id_list = names(platforms_list)
  # Create the working directory to your desired path
  if (!dir.exists(platform_folder)){
    dir.create(platform_folder, recursive = TRUE)
  }
  # Specify the GSE accession number and then download to a new subdirectory named after accession number
  for (geo_id in accession_id_list) {
    gse_accession = geo_id
    getGEOSuppFiles(gse_accession, makeDirectory = TRUE, baseDir = platform_folder)
    
    #only untar the affymetrix ones
    if (!is.na(platforms_list[[geo_id]])){
    
      geo_id_dir = sprintf("%s/%s/%s", getwd(), platform_folder, geo_id)  #directory where the tar file is
      
      tar_file = sprintf("%s/%s_RAW.tar", geo_id_dir, geo_id)
      untar(tar_file, exdir = geo_id_dir)
      
      file.remove(tar_file)
    }
  }
  
}