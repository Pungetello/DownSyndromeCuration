#-----------loading_libraries-----------
library(GEOquery)
source("PlatformsList.R")

#-----------functions-----------
DownloadData = function(accession_id_list, platform_folder){
  # Create and set the working directory to your desired path
  print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa")
  print(accession_id_list)
  if (!dir.exists(platform_folder)){
    dir.create(platform_folder, recursive = TRUE)
  }
  # Specify the GSE accession number and then download to a new subdirectory named after accession number
  for (accession_id in accession_id_list) {
    gse_accession = accession_id
    getGEOSuppFiles(gse_accession, makeDirectory = TRUE, baseDir = platform_folder)
  }
}
#TODO: refactor to do one at a time for space efficiency?