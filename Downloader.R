#-----------loading_libraries-----------
library(GEOquery)


#-----------downloading_data-----------
download_data = function(accession_id_list, platform_folder){
  original_wd = getwd()
  # Create and set the working directory to your desired path
  if (!dir.exists(platform_folder)){
    dir.create(platform_folder)
  }
  setwd(platform_folder)
  # Specify the GSE accession number and then download to a new subdirectory named after accession number
  for (accession_id in accession_id_list) {
    gse_accession = accession_id
    getGEOSuppFiles(gse_accession, makeDirectory = TRUE)
  }
  setwd(original_wd)
}


#yotestingline22