#-----------loading_libraries-----------
library(GEOquery)
source("PlatformsList.R")

#-----------downloading_data-----------
accession_id_list = names(platforms_list)       #TODO: rework so it's only for unique ones
platform_folder = "BrainArrayExamples"
# Create and set the working directory to your desired path
if (!dir.exists(platform_folder)){
  dir.create(platform_folder)
}
# Specify the GSE accession number and then download to a new subdirectory named after accession number
for (accession_id in accession_id_list) {
  gse_accession = accession_id
  getGEOSuppFiles(gse_accession, makeDirectory = TRUE, baseDir = platform_folder)
}