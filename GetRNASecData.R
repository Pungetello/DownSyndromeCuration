#change library path for supercomputer
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "" || file.access(user_lib, 2) != 0) {
  user_lib <- "~/R/library"
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
}


#-----------loading_libraries-----------
library(tidyverse)
library(RCurl)
library(biomaRt)

source("PlatformsList.R")

#-----------functions------------------

download_quality_output = function(geo_id){
  
  # DOWNLOADS FOR 101942,190053
  # NO NCBI GENERATED DATA FOR 109293,109294,184771,202938,210117
  
  link = paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=", geo_id, "&format=file&file=", geo_id, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  
  destination = paste0(getwd(), "/Data/NormalizedData/", geo_id, ".tsv.gz")
  
  if(!file.exists(destination)){
    if (RCurl::url.exists(link)){
      print(paste0("DATA EXISTS FOR ", geo_id))
      #check if link exists
      options(timeout = Inf)
      download.file(link, destination)
    } else {
      print(paste0("DATA FOR ", geo_id, " CANNOT BE DOWNLOADED"))
    }
  }
}


download_gene_data = function(){ #maybe add some parameters
  #download the human one  
  link = "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
  destination = file.path(getwd(), "Data/rna_gene_data.tsv.gz")
  if(!file.exists(destination)){
    #check if link exists
    options(timeout = Inf)
    download.file(link, destination)
  }
  return(destination)
}



#-----------Get_RNASec_Data-----------

if (!dir.exists("Data/NormalizedData")){
  dir.create("Data/NormalizedData", recursive = TRUE)
}

# download data for all RNAsec in platforms list
for (geo_id in names(platforms_list)){
  platform = platforms_list[[geo_id]]
  if (is.na(platform)){
    download_quality_output(geo_id)
  }
}

# download the genes needed
download_gene_data()
