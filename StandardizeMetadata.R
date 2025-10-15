#-----------loading_libraries-----------
library(GEOquery)
library(tidyverse)
source("PlatformsList.R")


#----------functions-------------

get_metadata = function(series_ID) {
  metadata = getGEO(series_ID)[[1]] # retrieves metadata using GEOquery function
  
  metadata = as_tibble(pData(metadata))
  return (metadata)
}


drop_cols = function(metadata, file_location, series_ID) {
  cols_to_drop = names(metadata)[grepl("contact", names(metadata)) |     #TODO: remove repetition
                                    grepl("library", names(metadata)) | 
                                    grepl("processing", names(metadata)) | 
                                    grepl("description", names(metadata)) |
                                    grepl("relation", names(metadata)) | 
                                    grepl("platform", names(metadata)) | 
                                    grepl("instrument", names(metadata)) |     #TODO: clarify that we are removing the correct ones
                                    grepl("protocol", names(metadata)) | 
                                    grepl("file", names(metadata)) | 
                                    grepl("date", names(metadata)) | 
                                    grepl("row", names(metadata)) | 
                                    grepl("status", names(metadata)) | 
                                    grepl("characteristics", names(metadata)) |
                                    grepl("time", names(metadata)) | 
                                    grepl("channel", names(metadata)) | 
                                    grepl("taxid", names(metadata))
  ]
  
  metadata_filtered = metadata %>% # drop columns that are not needed based on above criteria
  select(-cols_to_drop)
  
  write_tsv(metadata_filtered, paste0(file_location, series_ID, "_metadata.tsv"))
}


#--------------process_metadata-------------

file_location = "Data/Metadata/"

if (!dir.exists(file_location)){
  dir.create(file_location, recursive = TRUE)
}

for (geo_id in names(platforms_list)) { # loop through all series IDs
  metadata = get_metadata(geo_id) # get metadata and drop columns
  drop_cols(metadata, file_location, geo_id)
}

