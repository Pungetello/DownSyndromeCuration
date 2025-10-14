#-----------loading_libraries-----------
library(GEOquery)
library(tidyverse)
source("PlatformsList.R")


#----------functions-------------

get_metadata = function(series_ID, file_location) {
  metadata = getGEO(series_ID)[[1]] # retrieves metadata using GEOquery function
  
  metadata = as_tibble(pData(metadata))
  print(metadata)
  
  write_tsv(metadata, file_location) # writes metadata to a temporary file
}


drop_cols = function(temp_file_location, series_ID) {
  metadata = read_tsv(temp_file_location)                               #TODO: leave metadata as a variable instead of reading it in again?
  cols_to_drop = names(metadata)[grepl("contact", names(metadata)) |     #TODO: remove repetition
                                    grepl("library", names(metadata)) | 
                                    grepl("processing", names(metadata)) | 
                                    grepl("description", names(metadata)) |
                                    grepl("relation", names(metadata)) | 
                                    grepl("platform", names(metadata)) | 
                                    grepl("instrument", names(metadata)) | 
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
  
  # assign(paste0(series_ID, "_metadata_filtered"), metadata_filtered)
  # view(get(paste0(series_ID, "_metadata_filtered")))
  
  file.remove(temp_file_location) # remove temporary file
  
  file_location = "Data/Metadata/"
  write_tsv(metadata_filtered, paste0(file_location, series_ID, "_metadata.tsv"))
}


process_metadata = function(series_ID, temp_file_location) {
  if (file.exists(temp_file_location)) { # check if temporary file exists
    drop_cols(temp_file_location, series_ID) # if it does, drop columns
  } else {
    get_metadata(series_ID, temp_file_location) # if it doesn't, get metadata and drop columns
    drop_cols(temp_file_location, series_ID)
  }
}

#--------------Script-------------

temp_file_location = "Data/Metadata/temp_metadata_file.tsv"

for (geo_id in names(platforms_list)) { # loop through all series IDs
  process_metadata(geo_id, temp_file_location)
}

