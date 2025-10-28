#-----------loading_libraries-----------
library(GEOquery)
library(tidyverse)
source("PlatformsList.R")
source("MetadataAttributes.R")


#----------functions-------------
#TODO: expand characteristics column if we see it's needed
# retrieves metadata using GEOquery function
get_metadata = function(geo_ID) {
  metadata = getGEO(geo_ID)[[1]] 
  
  metadata = as_tibble(pData(metadata))
  return (metadata)
}


# returns version of metadata having removed columns with keywords
drop_cols = function(metadata, series_ID) {
  keywords = c("contact", "library", "processing", "description", "relation",
               "platform", "instrument", "protocol", "file", "date", "row",
               "status", "characteristics", "time", "channel", "taxid") # TODO: figure out what to do for characteristics column.
                                                                        # in some cases, characteristics should be kept. check to see if they got extracted. Add arg?
  cols_to_drop = names(metadata)[
    grepl(paste(keywords, collapse = "|"), names(metadata))
  ]
  
  # drop columns that are not needed based on above criteria
  metadata_filtered = metadata %>% 
  select(-cols_to_drop)
  
  return(metadata_filtered)
}


# goes through each row in the attributes tibble and searches the input for the attributes desired, creating a column for each in the output
standardize_tibble = function(geo_id, input_tbl, attr_tbl) {
  n = nrow(input_tbl)
  geo_col = rep(as.character(geo_id), n)
  
  # For each column, count how many regex pattern matches
  score_column = function(col, value_dict) {
    sum(map_int(value_dict, function(regex_pattern) {
      any(str_detect(col, regex_pattern), na.rm = TRUE)
    }))
  }
  
  # Iterate over rows of attr_tbl and return one tibble column per attribute
  cols_list = map(seq_len(nrow(attr_tbl)), function(i) { # 1:nrow(attr_tbl)
    attr = attr_tbl[i, ]
    attr_name = as.character(attr$attr_name)
    match_type = as.character(attr$match_type)
    
    if (match_type == "column") {
      pattern = as.character(attr$col_regex)
      col_match = names(input_tbl)[str_detect(names(input_tbl), pattern)]
      if (length(col_match) > 0) {
        # return the first column it matches the title of
        return(tibble(!!attr_name := input_tbl[[col_match[1]]]))
      } else {
        # error handling: if no column title matches, return column of all NA
        return(tibble(!!attr_name := rep(NA_character_, n)))
      }
    }
    
    else if (match_type == "value") {
      value_dict = attr$value_dict[[1]]
      
      # Keep only character or factor columns
      kept_names = names(input_tbl)[map_lgl(names(input_tbl), function(.x){is.character(input_tbl[[.x]]) || is.factor(input_tbl[[.x]])})]
      if (length(kept_names) == 0) {
        return(tibble(!!attr_name := rep(NA_character_, n)))
      }
      
      # create named list with each column stored as a vector
      kept_cols = map(kept_names, function(.x){input_tbl[[.x]]})
      names(kept_cols) = kept_names
      
      # Score each kept column
      scores = map_int(kept_cols, function(.x){score_column(.x, value_dict)})
      
      # If no column matched any pattern, return NA column
      if (all(scores == 0)) {
        return(tibble(!!attr_name := rep(NA_character_, n)))
      }
      
      # Choose best column (first with max score)
      best_index = which.max(scores)
      best_col_name = kept_names[best_index]
      best_col_values = kept_cols[[best_index]]
      
      # Build standardized vector (first-match-wins)
      standardized = rep(NA_character_, length(best_col_values))
      for (key in names(value_dict)) {
        pattern = value_dict[[key]]
        matches = str_detect(best_col_values, pattern)
        # Assign key only where standardized is still NA (first match wins)
        standardized[is.na(standardized) & matches] = key
      }
      
      return(tibble(!!attr_name := standardized))
    }
    else{
      stop("Unknown match type: ", match_type)
    }
  })
  
  # Combine all columns into one tibble
  result = bind_cols(tibble(GeoID = geo_col), bind_cols(cols_list))
  return(result)
}
#TODO: make bespoke solutions to file idiosyncrasies

#--------------process_metadata-------------

file_location = "Data/Metadata/"

if (!dir.exists(file_location)){
  dir.create(file_location, recursive = TRUE)
}


# loop through all series IDs
combined_output = tibble()
for (geo_id in names(platforms_list)) {
  # read metadata into a variable, drop unneeded columns, get standardized metadata
  metadata = get_metadata(geo_id) 
  filtered_metadata = drop_cols(metadata, geo_id)
  result = standardize_tibble(geo_id, filtered_metadata, attr_tbl)
  combined_output = bind_rows(combined_output, result)
}


rotated_standardized_metadata = pivot_longer(combined_output, !c("GeoID", "ID"), names_to = "Attribute", values_to = "Value")

write_tsv(rotated_standardized_metadata, paste0(file_location, "StandardizedMetadata.tsv"))
