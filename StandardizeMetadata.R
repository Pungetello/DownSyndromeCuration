#-----------loading_libraries-----------
library(GEOquery)
library(tidyverse)
library(janitor)
source("PlatformsList.R")
source("MetadataAttributes.R")


#----------functions-------------
#TODO: expand characteristics column if we see it's needed
# retrieves metadata using GEOquery function
get_metadata = function(geo_ID) {
  metadata = getGEO(geo_ID)[[1]] 
  
  metadata = as_tibble(pData(metadata))
  metadata = clean_names(metadata)
  metadata = fix_bespoke_issues(geo_ID, metadata) # debug
  return (metadata)
}


# fixes specific issues in the data for certain datasets
fix_bespoke_issues = function(geo_ID, metadata){
  if(geo_ID == "GSE1789"){      #TODO: fix it more
    metadata = mutate(metadata, age_ch1 = ifelse(is.na(age_ch1), substr(characteristics_ch1_1, start = 5, stop = nchar(characteristics_ch1_1) ), age_ch1)) %>%
      mutate(tissue_ch1 = ifelse(is.na(tissue_ch1), tissue_ch1[1], tissue_ch1))
  }else{
    return(metadata)
  }
}


# returns version of metadata having removed columns with keywords
drop_cols = function(metadata, series_ID) {
  keywords = c("contact", "library", "processing", "relation",
               "platform", "instrument", "file", "date", "row",
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


# splits the metadata into two tibbles: one containing only columns where there is a difference in values, one with only columns where they're all the same.
split_metadata = function(metadata) {
  diff_metadata = select(metadata, where(~n_distinct(.) > 1))
  print(diff_metadata)
  same_metadata = select(metadata, where(~n_distinct(.) == 1))
  print(same_metadata)
  return(c(diff_metadata, same_metadata))
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


#--------------process_metadata-------------

file_location = "Data/Metadata/"

if (!dir.exists(file_location)){
  dir.create(file_location, recursive = TRUE)
}


# loop through all series IDs
combined_output = tibble()
dataset_combined_output = tibble()

for (geo_id in names(platforms_list)) {
  
  # read metadata into a variable, drop unneeded columns, split into same and diff
  metadata = get_metadata(geo_id)
  metadata = drop_cols(metadata, geo_id)
  diff_metadata = select(metadata, where(~n_distinct(.) > 1))
  same_metadata = select(metadata, where(~n_distinct(.) == 1))
  print(same_metadata, n=Inf)
  print(diff_metadata, n = Inf)
  
  # get standardized metadata
  result = standardize_tibble(geo_id, diff_metadata, attr_tbl)
  
  # remove any attributes where every row is NA for this geoID, then rotate
  result = Filter(function(x)!all(is.na(x)), result)
  rotated_result = pivot_longer(result, !c("GeoID", "ID"), names_to = "Attribute", values_to = "Value")
  combined_output = bind_rows(combined_output, rotated_result)
  print(result)
  #get dataset metadata 
  dataset_result = standardize_tibble(geo_id, same_metadata, dataset_attr_tbl) #TODO: make dataset attribute tibble
  dataset_combined_output = bind_rows(dataset_combined_output)
}

write_tsv(combined_output, paste0(file_location, "StandardizedMetadata.tsv"))
write_tsv(dataset_combined_output, paste0(file_location, "DatasetMetadata.tsv"))
