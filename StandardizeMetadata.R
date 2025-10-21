#-----------loading_libraries-----------
library(GEOquery)
library(tidyverse)
source("PlatformsList.R")


#----------functions-------------

# retrieves metadata using GEOquery function
get_metadata = function(series_ID) {
  metadata = getGEO(series_ID)[[1]] 
  
  metadata = as_tibble(pData(metadata))
  return (metadata)
}


# returns version of metadata having removed columns with keywords
drop_cols = function(metadata, series_ID) {
  keywords = c("contact", "library", "processing", "description", "relation",
               "platform", "instrument", "protocol", "file", "date", "row",
               "status", "characteristics", "time", "channel", "taxid") # TODO: make sure we're removing the right things
                                                                                # in some cases, characteristics should be kept. check to see if they got extracted. Add arg?
  cols_to_drop = names(metadata)[
    grepl(paste(keywords, collapse = "|"), names(metadata))
  ]
  
  # drop columns that are not needed based on above criteria
  metadata_filtered = metadata %>% 
  select(-cols_to_drop)
  
  return(metadata_filtered)
}

#################################################################################

# goes through each row in the attributes tibble and searches the input for the attributes desired, creating a column for each in the output.
standardize_tibble <- function(geo_id, input_tbl, attr_tbl) {
  result <- map_dfc(seq_len(nrow(attr_tbl)), function(i) {
    attr <- attr_tbl[i, ]
    attr_name <- attr$attr_name
    match_type <- attr$match_type          #TODO: add GeoID column
    
    # ---- Column matching ----
    if (match_type == "column") {
      pattern <- attr$col_regex
      col_match <- names(input_tbl)[str_detect(names(input_tbl), pattern)]
      
      if (length(col_match) > 0) {
        return(tibble(!!attr_name := input_tbl[[col_match[1]]]))
      } else {
        return(tibble(!!attr_name := rep(NA_character_, nrow(input_tbl))))
      }
    }
    
    # ---- Value matching ----
    if (match_type == "value") {
      value_dict <- attr$value_dict[[1]] # because it's stored in list-column
      
      # For each column, count how many regex patterns it matches
      score_column <- function(col) {
        sum(map_int(value_dict, function(pat) {
          any(str_detect(tolower(as.character(col)), pat), na.rm = TRUE)
        }))
      }
      
      best_col <- names(input_tbl) %>%
        map_chr(~ .x) %>%
        keep(~ is.character(input_tbl[[.x]]) | is.factor(input_tbl[[.x]])) %>%
        map_int(~ score_column(input_tbl[[.x]])) %>%
        { names(input_tbl)[which.max(.)] }
      
      if (length(best_col) == 0 || score_column(input_tbl[[best_col]]) == 0) {
        return(tibble(!!attr_name := rep(NA_character_, nrow(input_tbl))))
      }
      
      col_values <- tolower(as.character(input_tbl[[best_col]]))
      standardized <- rep(NA_character_, length(col_values))
      
      for (key in names(value_dict)) {
        standardized[str_detect(col_values, value_dict[[key]])] <- key
      }
      
      return(tibble(!!attr_name := standardized))
    }
    
    stop("Unknown match type: ", match_type)
  })
  
  bind_cols(result)
  return(result)
}


##############################################################################
#TODO: understand what bot code is doing


#--------------process_metadata-------------

file_location = "Data/Metadata/"

if (!dir.exists(file_location)){
  dir.create(file_location, recursive = TRUE)
}


attr_tbl <- tibble(
  attr_name = c("ID", "Ploidy", "Sex"),
  match_type = c("column", "value", "value"),
  col_regex = c("geo", NA, NA),
  value_dict = list(
    NULL,
    list(
      trisomic = regex("trisomic|trisomy|down.syndrome|ts21", ignore_case = TRUE),
      disomic = regex("disomic|disomy|WT|normal|control|euploid", ignore_case = TRUE)
    ),
    #TODO: make a better way to construct this
    list(
      male = regex("male", ignore_case = TRUE),
      female = regex("female", ignore_case = TRUE)
    )
  )
)

output = tibble(geo_id, attr_tbl$attr_name)    #TODO: figure out how to initialize this so it will bind_rows correctly.
print(output)

# loop through all series IDs
for (geo_id in names(platforms_list)) {
  # read metadata into a variable, drop unneeded columns, get standardized metadata
  metadata = get_metadata(geo_id) 
  filtered_metadata = drop_cols(metadata, geo_id)
  result = standardize_tibble(geo_id, filtered_metadata, attr_tbl)
  print(result)
  bind_rows(output, result)
}

print(output)
#print(target_attributes_tibble, n = Inf)         #TODO: reformat: GSE_ID(dataset), geo_accession, Attribute, Value
write_tsv(output, paste0(file_location, "StandardizedMetadata.tsv"))
