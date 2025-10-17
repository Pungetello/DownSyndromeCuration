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

# adds a row to the target_attributes_tibble with the attributes needed
select_attributes = function(geo_id, metadata, target_attributes_tibble){
  
  # loop through all rows in metadata
  for (i in 1:nrow(metadata)){ 
    row = metadata[i, ]
    
    current_sex = NA
    current_ploidy = NA
    current_cell_type = NA
    ID = NA
    # TODO: could add more later, make code flexible for adding more
    
    
    #look for cell_type
    column_name = names(metadata)[grepl("cell type", names(metadata)) | # get column name for cell type
                                     grepl("cell_type", names(metadata))]
    if (length(column_name) == 0){
      column_name = names(metadata)[grepl("tissue", names(metadata))]
    }
    if (length(column_name) != 0) { # if column name exists
      current_cell_type = row[[column_name]] # current cell type is the value in that column
    }
    
    # look for ID
    column_name = names(metadata)[grepl("geo", names(metadata))] # get column name for ID
    if (length(column_name) != 0) { # if column name exists
      ID = row[[column_name]] # current ID is the value in that column
    }
    
    #look for sex
    if (any(str_detect(row, regex("female", ignore_case = TRUE)), na.rm = TRUE)){ # first check each cell in row for 'female' if found set current sex to female
      current_sex = "female" 
    } else if (any(str_detect(row, regex("male", ignore_case = TRUE)), na.rm = TRUE)){
      current_sex = "male"
    }
    
    #look for ploidy
    if (any(str_detect(row, regex("trisomic | trisomy | down syndrome | ts21", ignore_case = TRUE)), na.rm = TRUE)){ #TODO: make it match the first cell
      current_ploidy = "trisomic"
    } else if (any(str_detect(row, regex("disomic | disomy | WT | normal | control | euploid", ignore_case = TRUE)), na.rm = TRUE)){
      current_ploidy = "disomic"
    }
    
    target_attributes_tibble = add_row(target_attributes_tibble, geo_accession = ID , sex = current_sex, ploidy = current_ploidy, cell_type = current_cell_type) # append to tiblle with row contained retrieved values
  }
  
  return(target_attributes_tibble)
}         
#TODO: make this more tidyverse-ish


#--------------process_metadata-------------

file_location = "Data/Metadata/"

if (!dir.exists(file_location)){
  dir.create(file_location, recursive = TRUE)
}

# create tibble to store target attributes
target_attributes_tibble = tibble ( 
  geo_accession = character(), sex = character(), ploidy = character(), cell_type = character()   #TODO: clean up somewhow idk
  
)


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


# loop through all series IDs
for (geo_id in names(platforms_list)) {
  # read metadata into a variable, drop unneeded columns, and save it
  metadata = get_metadata(geo_id) 
  filtered_metadata = drop_cols(metadata, geo_id)
  #target_attributes_tibble = select_attributes(geo_id, filtered_metadata, target_attributes_tibble)
  result = standardize_tibble(geo_id, filtered_metadata, attr_tbl)
  print(result)
}

#print(target_attributes_tibble, n = Inf)         #TODO: reformat: GSE_ID(dataset), geo_accession, Attribute, Value
#write_tsv(target_attributes_tibble, paste0(file_location, "StandardizedMetadata.tsv"))
