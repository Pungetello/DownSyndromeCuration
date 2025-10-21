#-----------loading_libraries-----------
library(GEOquery)
library(tidyverse)
source("PlatformsList.R")


#----------functions-------------

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
library(tibble)
library(dplyr)
library(purrr)
library(stringr)
library(rlang)

standardize_tibble <- function(geo_id, input_tbl, attr_tbl) {
  # Basic checks
  n <- nrow(input_tbl)
  if (length(geo_id) == 1) {
    geo_col <- rep(as.character(geo_id), n)
  } else if (length(geo_id) == n) {
    geo_col <- as.character(geo_id)
  } else {
    stop("geo_id must be length 1 or length nrow(input_tbl).")
  }
  
  # Helper: score how many distinct dict-patterns appear anywhere in a column
  score_column <- function(col, value_dict) {
    # col: vector (already coerced to character/lower)
    sum(map_int(value_dict, function(pat) {
      pat_collapsed <- if (length(pat) > 1) paste(pat, collapse = "|") else pat
      any(str_detect(col, pat_collapsed), na.rm = TRUE)
    }))
  }
  
  # Iterate over rows of attr_tbl and return one tibble column per attribute
  cols_list <- map(seq_len(nrow(attr_tbl)), function(i) {
    attr <- attr_tbl[i, ]
    attr_name <- as.character(attr$attr_name)
    match_type <- as.character(attr$match_type)
    
    if (match_type == "column") {
      pattern <- as.character(attr$col_regex)
      col_match <- names(input_tbl)[str_detect(names(input_tbl), pattern)]
      if (length(col_match) > 0) {
        # Use first match; preserve original column (no forced coercion)
        return(tibble(!!attr_name := input_tbl[[col_match[1]]]))
      } else {
        return(tibble(!!attr_name := rep(NA_character_, n)))
      }
    }
    
    if (match_type == "value") {
      # value_dict is stored as a list-column (one element per attr row)
      value_dict <- attr$value_dict[[1]]
      if (is.null(value_dict) || length(value_dict) == 0) {
        return(tibble(!!attr_name := rep(NA_character_, n)))
      }
      
      # Keep only character or factor columns
      kept_names <- names(input_tbl)[map_lgl(names(input_tbl), ~ is.character(input_tbl[[.x]]) || is.factor(input_tbl[[.x]]))]
      if (length(kept_names) == 0) {
        return(tibble(!!attr_name := rep(NA_character_, n)))
      }
      
      # Precompute lowercased character versions of kept columns (for performance)
      kept_cols <- map(kept_names, ~ tolower(as.character(input_tbl[[.x]])))
      names(kept_cols) <- kept_names
      
      # Score each kept column
      scores <- map_int(kept_cols, ~ score_column(.x, value_dict))
      
      # If no column matched any pattern => return NA column
      if (all(scores == 0)) {
        return(tibble(!!attr_name := rep(NA_character_, n)))
      }
      
      # Choose best column (first with max score)
      best_idx <- which.max(scores)
      best_col_name <- kept_names[best_idx]
      col_values <- kept_cols[[best_idx]]
      
      # Build standardized vector (first-match-wins)
      standardized <- rep(NA_character_, length(col_values))
      for (key in names(value_dict)) {
        pat <- value_dict[[key]]
        pat_collapsed <- if (length(pat) > 1) paste(pat, collapse = "|") else pat
        matches <- str_detect(col_values, pat_collapsed)
        # Assign key only where standardized is still NA (first match wins)
        standardized[is.na(standardized) & matches] <- key
      }
      
      return(tibble(!!attr_name := standardized))
    }
    
    stop("Unknown match type: ", match_type)
  })
  
  # Bind all standardized attribute columns side-by-side; put GeoID first
  result <- bind_cols(tibble(GeoID = geo_col), bind_cols(cols_list))
  return(result)
}

# goes through each row in the attributes tibble and searches the input for the attributes desired, creating a column for each in the output.
# standardize_tibble <- function(geo_id, input_tbl, attr_tbl) {
#   result <- map_dfc(seq_len(nrow(attr_tbl)), function(i) {
#     attr <- attr_tbl[i, ]
#     attr_name <- attr$attr_name
#     match_type <- attr$match_type          #TODO: add GeoID column
#     
#     # ---- Column matching ----
#     if (match_type == "column") {
#       pattern <- attr$col_regex
#       col_match <- names(input_tbl)[str_detect(names(input_tbl), pattern)]
#       
#       if (length(col_match) > 0) {
#         return(tibble(!!attr_name := input_tbl[[col_match[1]]]))
#       } else {
#         return(tibble(!!attr_name := rep(NA_character_, nrow(input_tbl))))
#       }
#     }
#     
#     # ---- Value matching ----
#     if (match_type == "value") {
#       value_dict <- attr$value_dict[[1]] # because it's stored in list-column
#       
#       # For each column, count how many regex patterns it matches
#       score_column <- function(col) {
#         sum(map_int(value_dict, function(pat) {
#           any(str_detect(tolower(as.character(col)), pat), na.rm = TRUE)
#         }))
#       }
#       
#       best_col <- names(input_tbl) %>%
#         map_chr(~ .x) %>%
#         keep(~ is.character(input_tbl[[.x]]) | is.factor(input_tbl[[.x]])) %>%
#         map_int(~ score_column(input_tbl[[.x]])) %>%
#         { names(input_tbl)[which.max(.)] }
#       
#       if (length(best_col) == 0 || score_column(input_tbl[[best_col]]) == 0) {
#         return(tibble(!!attr_name := rep(NA_character_, nrow(input_tbl))))
#       }
#       
#       col_values <- tolower(as.character(input_tbl[[best_col]]))
#       standardized <- rep(NA_character_, length(col_values))
#       
#       for (key in names(value_dict)) {
#         standardized[str_detect(col_values, value_dict[[key]])] <- key
#       }
#       
#       return(tibble(!!attr_name := standardized))
#     }
#     
#     stop("Unknown match type: ", match_type)
#   })
#   
#   bind_cols(result)
#   return(result)
# }
# 

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



# loop through all series IDs
combined_output = tibble()
for (geo_id in names(platforms_list)) {
  # read metadata into a variable, drop unneeded columns, get standardized metadata
  metadata = get_metadata(geo_id) 
  filtered_metadata = drop_cols(metadata, geo_id)
  result = standardize_tibble(geo_id, filtered_metadata, attr_tbl)
  combined_output = bind_rows(combined_output, result)
}
print(combined_output, n=Inf) #debug

#print(target_attributes_tibble, n = Inf)         #TODO: reformat: GSE_ID(dataset), geo_accession, Attribute, Value
write_tsv(output, paste0(file_location, "StandardizedMetadata.tsv"))
