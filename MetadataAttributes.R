library(tidyverse)

# I'm not really sure which method is best. creating each row at once in the code, or letting it be a file that we read in?

#initialize tibble
attr_tbl = tibble(attr_name = character(), match_type = character(), col_regex = regex(character()), value_dict = vector("list", 0))

#add attribute rows one at a time
attr_tbl = add_row(attr_tbl, attr_name = "ID", match_type = "column", col_regex = regex("geo"), value_dict = NULL)

new_value_dict = list(
  trisomic = regex("trisomic|trisomy|down.syndrome|ts21", ignore_case = TRUE),
  disomic = regex("disomic|disomy|WT|normal|control|euploid", ignore_case = TRUE)
)
attr_tbl = add_row(attr_tbl, attr_name = "Ploidy", match_type = "value", col_regex = regex(""), value_dict = list(new_value_dict))

new_value_dict = list(
  male = regex("male", ignore_case = TRUE),
  female = regex("female", ignore_case = TRUE)
)
attr_tbl = add_row(attr_tbl, attr_name = "Sex", match_type = "value", col_regex = regex(""), value_dict = list(new_value_dict))

# rotated_standardized_metadata = rotate_longer(combined_output, !c("GeoID", "ID"), names_to = "Attribute", values_to = "Value")