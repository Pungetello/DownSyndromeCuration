library(tidyverse)
# TODO: figure out for the three data sets how many metadata attributes are useful
# I'm not really sure which method is best. creating each row at once in the code, or letting it be a file that we read in?

#initialize tibble
attr_tbl = tibble(attr_name = character(), match_type = character(), col_regex = list(regex(character())), value_dict = vector("list", 0))

#add attribute rows one at a time
attr_tbl = add_row(attr_tbl, attr_name = "ID", match_type = "column", col_regex = list(regex("geo")), value_dict = NULL)

new_value_dict = list(
  affected_group = regex("trisomic|trisomy|down.syndrome|ts21|TcMAC21|Ts65Dn|Ts65Dn|1TybEmcf|treatment: untreated", ignore_case = TRUE),
  control_group = regex("disomic|disomy|WT|normal|control|euploid|wild-type|Ring chromosome 21|treatment: beta-estradiol|healthy donor (HD)", ignore_case = TRUE),
  na = regex("treatment: beta-estradiol", ignore_case = TRUE)
)
attr_tbl = add_row(attr_tbl, attr_name = "Ploidy", match_type = "value", col_regex = list(regex("")), value_dict = list(new_value_dict))

new_value_dict = list(
  female = regex("female", ignore_case = TRUE),
  male = regex("male", ignore_case = TRUE)
)
# #Sex
# attr_tbl = add_row(attr_tbl, attr_name = "Sex", match_type = "value", col_regex = list(regex("")), value_dict = list(new_value_dict))
# 
# #Tissue
# attr_tbl = add_row(attr_tbl, attr_name = "Tissue", match_type = "column", col_regex = list(regex("tissue|cell.type|source", ignore_case = TRUE)), value_dict = NULL)
# print(attr_tbl)

#========================================Dataset Metadata======================================================

dataset_attr_tbl = tibble(attr_name = character(), match_type = character(), col_regex = list(regex(character())), value_dict = vector("list", 0))

