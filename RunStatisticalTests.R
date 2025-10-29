#-----------loading_libraries-----------
source("PlatformsList.R")
library(arrayQualityMetrics)
library(SCAN.UPC)
library(tidyverse)


#-----------functions------------------

#line 286 uses SCAN with the name of the Brainarray package to normalize the data

quality_control_removal = function(file_list, platform, geo_id){

  cel_files = read.celfiles(file_list) # TODO: perhaps refactor so it reads in one cel file at a time?
  
  test_results = arrayQualityMetrics(expressionset = cel_files, force = TRUE, outdir = "QualityOutput")
  
  unlink("QualityOutput", recursive = TRUE) #TODO: figure out how to not download in first place
  
  
  outlierIndices = c(test_results$modules$heatmap@outliers@which, test_results$modules$boxplot@outliers@which, test_results$modules$maplot@outliers@which)
  
  # Get the assay names corresponding to the indices above from the arrayTable data frame.
  outlierNames = test_results$arrayTable$sampleNames[outlierIndices]
  
  # Now we have to count occurences of each assay name in the outlierNames
  # vector. Any that are in three times and were thus classified as outliers
  # by all three methods should be rejected.
  # This part uses sapply() to go through the unique assay names found in the
  # outlierNames vector, and use length() and which() to count how many times each
  # one appears. #TODO: make this comment shorter
  counts = sapply(unique(outlierNames), function(x) {
    # len is the number of times the assay name (x) appears in the outlierNames
    # vector.
    len = length(which(outlierNames == x))
  })
  
  
    
  print("REJECTED ASSAYS:\t") #debug
  print(paste(names(which(counts == 3)), collapse="\t")) #debug
  
  # Keep files that do not fail all three tests
  accepted_indices = which(counts != 3)
  files_to_keep = c(file_list[accepted_indices])

  file_vector = unlist(file_list)
  gsm_ids = str_extract(file_vector, "GSM\\d+")
  pass_reject_tibble = as_tibble(gsm_ids)
  # pass_reject_tibble = mutate(pass_reject_tibble, cel_file = gsm_ids) #TODO: what does this do? Delete?
  num_samples = length(gsm_ids)

  current_platforms = rep(platform, num_samples)
  current_platforms = unlist(current_platforms)
  current_platforms = gsub("\t", " ", current_platforms)
  pass_reject_tibble = add_column(pass_reject_tibble, platform = current_platforms)
  
  current_geo_id <- rep(geo_id, num_samples)
  current_geo_id <- unlist(current_geo_id)
  pass_reject_tibble <- add_column(pass_reject_tibble, geo_id = current_geo_id)

  rejected_gsms <- names(which(counts == 3))
  rejected_gsms <- str_extract(rejected_gsms, "GSM\\d+")
  # print(paste(rejected_gsms, '3333'))
  pass_status <- !(gsm_ids %in% rejected_gsms)
  # print(gsm_ids)
  # print(pass_status)
  pass_reject_tibble <- add_column(pass_reject_tibble, "Pass?" = pass_status)
  
  write_tsv(pass_reject_tibble, "Data/quality_output_file.tsv", append = TRUE, col_names = FALSE)
  
  return(files_to_keep)
}

#-----------run statistical tests-----------

#TODO: add platforms to package list for normalization?


for (geo_id in names(platforms_list)){
  file_data = sprintf("%s/Data/Files/%s", getwd(), geo_id)
  file_list = list.files(path = file_data, pattern="^[^.]*\\.CEL\\.gz$", full.names= TRUE, ignore.case = TRUE)
  filtered_file_list = quality_control_removal(file_list, platforms_list[geo_id], geo_id)
  
  
  
}

