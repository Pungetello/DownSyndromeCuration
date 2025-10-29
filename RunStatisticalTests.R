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
  # count occurences of each 
  counts = sapply(unique(outlierNames), function(x) {
    tests_failed = length(which(outlierNames == x))
  })
  
  
    
  print("REJECTED ASSAYS:\t") #debug
  print(paste(names(which(counts == 3)), collapse="\t")) #debug
  
  # Keep files that do not fail all three tests
  accepted_indices = which(counts != 3)
  files_to_keep = c(file_list[accepted_indices])

  file_vector = unlist(file_list)
  gsm_ids = str_extract(file_vector, "GSM\\d+")
  pass_reject_tibble = as_tibble(gsm_ids)
  pass_reject_tibble = rename(pass_reject_tibble, GSM_id = value)
  num_samples = length(gsm_ids)

  current_platforms = rep(platform, num_samples)
  #current_platforms = gsub("\t", " ", current_platforms) 
  pass_reject_tibble = add_column(pass_reject_tibble, platform = current_platforms)
  
  current_geo_id = rep(geo_id, num_samples)
  pass_reject_tibble = add_column(pass_reject_tibble, geo_id = current_geo_id)

  rejected_gsms = names(which(counts == 3))
  rejected_gsms <- str_extract(rejected_gsms, "GSM\\d+")
  pass_status <- !(gsm_ids %in% rejected_gsms)
  pass_reject_tibble <- add_column(pass_reject_tibble, pass_status)
  print(pass_reject_tibble) # debug
  
  write_tsv(pass_reject_tibble, "Data/quality_output_file.tsv", append = TRUE, col_names = FALSE)
  
  return(files_to_keep)
}


# TODO: go over this function
save_normalized_file <- function(geo_id, platform, normalized){
  normalized_tibble = as_tibble(normalized)
  normalized_tibble = normalized_tibble %>%
    rename_with(
      ~make.unique(coalesce(str_extract(., "\\d+"), .), sep = "_"),
      everything()
    )
  
  test_dataframe = as.data.frame(normalized)
  normalized_row_names = rownames(test_dataframe)
  new_names = str_extract(normalized_row_names, "GSM\\d+")
  
  final_tibble = normalized_tibble %>%
    add_column("Sample_ID" = new_names, .before = 1)
  
  tibble_file_location = paste0("Data/", geo_id, platform, ".tsv.gz")
  write_tsv(final_tibble, tibble_file_location)
}

# TODO: go over this function
get_scan_upc_files <- function(cel_files_id, platform_to_package_list, platform){
  
  # Sets the file pattern to .CEL, so scan pulls everything with that ending
  celFilePattern <- file.path(tar_file_output_f, "*.CEL.gz")
  
  # formated string for the SCAN output
  scan_output_file_f = sprintf("affymetrix_data/%s_SCAN", geo_id)
  
  # This cleans up the data and removes outliers
  platform = unlist(platform)
  pkgName = platform_to_package_list[[platform]]
  
  # last step to converting the information
  print('test2')
  normalized = SCAN(celFilePattern, convThreshold = .9, probeLevelOutDirPath = NA, probeSummaryPackage=pkgName)
  print('test3')
  print(normalized)
  print(typeof(normalized))
  return (normalized)
}


#-----------run statistical tests-----------

#TODO: add platforms to package list for normalization: waiting on website


#for (geo_id in names(platforms_list)){
geo_id = "GSE1789"
  file_data = sprintf("%s/Data/Files/%s", getwd(), geo_id)
  file_list = list.files(path = file_data, pattern="^[^.]*\\.CEL\\.gz$", full.names= TRUE, ignore.case = TRUE)
  filtered_file_list = quality_control_removal(file_list, platforms_list[[geo_id]], geo_id)
  
  #normalized = get_scan_upc_files(file_list, platform_to_package_list, platform)
  #save_normalized_file(geo_id, platform, normalized)
  
#}

