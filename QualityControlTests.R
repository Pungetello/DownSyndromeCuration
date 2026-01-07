#-----------loading_libraries-----------
source("PlatformsList.R")
library(arrayQualityMetrics)
library(SCAN.UPC)
library(tidyverse)


#-----------functions------------------

#Runs statistical tests on the data, determining which ones should be removed due to low quality
quality_control_removal = function(file_list, platform, geo_id){
  
  cel_files = read.celfiles(file_list) # TODO: perhaps refactor so it reads in one cel file at a time?
  
  test_results = arrayQualityMetrics(expressionset = cel_files, force = TRUE, outdir = "QualityOutput")
  
  unlink("QualityOutput", recursive = TRUE)
  
  outlierIndices = c(test_results$modules$heatmap@outliers@which, test_results$modules$boxplot@outliers@which, test_results$modules$maplot@outliers@which)
  # Get the assay names corresponding to the indices above from the arrayTable data frame.
  outlierNames = test_results$arrayTable$sampleNames[outlierIndices]
  
  # count occurrences of each assay name in outlierNames: if it appears 3 times, it is an outlier and should be removed
  counts = sapply(unique(outlierNames), function(x) {
    tests_failed = length(which(outlierNames == x))
  })
  
  
  # Keep files that do not fail all three tests
  accepted_indices = which(counts != 3)
  files_to_keep = c(file_list[accepted_indices])
  
  file_vector = unlist(file_list)
  gsm_ids = str_extract(file_vector, "GSM\\d+")
  pass_reject_tibble = as_tibble(gsm_ids)
  pass_reject_tibble = dplyr::rename(pass_reject_tibble, GSM_id = value)
  num_samples = length(gsm_ids)
  
  current_platforms = rep(platform, num_samples)
  pass_reject_tibble = add_column(pass_reject_tibble, platform = current_platforms)
  
  current_geo_id = rep(geo_id, num_samples)
  pass_reject_tibble = add_column(pass_reject_tibble, geo_id = current_geo_id)
  
  rejected_gsms = names(which(counts == 3))
  rejected_gsms = str_extract(rejected_gsms, "GSM\\d+")
  pass_status = !(gsm_ids %in% rejected_gsms)
  pass_reject_tibble = add_column(pass_reject_tibble, pass_status)
  colnames(pass_reject_tibble) = c("GSM_ID", "Platform", "GSE_ID", "DataQuality")
  pass_reject_tibble = mutate(pass_reject_tibble, DataQuality = ifelse(DataQuality == "TRUE", "PASS", "FAIL"))
  
  
  if (!file.exists("Data/quality_output_file.tsv")) {
    write_tsv(pass_reject_tibble, "Data/quality_output_file.tsv")
  } else {
    write_tsv(pass_reject_tibble, "Data/quality_output_file.tsv", append = TRUE, col_names = FALSE)
  }
  
  return(files_to_keep)
}


#TODO: fix NA issue 
#Uses the SCAN function to normalize the data
get_scan_upc_files = function(cel_files_id, pkgName, geo_id){
  convThreshold = 0.9
  
  # Set the file pattern to .CEL, so scan pulls everything with that ending
  celFilePattern = file.path(sprintf("%s/Data/Files/%s", getwd(), geo_id), "*.CEL*") 
  
  # clean up the data and removes outliers
  #normalized = SCAN(geo_id, convThreshold = convThreshold) #1
  #print(normalized)
  #normalized = exprs(normalized)
  normalized = SCAN(celFilePattern, convThreshold = .9, probeLevelOutDirPath = NA, probeSummaryPackage=pkgName) #3 only one that only does some NA's instead of all
  #normalized = SCAN(celFilePattern, convThreshold = .9, probeLevelOutDirPath = NA) #should be 0.01, 0.9 for debug #4
  return (normalized)
}



#changes normalized data to matrix, creates new directory for these files and saves them gzipped
save_normalized_file = function(geo_id, platform, normalized) {
  
  normalized = as.matrix(normalized)
  output_dir = file.path(getwd(), "Data", "NormalizedData")
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
  
  file_path = file.path(output_dir, paste0(geo_id, "_", platform, ".tsv.gz"))
  if (!file.exists(file_path)){
  
    gsm_pattern = "ENSG\\d+_at"
    sample_ids = stringr::str_extract(rownames(normalized), gsm_pattern)
    
    final_tibble = as_tibble(normalized, rownames = NULL) %>%
      bind_cols(Sample_ID = sample_ids, .)
    
    write_tsv(final_tibble, file_path)
  }
}



#removes samples from the SampleMetadata file if they were rejected in the previous step
remove_fails_from_metadata = function(pass_reject_tibble){
  failed_gsm_ids = filter(pass_reject_tibble, DataQuality == "FAIL")%>%
    pull(GSM_ID)
  
  metadata_location = paste0(getwd(), "/Data/Metadata/SampleMetadata.tsv")
  metadata_tibble = read_tsv(metadata_location)%>%
    filter(! ID %in% failed_gsm_ids)
  
  write_tsv(metadata_tibble, metadata_location)
}



#-----------run statistical tests-----------

for (geo_id in names(platforms_list)){
  platform = platforms_list[[geo_id]]
  if (!is.na(platform)){
  
    file_data = sprintf("%s/Data/Files/%s", getwd(), geo_id)
    file_list = list.files(path = file_data, pattern="^[^.]*\\.CEL\\.gz$", full.names= TRUE, ignore.case = TRUE)
    filtered_file_list = quality_control_removal(file_list, platform, geo_id)
    
    #should be working, once we can run InstallArrayPackages again.
    normalized = get_scan_upc_files(file_list, platform, geo_id)
    save_normalized_file(geo_id, platform, normalized)
  }
}

quality_output_tibble = read_tsv(paste0(getwd(), "/Data/quality_output_file.tsv"))
remove_fails_from_metadata(quality_output_tibble)


#TODO: make it not make a new normalized data file if it already exists