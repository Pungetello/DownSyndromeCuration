
#-----------loading_libraries-----------



#----------functions-------------

make_sample_metadata = function(geo_id, sample_metadata){
  DE_filename = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_DE.tsv")
  if(!file.exists(DE_filename)){
    print("NO DifExpAnalysis RESULTS")
    return()
  }
  
  DE_file = read_tsv(DE_filename)
  
  
  
}





#--------------metadata_tables-------------

#read in sample and dataset metadata, along with RPKM files, to make the tables

sample_metadata = read_tsv(paste0(getwd(), "/Data/Metadata/SampleMetadata.tsv"))

for (geo_id in names(platforms_list)) {
  
  #make Sample_metadata
  make_sample_metadata(geo_id, sample_metadata)
  
  #make Abundance_data
  
  #make Differential_analysis_results
  
  
}

