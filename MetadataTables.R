
#-----------loading_libraries-----------

library(tidyverse)
library(readxl)
library(janitor)
library(DBI)

#----------functions-------------

# retrieves metadata using GEOquery function
get_metadata = function(geo_ID) {
  metadata = getGEO(geo_ID)[[1]] 
  
  metadata = as_tibble(pData(metadata))
  metadata = clean_names(metadata)
  metadata = fix_bespoke_issues(geo_ID, metadata) # debug
  return (metadata)
}
#TODO: these kinda repeat each other, refactor
get_gse_metadata = function(gse_id) {
  # Get GEO dataset (GSEMatrix = FALSE returns the main metadata structure)
  gse = getGEO(gse_id, GSEMatrix = FALSE)
  
  # Extract metadata list
  meta = Meta(gse)
  
  # Collapse any vector elements into single strings
  meta_collapsed = map(meta, function(x) {
    if (length(x) > 1) paste(x, collapse = "; ") else x
  })
  
  # Create a single-row tibble
  tibble::as_tibble_row(meta_collapsed)
}



make_sample_metadata = function(geo_id, sample_metadata, model){
  DE_filename = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_DE.tsv")
  if(!file.exists(DE_filename)){
    print("NO DifExpAnalysis RESULTS")
    return()
  }
  attributes = filter(model,`Attribute tier`=="Tier1 (applies to all dataset)")%>%
    filter(!endsWith(`Required`, "Optional"))%>%
    pull(`Attribute name`)
  
  #DE_file = read_tsv(DE_filename)
  metadata = get_metadata(geo_id)
  #print(metadata, width=Inf)#debug
  #more_metadata = get_gse_metadata(geo_id)
  #print(more_metadata, width=Inf)
  
  #same for all GSMs
  Date_exported = format(Sys.Date(), "%m%d%Y") #Was this when I downloaded it, or when I make this?
  Data_contact = "Katherine_McKinney"
  
  
  table = tibble::as_tibble(setNames(rep(list(character()), length(attributes)), attributes))
  
  GSMs = filter(sample_metadata, Dataset_ID == geo_id)%>%
    pull(ID)
  
  
  for(GSM in GSMs){
    DatasetID = geo_id #table says get EMODS ID and name, but chat says NCBI doesn't have one, so from where?
    Dataset_name = NA #how is this different from the id?
    SampleID = GSM
    Data_model_version = NA #Where to find this?
    Script = NA #hwat
    
    if(has_name(metadata, "strain_ch1")){
      X__Sample_Genotype = metadata["strain_ch1"][[1]]
    }else{
      X__Sample_Genotype = NA
    }
    
    status = filter(sample_metadata, ID==GSM)%>%
      pull(Value)
    if(status=="affected_group"){
      X__Sample_Karyotype = "T21"
    }else if(status=="control_group"){
      X__Sample_Karyotype = "Control"
    }else{
      X__Sample_Karyotype = NA
    }
    
    
    table = add_row(table, DatasetID=DatasetID, Dataset_name=Dataset_name, SampleID=SampleID, Data_model_version=Data_model_version, Date_exported=Date_exported, Data_contact=Data_contact, Script=Script, X__Sample_Genotype=X__Sample_Genotype, X__Sample_Karyotype=X__Sample_Karyotype)
  }
  #print(table)
  write_tsv(table, paste0(getwd(), "/Data/Metadata/", geo_id, "_Sample_metadata.tsv"))
}





#--------------metadata_tables-------------

#read in sample and dataset metadata, along with RPKM files, to make the tables

sample_metadata = read_tsv(paste0(getwd(), "/Data/Metadata/SampleMetadata.tsv"))
model = read_excel(paste0(getwd(), "/EMODS_data_model_v0.5.2_dictionary_bulk_RNASeq.xlsx"), sheet=1, skip=2)
print(model)#debug

for (geo_id in names(platforms_list)) {
geo_id = "GSE109293"
  
  #make Sample_metadata
  make_sample_metadata(geo_id, sample_metadata, model)
  
  #make Abundance_data
  
  #make Differential_analysis_results
  
  
}

