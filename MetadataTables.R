
#-----------loading_libraries-----------

library(tidyverse)
library(readxl)
library(janitor)
library(DBI)
library(GEOquery)
source("Datasets.R")

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


find_value_from_keys = function(metadata, i, keys){
  result = NA
  for (k in keys) {
    if (has_name(metadata, k)) {
      result = metadata[[k]][[i]]
      break
    }
  }
  return(result)
}



make_sample_metadata = function(geo_id, sample_metadata, model){

  #make results table
  attributes = filter(model,`Attribute tier`=="Tier1 (applies to all dataset)")%>%
    pull(`Attribute name`)
  
  table = tibble::as_tibble(setNames(rep(list(character()), length(attributes)), attributes))
  
  GSMs = filter(sample_metadata, Dataset_ID == geo_id)%>%
    pull(ID)
  
  
  
  #get geo_id specific sources
  #DE_file = read_tsv(DE_filename)
  metadata = get_metadata(geo_id)
  print(metadata, width=Inf)#debug
  #more_metadata = get_gse_metadata(geo_id)
  #print(more_metadata, width=Inf)
  
  
  #define variables for columns that are the same in all rows
  Date_exported = format(Sys.Date(), "%m%d%Y") #Was this when I downloaded it, or when I make this?
  Data_contact = "Stephen Piccolo"
  Additional_details = NA
  
  
  for(i in 1:length(GSMs)){
    GSM = GSMs[i]
    DatasetID = geo_id #table says get EMODS ID and name, but chat says NCBI doesn't have one, so from where?
    Dataset_name = NA
    SampleID = GSM
    Data_model_version = NA
    Script = "https://github.com/Pungetello/DownSyndromeCuration/blob/main/Main.R"
    
    
    #add values for other tibble's genotype columns
    X__Sample_Genotype = find_value_from_keys(metadata, i, c("strain_ch1", "genotype_ch1", "background_strain_ch1"))#GSE109293&4, GSE202938, GSE210117
    
    status = filter(sample_metadata, ID==GSM)%>%
      pull(Value)
    if(status=="affected_group"){
      X__Sample_Karyotype = "T21"
    }else if(status=="control_group"){
      X__Sample_Karyotype = "Control"
    }else{
      X__Sample_Karyotype = NA
    }
    
    X__Sample_Treatment = find_value_from_keys(metadata, i, c())#GSE109293&4, GSE202938, GSE210117
    X__Sample_Sex = find_value_from_keys(metadata, i, c())
    X__Sample_Age_group = find_value_from_keys(metadata, i, c())
    X__Sample_age_in_days_post_birth = find_value_from_keys(metadata, i, c())
    X__Sample_age_in_days_post_conception = find_value_from_keys(metadata, i, c("source_name_ch1"))
    X__Sample_Age_in_weeks = find_value_from_keys(metadata, i, c())
    X__Sample_Harvest_batch = find_value_from_keys(metadata, i, c())
    X__Sample_Cell_type = find_value_from_keys(metadata, i, c("cell_type_ch1"))
    X__Sample_Cell_line = find_value_from_keys(metadata, i, c())
    X__Sample_DonorID = find_value_from_keys(metadata, i, c())
    X__Sample_Batch = find_value_from_keys(metadata, i, c())
    X__Sample_Comparison_group = find_value_from_keys(metadata, i, c())#is this anything?
    X__Sample_tatoo = find_value_from_keys(metadata, i, c())
    X__Sample_CageID = find_value_from_keys(metadata, i, c())
    X__Sample_Donor_cell_type = find_value_from_keys(metadata, i, c())
    X__Sample_Pre_differentation_passage_number = find_value_from_keys(metadata, i, c())
    X__Sample_Post_differentation_passage_number = find_value_from_keys(metadata, i, c())
    X__Sample_Differentiation_factors = find_value_from_keys(metadata, i, c())
    X__Sample_Days_in_differentiation = find_value_from_keys(metadata, i, c())
    X__Sample_Oxygen_percent = find_value_from_keys(metadata, i, c())
    X__Sample_Coating_matrix = find_value_from_keys(metadata, i, c())
    
    
    table = add_row(table, DatasetID = DatasetID,
                    Dataset_name = Dataset_name,
                    SampleID = SampleID,
                    Additional_details = Additional_details,
                    Data_model_version = Data_model_version,
                    Date_exported = Date_exported,
                    Data_contact = Data_contact,
                    Script = Script,
                    X__Sample_Genotype = X__Sample_Genotype,
                    X__Sample_Karyotype = X__Sample_Karyotype,
                    X__Sample_Treatment = X__Sample_Treatment,
                    X__Sample_Sex = X__Sample_Sex,
                    X__Sample_Age_group = X__Sample_Age_group,
                    X__Sample_age_in_days_post_birth = X__Sample_age_in_days_post_birth,
                    X__Sample_age_in_days_post_conception = X__Sample_age_in_days_post_conception,
                    X__Sample_Age_in_weeks = X__Sample_Age_in_weeks,
                    X__Sample_Harvest_batch = X__Sample_Harvest_batch,
                    X__Sample_Cell_type = X__Sample_Cell_type,
                    X__Sample_Cell_line = X__Sample_Cell_line,
                    X__Sample_DonorID = X__Sample_DonorID,
                    X__Sample_Batch = X__Sample_Batch,
                    X__Sample_Comparison_group = X__Sample_Comparison_group,
                    X__Sample_tatoo = X__Sample_tatoo,
                    X__Sample_CageID = X__Sample_CageID,
                    X__Sample_Donor_cell_type = X__Sample_Donor_cell_type,
                    X__Sample_Pre_differentation_passage_number = X__Sample_Pre_differentation_passage_number,
                    X__Sample_Post_differentation_passage_number = X__Sample_Post_differentation_passage_number,
                    X__Sample_Differentiation_factors = X__Sample_Differentiation_factors,
                    X__Sample_Days_in_differentiation = X__Sample_Days_in_differentiation,
                    X__Sample_Oxygen_percent = X__Sample_Oxygen_percent,
                    X__Sample_Coating_matrix = X__Sample_Coating_matrix)
                    
  }
  #remove any attributes that are all NA's
  table = table[, colSums(!is.na(table)) > 0]
  #print(table)
  write_tsv(table, paste0(getwd(), "/Data/Metadata/", geo_id, "_Sample_metadata.tsv"))
}



make_abundance_data = function(geo_id, model, mouse_genes){
  
  #make results table
  attributes = filter(model,`Attribute tier`=="Tier1 (applies to all dataset)")%>%
    filter(!endsWith(`Required`, "Optional"))%>%
    pull(`Attribute name`)
  
  table = tibble::as_tibble(setNames(rep(list(character()), length(attributes)), attributes))
  table$Value = as.double(table$Value)
  
  GSMs = filter(sample_metadata, Dataset_ID == geo_id)%>%
    pull(ID)
  
  #get geo_id specific sources
  RPKM_file = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_RPKM.tsv")
  RPKM = read_tsv(RPKM_file)
  SRRs = colnames(RPKM)[-1]
  gene_metadata_file = read_tsv(paste0(getwd(), "/Data/Metadata/GeneMetadata/", geo_id, ".tsv.gz")) #from GeneMetadata.R
  
  #define variables for columns that are the same in all rows
  Date_exported = format(Sys.Date(), "%m%d%Y") #Was this when I downloaded it, or when I make this?
  Data_contact = "Stephen Piccolo"
  FeatureID_type = "Ensembl"
  DatasetID = geo_id #table says get EMODS ID and name, but chat says NCBI doesn't have one, so from where?
  Dataset_name = NA #how is this different from the id?
  Data_model_version = NA
  Script = "https://github.com/Pungetello/DownSyndromeCuration/blob/main/Main.R"
  
  for(SRR in SRRs){
    SampleID = SRR#should I convert these all to the GSMs to make it consistent?
    
    #Do I need a row for each gene and its count for each GSM? This will be very big.
    FeatureID = pull(RPKM, "gene_id") #Gene/protein/metab identifier. Vector!
    Value = pull(RPKM, SRR)#Feature abundance in sample. Should be vector of same length!
    Units = "RPKM"
    
    SRR_tibble = inner_join(tibble(DatasetID,Dataset_name,SampleID,FeatureID,FeatureID_type,Value,Units,Data_model_version,Date_exported,Data_contact,Script), gene_metadata_file, by=join_by("FeatureID"=="ensembl_gene_id"))%>%
      rename("X__Feature_chromosome"="chromosome_name", "X__Feature_gene_type"="gene_biotype", "Feature_name"="hgnc_symbol")%>%
      dplyr::select(!c("entrezgene_id","start_position","end_position"))
    
    #SRR_tibble = tibble(DatasetID=DatasetID, Dataset_name=Dataset_name, SampleID=SampleID, FeatureID = FeatureID, FeatureID_type = FeatureID_type, Feature_name = Feature_name, Value = Value, Units = Units, Data_model_version=Data_model_version, Date_exported=Date_exported, Data_contact=Data_contact, Script=Script)
    table = bind_rows(table, SRR_tibble)
  }
  
  #print(table)
  write_tsv(table, paste0(getwd(), "/Data/Metadata/", geo_id, "_Abundance_data.tsv.gz"))
}



make_differential_analysis_results = function(geo_id, model){
  
  #make results table
  attributes = filter(model,`Attribute tier`=="Tier1 (applies to all dataset)")%>%
    filter(!endsWith(`Required`, "Optional"))%>%
    pull(`Attribute name`)
  
  table = tibble::as_tibble(setNames(rep(list(character()), length(attributes)), attributes))
  
  #get geo_id specific sources
  DE_file = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_DE.tsv")
  DE_results = read_tsv(DE_file)
  
  #define variables for columns that are the same in all rows
  Date_exported = format(Sys.Date(), "%m%d%Y") #Was this when I downloaded it, or when I make this?
  Data_contact = "Stephen Piccolo"
  DatasetID = geo_id #table says get EMODS ID and name, but chat says NCBI doesn't have one, so from where?
  Dataset_name = NA #how is this different from the id?
    
    
    #Do I need a row for each gene and its count?
    FeatureID = NA#Gene/protein/metab identifier
    FeatureID_type = NA#Identifier type - Entity|Database
    Feature_name = NA#Feature name/symbol
    Value = NA#Feature abundance in sample
    Units = NA#Feature abundance metric unit
    
    Data_model_version = NA
    Script = "https://github.com/Pungetello/DownSyndromeCuration/blob/main/Main.R"
    
    Statistical_method = "DESeq2"
    Comparison = "Genotype|Dp16_vs_WT"
    Model_specification = NA
    FeatureID = pull(DE_results, "gene")
    FeatureID_type = "Ensembl"
    Feature_name = NA#gene symbol?
    FoldChange = pull(DE_results, "log2FoldChange")
    pvalue = pull(DE_results, "pvalue")
    padj = pull(DE_results, "padj")
    padj_type = "BH FDR"
    X__Feature_chromosome = NA#see above
    X__Feature_gene_type = NA#

    table = tibble(DatasetID=DatasetID, Dataset_name=Dataset_name, Statistical_method=Statistical_method, Comparison=Comparison, Model_specification=Model_specification, FeatureID = FeatureID, FeatureID_type = FeatureID_type, Feature_name = Feature_name, FoldChange=FoldChange, pvalue=pvalue, padj=padj, padj_type=padj_type, Data_model_version=Data_model_version, Date_exported=Date_exported, Data_contact=Data_contact, Script=Script)

  #print(table)
  write_tsv(table, paste0(getwd(), "/Data/Metadata/", geo_id, "_differential_analysis_results.tsv"))
}




#--------------metadata_tables-------------

#read in sample and dataset metadata, along with RPKM files, to make the tables

sample_metadata = read_tsv(paste0(getwd(), "/Data/Metadata/SampleMetadata.tsv"))
#rna_genes = read_tsv(paste0(getwd(), "/Data/rna_gene_data.tsv.gz"))#TODO: this is human one, need to get mouse

SM_model = read_excel(paste0(getwd(), "/EMODS_data_model_v0.5.2_dictionary_bulk_RNASeq.xlsx"), sheet=1, skip=2)
AD_model = read_excel(paste0(getwd(), "/EMODS_data_model_v0.5.2_dictionary_bulk_RNASeq.xlsx"), sheet=2, skip=2)
DAR_model = read_excel(paste0(getwd(), "/EMODS_data_model_v0.5.2_dictionary_bulk_RNASeq.xlsx"), sheet=3, skip=2)
#print(AD_model)#debug

for (geo_id in pull(Datasets, Name)){
  #geo_id = "GSE109293"
  RPKM_filename = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_RPKM.tsv")
  DE_filename = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_DE.tsv")
  if(!file.exists(RPKM_filename) || !file.exists(DE_filename)){
    print("NO RPKM or NO DE RESULTS")
    next()
  }
  
  
  #make Sample_metadata
  #make_sample_metadata(geo_id, sample_metadata, SM_model)
  
  #make Abundance_data
  make_abundance_data(geo_id, AD_model)
  
  #make Differential_analysis_results
  #make_differential_analysis_results(geo_id, DAR_model)
  
  
}

