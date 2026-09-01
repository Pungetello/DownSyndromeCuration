
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
  metadata = getGEO(geo_ID)[[1]]%>%
    pData()%>%
    as_tibble()%>%
    clean_names()#%>%
    # rename_with(function(Attribute){Attribute = ifelse(endsWith(Attribute, "_ch1"), str_sub(Attribute, 1, -5), Attribute)})
  #metadata = fix_bespoke_issues(geo_ID, metadata) # add back in if any bespoke issues are found
  return (metadata)
}


#TODO: these kinda repeat each other, refactor. But first finish the rest so you can know if/how you actually need them.
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



#in order, searches the metadata for a column where the title matches a key, then returns the value in row i
find_value_from_keys = function(metadata, i, keys){
  for(key in keys){
    if (has_name(metadata, key)) {
        result = metadata[[key]][[i]]
        return(result)
    }
  }
  print("NO MATCH FOUND!")
  return(NA)
}



#makes the study metadata table for the geo_id given according to the model provided
make_study_metadata = function(geo_id, model){
  #make results table
  attributes = dplyr::filter(model,`Table`=="Study_metadata")%>%
    pull(`Attribute name`)
  table = tibble::as_tibble(setNames(rep(list(character()), length(attributes)), attributes))
  
  #prepare metadata
  website_metadata = get_gse_metadata(geo_id)
  
  
  #for each column in results table, search metadata for the proper match or set manually
  StudyID = geo_id
  Study_name = find_value_from_keys(website_metadata, 1, "title")
  Study_description = find_value_from_keys(website_metadata, 1, "overall_design") #or should this be summary? TODO: do both!
  PMID = find_value_from_keys(website_metadata, 1, "pubmed_id")
  Data_model_version = "v0.5.3"
  Date_exported = Date_exported = format(Sys.Date(), "%m%d%Y")
  Data_contact = "Stephen Piccolo" #correct? Not the study's contact?
  Script = "https://github.com/Pungetello/DownSyndromeCuration/blob/main/Main.R"
 
  add_row(table, StudyID, Study_name, Study_description, PMID, Data_model_version, Date_exported, Data_contact, Script)%>%#FIX
    #print()%>%
    return()
}



#makes the dataset metadata for the given geo_id's given dataset. Some repetition from last function, refactor?
make_dataset_metadata = function(geo_id, dataset, metadata, model, GSE_to_SRR){
  #make results table
  attributes = dplyr::filter(model,`Table`=="Dataset_metadata")%>%
    filter(!endsWith(`Required`, "Optional"))%>%
    pull(`Attribute name`)
  table = tibble::as_tibble(setNames(rep(list(character()), length(attributes)), attributes))
  
  #prepare metadata #see if we actually need this later
  GSMs = filter(GSE_to_SRR, GSE == geo_id, Dataset == dataset)%>%
    pull(GSM)
  dataset_metadata = dplyr::filter(metadata, geo_accession %in% GSMs)%>%
    print()
  
  
  #for each column in results table, search metadata for the proper match or set manually
  StudyID = geo_id
  DatasetID = paste0(geo_id, "_", dataset)
  Model_system = find_value_from_keys(dataset_metadata, 1, "organism_ch1")#should this also include the model type, like for the different mice?
  Tissue_cell_type = find_value_from_keys(dataset_metadata, 1, "source_name_ch1")#issue: this is different for each sample. Also doesn't provide enough information, some is in 'title' as well
  Sample_type = find_value_from_keys(dataset_metadata, 1, "molecule_ch1")
  Data_type = "bulk_RNASeq"
  Protocol_details = find_value_from_keys(dataset_metadata, 1, c("treatment_protocol_ch1", "extract_protocol_ch1"))#some have multiple columns of various protocols. Should I combine and use them all?
  Visibility = "default" #indicates whether to display in portal
  
  #These four are the same for all metadata tables. Make global variables?
  Data_model_version = "v0.5.3"
  Date_exported = Date_exported = format(Sys.Date(), "%m%d%Y")
  Data_contact = "Stephen Piccolo" #correct? Not the study's contact?
  Script = "https://github.com/Pungetello/DownSyndromeCuration/blob/main/Main.R"
  
  
  add_row(table, StudyID, DatasetID, Model_system, Tissue_cell_type, Sample_type, Data_type, Protocol_details, Visibility, Data_model_version, Date_exported, Data_contact, Script)%>%#FIX
    #print(width=Inf)%>%
    return()
}



#TODO: REWORK THIS A BUNCH!
#Makes the sample metadata for all of the samples within a given study and dataset, according to the provided model.
make_sample_metadata = function(geo_id, dataset, metadata, model, GSE_to_SRR, sample_metadata){

  #make results table
  attributes = dplyr::filter(model,`Table`=="Sample_metadata")%>%
    filter(!endsWith(`Required`, "Optional"))%>%
    pull(`Attribute name`)
  table = tibble::as_tibble(setNames(rep(list(character()), length(attributes)), attributes))
  
  GSMs = filter(GSE_to_SRR, GSE == geo_id, Dataset == dataset)%>%
    pull(GSM)
  dataset_metadata = dplyr::filter(metadata, geo_accession %in% GSMs)
  
  DatasetID = paste0(geo_id, "_", dataset)
  #these four again
  Data_model_version = "v0.5.3"
  Date_exported = format(Sys.Date(), "%m%d%Y")
  Data_contact = "Stephen Piccolo" #correct? Not the study's contact?
  Script = "https://github.com/Pungetello/DownSyndromeCuration/blob/main/Main.R"
  
  for(i in 1:nrow(metadata)){
    SampleID = find_value_from_keys(metadata, i, "title")

    if(Datasets$Organism[Datasets$Name == geo_id] == "human"){ #might need sample metadata after all? rework this.
      status = filter(sample_metadata, ID==GSM)%>%
        pull(Value)
      if(status=="affected_group"){
        X__Sample_Karyotype = "T21"
      }else if(status=="control_group"){
        X__Sample_Karyotype = "Control"
      }
      table = add_row(table, DatasetID,SampleID,
                      Data_model_version,Date_exported,Data_contact,Script,
                      X__Sample_Genotype)
    }else{
      X__Sample_Genotype = find_value_from_keys(metadata, i, "genotype_ch1")#GSE109293&4, GSE202938, GSE210117
      table = add_row(table, DatasetID,SampleID,
                      Data_model_version,Date_exported,Data_contact,Script,
                      X__Sample_Karyotype)
    }
  }
  print(table)%>%
    return()
  # write_tsv(table, paste0(getwd(), "/Data/Metadata/", geo_id, "_Sample_metadata.tsv"))
}



make_abundance_data = function(geo_id, dataset, metadata, model, GSE_to_SRR){
  
  #make results table
  attributes = filter(model,`Attribute tier`=="Tier1 (applies to all dataset)")%>%
    filter(!endsWith(`Required`, "Optional"))%>%
    pull(`Attribute name`)
  
  table = tibble::as_tibble(setNames(rep(list(character()), length(attributes)), attributes))
  table$Value = as.double(table$Value)
  
  GSMs = filter(GSE_to_SRR, GSE == geo_id, Dataset == dataset)%>%
    pull(GSM)
  
  #get geo_id specific sources
  RPKM_file = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_MAC_RPKM.tsv")
  if(!file.exists(RPKM_file)){
    RPKM_file = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_RPKM.tsv")
  }
  SRRs = read_tsv(RPKM_file)%>%
    colnames()[-1]
  gene_metadata_file = read_tsv(paste0(getwd(), "/Data/Metadata/GeneMetadata/", geo_id, ".tsv.gz")) #from GeneMetadata.R
  
  #define variables for columns that are the same in all rows
  DatasetID = paste0(geo_id, "_", dataset)
  FeatureID_type = "Ensembl"
  Units = "RPKM"
  
  Data_model_version = "v0.5.3"
  Date_exported = format(Sys.Date(), "%m%d%Y")
  Data_contact = "Stephen Piccolo"
  Script = "https://github.com/Pungetello/DownSyndromeCuration/blob/main/Main.R"
  
  for(i in 1:length(GSMs)){
    GSM_i = GSMs[i]
    SRR = filter(GSE_to_SRR, GSM == GSM_i)%>%
      pull("SRR")
    SampleID = filter(metadata, geo_accession == GSM_i)%>%
      pull("title") # different method than in above function, but needed because we can't assume GSMs are in same order
    
    #I need a row for each gene and its count for each GSM. This will be very big.
    FeatureID = pull(RPKM, "gene_id") #Gene/protein/metab identifier. Vector!
    Value = pull(RPKM, SRR)#Feature abundance in sample. Should be vector of same length!
    
    #use gene_metadata_file to get the Feature_name
    GSM_tibble = inner_join(tibble(DatasetID,SampleID,FeatureID,FeatureID_type,Value,Units,Data_model_version,Date_exported,Data_contact,Script), gene_metadata_file, by=join_by("FeatureID"=="ensembl_gene_id"))%>%
      dplyr::rename("Feature_name"="external_gene_name")%>%
      dplyr::select(!c("entrezgene_id","start_position","end_position","chromosome_name","gene_biotype"))
    
    #SRR_tibble = tibble(DatasetID=DatasetID, SampleID=SampleID, FeatureID = FeatureID, FeatureID_type = FeatureID_type, Feature_name = Feature_name, Value = Value, Units = Units, Data_model_version=Data_model_version, Date_exported=Date_exported, Data_contact=Data_contact, Script=Script)
    table = bind_rows(table, GSM_tibble)
  }
  
  print(table)%>%
    return()
  # write_tsv(table, paste0(getwd(), "/Data/Metadata/", geo_id, "_Feature_abundance.tsv.gz"))
}



make_differential_analysis_results = function(geo_id, dataset, model, GSE_to_SRR){
  
  #make results table
  attributes = dplyr::filter(model,`Table`=="Results_DE")%>%
    dplyr::filter(!endsWith(`Required`, "Optional"))%>%
    pull(`Attribute name`)
  
  table = tibble::as_tibble(setNames(rep(list(character()), length(attributes)), attributes))
  #rework: table is made only to be overwritten later
  
  #get geo_id specific sources
  DE_results = read_tsv(paste0(getwd(), "/Data/NormalizedData/", geo_id, "_DE.tsv"))
  
  #define variables for columns that are the same in all rows
  Data_model_version = "v0.5.3"
  Date_exported = format(Sys.Date(), "%m%d%Y")
  Data_contact = "Stephen Piccolo"
  Script = "https://github.com/Pungetello/DownSyndromeCuration/blob/main/Main.R"
  
  DatasetID = paste0(geo_id, "_", dataset)
    
    Statistical_method = "DESeq2"
    Comparison = "Genotype|Dp16_vs_WT"#REWORK!
    Model_specification = NA #"Mathematical model supplied to statistical method"
    FeatureID = pull(DE_results, "gene")
    FeatureID_type = "Ensembl"
    FoldChange = pull(DE_results, "log2FoldChange")
    pvalue = pull(DE_results, "pvalue")
    padj = pull(DE_results, "padj")
    padj_type = "BH FDR"

    table = tibble(DatasetID=DatasetID, Dataset_name=Dataset_name, Statistical_method=Statistical_method, Comparison=Comparison, Model_specification=Model_specification, FeatureID = FeatureID, FeatureID_type = FeatureID_type, FoldChange=FoldChange, pvalue=pvalue, padj=padj, padj_type=padj_type, Data_model_version=Data_model_version, Date_exported=Date_exported, Data_contact=Data_contact, Script=Script)

    gene_metadata_file = read_tsv(paste0(getwd(), "/Data/Metadata/GeneMetadata/", geo_id, ".tsv.gz"))%>% #from GeneMetadata.R
      dplyr::rename("Feature_name"="external_gene_name", "FeatureID"="ensembl_gene_id")%>%
      select("Feature_name","FeatureID")
    
    table = inner_join(table, gene_metadata_file, by = "FeatureID")%>%
      relocate("Feature_name", .after = "FeatureID_type")%>%
      print()
    
    # DatasetID
    # Statistical_method
    # Comparison
    # Model_specification
    # FeatureID
    # FeatureID_typeX
    # Feature_nameX
    # FoldChangeX
    # pvalueX
    # padjX
    # padj_typeX
    # Data_model_version
    # Date_exported
    # Data_contact
    # Script
    
    
  #print(table)
  write_tsv(table, paste0(getwd(), "/Data/Metadata/", geo_id, "_differential_analysis_results.tsv"))
}




#--------------metadata_tables-------------

#read in sample and dataset metadata, along with RPKM files, to make the tables

sample_metadata = read_tsv(paste0(getwd(), "/Data/Metadata/SampleMetadata.tsv"))
#rna_genes = read_tsv(paste0(getwd(), "/Data/rna_gene_data.tsv.gz"))#TODO: this is human one, need to get mouse

model = read_excel(paste0(getwd(), "/EMODS_data_model_v0.5.3_dictionary.xlsx"), sheet=1, skip=0)
GSE_to_SRR = read_tsv(paste0(getwd(), "/Data/RNA_GSE_to_SRR.tsv"))

for (geo_id in pull(Datasets, Name)){
  print(geo_id)
  
  #make Study_metadata. Each GEO = one study. Studies each take only one row, so these are all one-row tsvs.
  make_study_metadata(geo_id, model)%>%
    write_tsv(paste0(getwd(), "/Data/Metadata/", geo_id, "_study_metadata.tsv"))
  
  metadata = get_metadata(geo_id)  
  datasets = filter(GSE_to_SRR, GSE == geo_id)%>%
    pull(Dataset)%>%
    unique()
  
  
  #TODO: do I have the dataset loop out here, or inside each function? Do I append to the file with each new dataset, 
  #or keep massive tibble variables in memory until all datasets are gone through? I think it'll be fine as the latter.
  Dataset_metadata = c()
  Sample_metadata = c()
  Feature_abundance = c()
  Results_DE = c()
  for(dataset in datasets){ #combine all dataset results into one file per GSE
    print(paste0("DATASET ", dataset))
    
    # RPKM_filename = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_", dataset, "_RPKM.tsv")#new implementation, update files to match
    # DE_filename = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_", dataset, "_DE.tsv")
    # if(!file.exists(RPKM_filename) || !file.exists(DE_filename)){
    #   print("NO RPKM or NO DE RESULTS")
    #   next()
    # }
    #TODO: maybe make the dataset loops inside each function, instead of having to combine them all before writing
    
    #make Dataset_metadata
    next_dataset = make_dataset_metadata(geo_id, dataset, metadata, model, GSE_to_SRR)
    Dataset_metadata = c(Dataset_metadata, next_dataset)
    
    #make Sample_metadata
    next_dataset = make_sample_metadata(geo_id, dataset, metadata, model, GSE_to_SRR, sample_metadata)
    Sample_metadata = c(Sample_metadata, next_dataset)
    
    #make Abundance_data
    next_dataset = make_abundance_data(geo_id, dataset, metadata, model, GSE_to_SRR)
    Feature_abundance = c(Feature_abundance, next_dataset)

    #make Differential_analysis_results
    next_dataset = make_differential_analysis_results(geo_id, dataset, model, GSE_to_SRR)
    Results_DE = c(Results_DE, next_dataset)
    
    
  }
  write_tsv(Dataset_metadata, paste0(getwd(), "/Data/Metadata/", geo_id, "_dataset_metadata.tsv"))
  write_tsv(Sample_metadata, paste0(getwd(), "/Data/Metadata/", geo_id, "_sample_metadata.tsv"))
  write_tsv(Feature_abundance, paste0(getwd(), "/Data/Metadata/", geo_id, "_feature_abundance.tsv"))
  write_tsv(Results_DE, paste0(getwd(), "/Data/Metadata/", geo_id, "_results_DE.tsv"))
  #make Results_Pathway_metadata?? Ask Zenitha.
}

