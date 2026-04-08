#-----------loading_libraries-----------
library(tidyverse)
library(biomaRt)
source("Datasets.R")



#-----------functions------------------

# returns a vector with the Gene ID's present in the quality output file
get_genes_from_quality_output = function(in_file, column) {
  
  file_tibble = read_tsv(in_file)%>%
    dplyr::pull(column)
  geneIDs = c(file_tibble)
  
  unique(geneIDs)%>%
    return()
}



#use biomaRt to get corresponding data
get_gene_metadata = function(values, filters, organism){
  
  if(organism == "hs"){
    mart = useEnsembl(biomart = "genes",
                      dataset = "hsapiens_gene_ensembl",
                      mirror = "useast")
  }else{
    mart = useEnsembl(biomart = "genes",
                      dataset = "mmusculus_gene_ensembl",
                      mirror = "useast")
  }
  
  
  if(length(values) > 0){
    gene_metadata = getBM(
      attributes = c("entrezgene_id","hgnc_symbol","ensembl_gene_id","chromosome_name","start_position","end_position","gene_biotype"),
      filters = filters,
      values = values,
      mart = mart
    )
  }
  return(gene_metadata)
}



#get gene ID's from file, get metadata, combine and save
save_metadata_file = function(in_file, destination, column_title, id_type, organism){
  if (!file.exists(in_file)){
    print(paste0("NORMALIZED DATA FILE NOT FOUND"))
    return()
  }
  
  #get vector of all Gene ID's from that file
  gene_ids = get_genes_from_quality_output(in_file, column_title)
  
  if(length(gene_ids) > 0){
    
    gene_metadata = get_gene_metadata(gene_ids, id_type, organism)
    write_tsv(gene_metadata, destination)
    
  } else {
    print(paste0("NO GENE IDS IN FILE ", in_file))
  }
  
}



#replace the id column in the normalized data with ensembl id's so it's same as affymetrix
replace_id_col = function(normalized_file, metadata_file){
  ids_tibble = read_tsv(metadata_file)%>%
    dplyr::select(entrezgene_id, ensembl_gene_id)
  
  normalized_data = read_tsv(normalized_file)
  
  fixed_data = inner_join(ids_tibble, normalized_data, by = c("entrezgene_id" = "GeneID"))%>%
    dplyr::select(!"entrezgene_id")%>%
    rename("SampleID" = "ensembl_gene_id")
  
  write_tsv(fixed_data, normalized_file)
}



#-----------script------------------

output_file_location = "Data/Metadata/GeneMetadata"

if (!dir.exists(output_file_location)){
  dir.create(output_file_location, recursive = TRUE)
}

for (geo_id in pull(Datasets, Name)){
  print(geo_id) #debug
  destination = paste0(getwd(), "/Data/Metadata/GeneMetadata/", geo_id, ".tsv.gz")
  
  #makes sure it doesn't remake it
  if (file.exists(destination)){next}
  
  
  if(Datasets$Organism[Datasets$Name == geo_id] == "human"){
    #Human
    in_file = paste0(getwd(), "/Data/NormalizedData/", geo_id, ".tsv.gz")
    if(Datasets$Type[Datasets$Name == geo_id] == "RNA"){
      save_metadata_file(in_file, destination, "SampleID", "ensembl_gene_id", "hs")
    }else{
      save_metadata_file(in_file, destination, "Sample_ID", "ensembl_gene_id", "hs")
    }
    

  }else{
    #Mouse
    in_file = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_gene_counts.csv")#could also use DE or RPKM or any of the others
    save_metadata_file(in_file, destination, "gene_id", "ensembl_gene_id", "mm")

    #replace_id_col(in_file, destination)

  }
}