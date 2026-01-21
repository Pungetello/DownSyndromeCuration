#-----------loading_libraries-----------
library(tidyverse)
library(biomaRt)

source("PlatformsList.R")

#-----------functions------------------

# returns a vector with the Gene ID's present in the quality output file
get_genes_from_quality_output = function(file, column) {
  
  file_tibble = read_tsv(file)%>%
    dplyr::pull(column)
  geneIDs = c(file_tibble)
  
  unique(geneIDs)%>%
    return()
}



#use biomaRt to get corresponding data
get_gene_metadata = function(values, filters){
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  if(length(values) > 0){
    gene_metadata = getBM(
      attributes = c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
      filters = filters,
      values = values,
      mart = mart
    )
  }
  return(gene_metadata)
}



#-----------script------------------

output_file_location = "Data/Metadata/GeneMetadata"

if (!dir.exists(output_file_location)){
  dir.create(output_file_location, recursive = TRUE)
}

filePath = paste0(getwd(), "/Data/NormalizedData")

for (geo_id in names(platforms_list)){
  destination = paste0(getwd(), "/Data/Metadata/GeneMetadata/", geo_id, ".tsv.gz")
  
  #makes sure it doesn't redownload
  if (file.exists(destination)){next}
    
  platform = platforms_list[[geo_id]]
  
  if (!is.na(platform)){
    #Affymetrix
    file = paste0(getwd(), "/Data/NormalizedData/", geo_id, ".tsv.gz")
    
    if (!file.exists(file)){
      print("NORMALIZED DATA FILE DOES NOT EXIST")
      next
    }

    #get vector of all Gene ID's from that file
    ensembl_ids = get_genes_from_quality_output(file, "Sample_ID")
    if(length(gene_symbols) > 0){
    
      gene_metadata = get_gene_metadata(ensembl_ids, "ensembl_gene_id")
      write_tsv(gene_metadata, )
      
    } else {
      print(paste0("NO GENE IDS IN FILE ", file))
    }

  }else{
    #RNA Sec
    file = paste0(getwd(), "/Data/rna_gene_data.tsv.gz") #TODO: should we be reading the list from the data? if so, how?
    
    if (!file.exists(file)){next}
    
    gene_symbols = get_genes_from_quality_output(file, "Symbol")
    if(length(gene_symbols) > 0){
      
      gene_metadata = get_gene_metadata(gene_symbols, "hgnc_symbol")
      write_tsv(gene_metadata, paste0(getwd(), "/Data/Metadata/GeneMetadata/", geo_id, ".tsv.gz"))
      
    } else {
      print(paste0("NO GENE IDS IN FILE ", file))
    }
  }
}