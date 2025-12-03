#-----------loading_libraries-----------
library(tidyverse)
library(biomaRt)

#-----------functions------------------

# returns a vector with the Gene ID's present in the quality output file
get_gene_ids_from_quality_output = function(file) {
    
  file_tibble = read_tsv(file)%>%
    drop_na()%>%
    mutate(Sample_ID = str_sub(Sample_ID, 1, -4))%>%
    dplyr::pull(Sample_ID)
  geneIDs = c(file_tibble)
  
  unique(geneIDs)%>%
    return()
}

#-----------script------------------

output_file_location = "Data/Metadata/GeneMetadata"

if (!dir.exists(output_file_location)){
  dir.create(output_file_location, recursive = TRUE)
}

filePath = paste0(getwd(), "/Data/NormalizedData")
file_list = list.files(path = filePath, pattern="\\S+.tsv.gz", full.names= TRUE, ignore.case = TRUE)

for(file in file_list){

  #get vector of all Gene ID's from that file
  ensembl_ids = get_gene_ids_from_quality_output(file)
  
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  #use biomaRt to get corresponding data
  if(length(ensembl_ids) > 0){
    gene_metadata = getBM(
      attributes = c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = mart
    )
    
    dataset_name = str_match(file, "\\/([\\w]+[\\d]+)_")[,2]
    write_tsv(gene_metadata, paste0(getwd(), "/Data/Metadata/GeneMetadata/", dataset_name, "_genes.tsv"))
  } else {
    print(paste0("NO GENE IDS IN FILE ", file))
  }

}