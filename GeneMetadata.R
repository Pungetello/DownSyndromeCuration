#-----------loading_libraries-----------
library(tidyverse)
library(biomaRt)

#-----------functions------------------

# returns a vector with the Gene ID's present in all the quality output files
get_gene_ids_from_quality_output = function(filePath) {
  file_list = list.files(path = filePath, pattern="\\S+.tsv.gz", full.names= TRUE, ignore.case = TRUE)
  geneIDs = c()
  
  for (file in file_list){
    file_tibble = read_tsv(file)%>%
      drop_na()%>%
      mutate(Sample_ID = str_sub(Sample_ID, 1, -4))%>%
      dplyr::pull(Sample_ID)
    geneIDs = c(geneIDs, file_tibble)
  }
  
  unique(geneIDs)%>%
    return()
}

#-----------script------------------

#get vector of all Gene ID's from NormalizedData
ensembl_ids = get_gene_ids_from_quality_output(paste0(getwd(), "/Data/NormalizedData"))

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#use biomaRt to get corresponding data
gene_metadata = getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

write_tsv(gene_metadata, paste0(getwd(), "/Data/Metadata/GeneMetadata.tsv"))