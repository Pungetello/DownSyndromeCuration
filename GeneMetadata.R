#-----------loading_libraries-----------
library(tidyverse)
library(biomaRt)

#-----------functions------------------

# TODO: write this
# returns a vector with the Gene ID's present in all the quality output files
get_gene_ids_from_quality_output = function() {
  file_list = list.files(path = file_data, pattern="^[^.]*\\.CEL\\.gz$", full.names= TRUE, ignore.case = TRUE)
  # output_files = read.celfiles(file_list) TODO: find command to read that file type
  
}

#-----------script------------------

#get tibble that's just one column of all Gene ID's from NormalizedData
entrez_ids = get_gene_ids_from_quality_output()

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_metadata = getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "entrezgene_id",
  values = entrez_ids,
  mart = mart
)



