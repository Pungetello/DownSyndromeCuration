#-----------loading_libraries-----------
library(tidyverse)
library(biomaRt)

#-----------functions------------------

# TODO: write this
# returns a vector with the Gene ID's present in all the quality output files
get_gene_ids_from_quality_output = function(filePath) {
  file_list = list.files(path = filePath, pattern="\\S+.tsv.gz", full.names= TRUE, ignore.case = TRUE)
  geneIDs = c()
  
  for (file in file_list){
    file_tibble = read_tsv(file)%>%
      select(Sample_ID)%>%
      drop_na()
    geneIDs = geneIDs + file_tibble
  }
  
}

#-----------script------------------

#get vector of all Gene ID's from NormalizedData
entrez_ids = get_gene_ids_from_quality_output(paste0(getwd(), "/Data/NormalizedData"))

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_metadata = getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "entrezgene_id",
  values = entrez_ids,
  mart = mart
)



