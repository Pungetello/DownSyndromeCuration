#-----------loading_libraries-----------
library(tidyverse)
library(RCurl)
library(biomaRt)

source("PlatformsList.R")

#-----------functions------------------

download_quality_output = function(geo_id){
  
  # DOWNLOADS FOR 101942,190053
  # NO NCBI GENERATED DATA FOR 109293,109294,184771,202938,210117
  
  link = paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=", geo_id, "&format=file&file=", geo_id, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  
  destination = paste0(getwd(), "/Data/NormalizedData/RNA/", geo_id, ".tsv.gz")
  
  if(!file.exists(destination) && RCurl::url.exists(link)){
    print(paste0("FILE EXISTS FOR ", geo_id))
    #check if link exists
    options(timeout = Inf)
    download.file(link, destination)
  }
}


download_gene_data = function(){ #maybe add some parameters
  #download the human one  
  link = "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
  destination = file.path(getwd(), "Data/rna_gene_data.tsv.gz")
  if(!file.exists(destination)){
    #check if link exists
    options(timeout = Inf)
    download.file(link, destination)
  }
  return(destination)
}

#TODO: refactor so you're not repeating from GeneMetadata!
#TODO: save gene file in all the places!

# returns a vector with the Gene symbols present in the file
get_gene_symbols_from_quality_output = function(file) {
  
  file_tibble = read_tsv(file)%>%
    drop_na()%>%
    dplyr::pull(Symbol)
  geneIDs = c(file_tibble)
  
  unique(geneIDs)%>%
    return()
}


#-----------Get_RNASec_Data-----------

if (!dir.exists("Data/NormalizedData/RNA")){
  dir.create("Data/NormalizedData/RNA", recursive = TRUE)
}

for (geo_id in names(platforms_list)){
  platform = platforms_list[[geo_id]]
  if (is.na(platform)){
    download_quality_output(geo_id)
  }
}

gene_file = download_gene_data()

#get vector of all Gene symbols from that file
gene_symbols = get_gene_symbols_from_quality_output(gene_file)

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#use biomaRt to get corresponding data
if(length(gene_symbols) > 0){
  gene_metadata = getBM(
    attributes = c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"), #cols are in a different order for some reason?
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart
  )
  
  dataset_name = "rna"#str_match(file, "\\/([\\w]+[\\d]+)_")[,2]
  write_tsv(gene_metadata, paste0(getwd(), "/Data/Metadata/GeneMetadata/", dataset_name, "_genes.tsv"))
} else {
  print(paste0("NO GENE IDS IN FILE ", file))
}
