#change library path for supercomputer
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "" || file.access(user_lib, 2) != 0) {
  user_lib <- "~/R/library"
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
}


#-----------loading_libraries-----------
library(GEOquery)
library(rentrez)
library(tidyverse)
source("Datasets.R")
#source ~/.bashrc


#----------functions-------------

#command locations

#Optimus Prime locations
# fasterq = normalizePath(paste0(getwd(),"/sratoolkit.current-win64/sratoolkit.3.3.0-win64/bin/fasterq-dump.exe"))
# prefetch = normalizePath(paste0(getwd(),"/sratoolkit.current-win64/sratoolkit.3.3.0-win64/bin/prefetch.exe"))

#Supercomputer locations
fasterq = "fasterq-dump"
prefetch = "prefetch"

#TODO: add script to download and install SRA toolkit in location used.


#make sure the SRA toolkit has been downloaded by the user
check_sra = function(){
  if (Sys.which("fasterq-dump") == "") {
    stop(
      "SRA Toolkit (fasterq-dump) not found.\n",
      "Please install it:\n",
      "  macOS: brew install sra-tools\n",
      "  Linux: conda install -c bioconda sra-tools\n",
      "  Windows: https://github.com/ncbi/sra-tools"
    )
  }
}


#Search SRA database for the SRX ID and return the SRR ID from its metadata
get_srr_from_srx = function(srx_id) {
  
  search_result = entrez_search(db = "sra", term = srx_id)
  summary_metadata = entrez_summary(db = "sra", id = search_result$ids)
  
  srr = regmatches(summary_metadata$runs, regexpr("SRR\\d+", summary_metadata$runs))
  
  return(srr)
}



#creates the file mapping all RNA GSE's to their GSM's, SRX's and SRR's. Assumes each GSM has 1 SRX and 1 SRR.
create_GSE_to_SRR = function(datasets_table){
  GSE_to_SRR = tibble(GSE = character(), GSM = character(), SRX = character(), SRR = character())
  
  for (geo_id in pull(datasets_table, Name)){
    if(Datasets$Type[Datasets$Name == geo_id] == "RNA"){
      gse = getGEO(geo_id, GSEMatrix = FALSE)
      
      #get SRA line for each GSM for the geo_id
      for(gsm in GSMList(gse)) {
        relations = Meta(gsm)$relation
        sra_line = grep("SRA", relations, value = TRUE)
        
        #extract SRX from line, convert to SRR
        srx = strsplit(sra_line, '=')[[1]][2]
        srr = get_srr_from_srx(srx)
        
        GSE_to_SRR = add_row(GSE_to_SRR, GSE=geo_id, GSM=Meta(gsm)$geo_accession, SRX=srx, SRR=srr)
      }
    }
  }
  
  #write dataframe of GSE mapped to each SRR to file
  write_tsv(GSE_to_SRR, "Data/RNA_GSE_to_SRR.tsv")
}



#get the true raw data using the SRA toolkit
download_raw = function(geo_id){
  
  GSE_to_SRR = read_tsv(paste0(getwd(), "/Data/RNA_GSE_to_SRR.tsv"))
  
  srrs = filter(GSE_to_SRR, GSE==geo_id)%>%
    pull(SRR)
  
  for(srr in srrs){
    
    if(!file.exists(paste0(getwd(), "/Data/RawRNA/", srr, "/", srr, ".sra"))){
      print(paste0(getwd(), "/Data/RawRNA/", srr, "NOT FOUND, DOWNLOADING RAW DATA FOR ", srr))
    
      #download raw data by SRR ID if not done already
      system2(
        prefetch,
        args = c(srr, "-O", paste0(getwd(), "/Data/RawRNA")))
    }
  }
}



safe_download = function(link, destination){
  #download file from link only if nothing at destination yet
  if(!file.exists(destination)){
    options(timeout = Inf)
    download.file(link, destination)
  }
}



#download reference genomes for GRCm39, C57BL_6J, and DBA_2J
download_reference = function(){
  if (!dir.exists("RefGenomes")){
    dir.create("RefGenomes", recursive = TRUE)
  }
  
  #download reference genome
  link = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz"
  destination = paste0(getwd(), "/RefGenomes/GRCm39_ref.fna.gz")
  safe_download(link, destination)
  
  #download annotation table
  link = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz"
  destination = paste0(getwd(), "/RefGenomes/M38_ann.gtf.gz")
  safe_download(link, destination)
  
  #human files
  #https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
  #https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
}



#--------------Download_RNA_data-------------

#create a file mapping all GSE's in platforms_list to their respective GSM's, SRX's and SRR's.
create_GSE_to_SRR(Datasets)


#filter to geo_ids for RNAsec that do not have NormalizedData downloaded. Make sure to run GetRNASecData before this.
for (geo_id in pull(Datasets, Name)){
  if(Datasets$Type[Datasets$Name == geo_id] == "RNA"){
    print(geo_id)
    destination = paste0(getwd(), "/Data/NormalizedData/", geo_id, ".tsv.gz")
    #if(!file.exists(destination)){ #for filtering out human RNA

      #make sure SRA toolkit is downloaded
      check_sra()

      #prefetch the raw data
      download_raw(geo_id)

      #download reference genome needed
      download_reference()

    }
  }
#}
