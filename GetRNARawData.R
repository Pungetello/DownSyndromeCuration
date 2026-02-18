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

source("PlatformsList.R")
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
  
  srr = regmatches(summary_metadata$runs, regexpr("SRR\\d+", summ$runs))
  
  return(srr)
}



#get the true raw data using the SRA toolkit
download_raw = function(geo_id){
  
  gse = getGEO(geo_id, GSEMatrix = FALSE)
  
  #get SRA line for each GSM for the geo_id
  srrs = lapply(GSMList(gse), function(gsm) {
    relations = Meta(gsm)$relation
    sra_line = grep("SRA", relations, value = TRUE)
    return(sra_line)
  })
  
  for(link in srrs){
    
    #extract SRX from line, convert to SRR
    srx = strsplit(link, '=')[[1]][2]
    srr = get_srr_from_srx(srx)
    
    #download raw data by SRR ID
    system2(
      prefetch,
      args = c(srr, "-O", paste0(getwd(), "/Data/RawRNA/", srr)))
  }
}



#download reference genomes for GRCm39, since that is likely all we need
download_reference = function(){
  if (!dir.exists("RefGenomes")){
    dir.create("RefGenomes", recursive = TRUE)
  }
  
  #download reference genome
  link = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_genomic.gtf.gz"
  destination = paste0(getwd(), "/RefGenomes/GRCm39_ann.gtf.gz")
  if(!file.exists(destination)){
    options(timeout = Inf)
    download.file(link, destination)
  }
  
  #download annotation table
  link = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_rna_from_genomic.fna.gz"
  destination = paste0(getwd(), "/RefGenomes/GRCm39_ref.fna.gz")
  if(!file.exists(destination)){
    options(timeout = Inf)
    download.file(link, destination)
  }
}#TODO: make these functions?



#--------------Download_RNA_data-------------

#filter to geo_ids for RNAsec that do not have NormalizedData downloaded. Make sure to run GetRNASecData before this.
for (geo_id in names(platforms_list)){
  platform = platforms_list[[geo_id]]
  
  if(is.na(platform)){
    destination = paste0(getwd(), "/Data/NormalizedData/", geo_id, ".tsv.gz")
    if(!file.exists(destination)){
      
      #make sure SRA toolkit is downloaded
      check_sra()
      
      #prefetch the raw data
      download_raw(geo_id)
      
      #download reference genome needed
      download_reference()

    }
  }
}
