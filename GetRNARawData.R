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



#fasterq-dump command location

#Optimus Prime locations
# fasterq = normalizePath(paste0(getwd(),"/sratoolkit.current-win64/sratoolkit.3.3.0-win64/bin/fasterq-dump.exe"))
# prefetch = normalizePath(paste0(getwd(),"/sratoolkit.current-win64/sratoolkit.3.3.0-win64/bin/prefetch.exe"))

#Supercomputer locations
fasterq = "fasterq-dump"
prefetch = "prefetch"

#TODO: add script to download and install SRA toolkit in location used.



#make sure the SRA toolkit has been downloaded by the user
check_sra = function(){
  if (Sys.which(fasterq) == "") {
    print("SRA Toolkit (fasterq-dump) not found. Please install.")
    
    #TODO: the following will download it but not for all machines so it should probably be removed
    
    #Maybe reconfigure so it can work with different operating systems and SRA toolkit versions?
    # system2("wget", args = c("--output-document", "sratoolkit.tar.gz", "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz"))
    # 
    # system2("tar", args = c("-vxzf", "sratoolkit.tar.gz"))
    # 
    # system2("export", c("PATH=$PATH:~/sratoolkit/sratoolkit.3.3.0-platform/bin"))
    # 
    # system2("source", args = "~/.bashrc")
    # 
    # #test
    # system2("fastrq-dump", args = "--version")
  }
}



get_srr_from_srx <- function(srx_id) {
  # 1. Search the SRA database for the SRX ID
  search <- entrez_search(db = "sra", term = srx_id)
  
  # 2. Get the summary metadata for that ID
  summ <- entrez_summary(db = "sra", id = search$ids)
  
  # 3. Extract the SRR ID (the Run) from the summary
  # The run info is usually stored in the 'runs' string inside the summary
  # This regex grabs anything starting with 'SRR' followed by digits
  srr <- regmatches(summ$runs, regexpr("SRR\\d+", summ$runs))
  return(srr)
}



#get the true raw data using the SRA toolkit
download_raw = function(geo_id){
  print(paste0("PREFETCHING RAW DATA FOR ", geo_id))
  
  gse = getGEO(geo_id, GSEMatrix = FALSE)
  
  srrs <- lapply(GSMList(gse), function(gsm) {
    relations <- Meta(gsm)$relation
    sra_line <- grep("SRA", relations, value = TRUE)
    return(sra_line)
  })
  
  for(link in srrs){
    
    srx = strsplit(link, '=')[[1]][2]
    srr = get_srr_from_srx(srx)
    
    # print(paste(
    #   fasterq,
    #   srr
    #   #"--split-files",
    #   #"--outdir fastq"
    # ))
    
    system2(
      prefetch,
      args = srr)
  }
  
}



#--------------Download_RNA_data-------------

#filter to geo_ids for RNAsec that do not have NormalizedData downloaded. Make sure to run GetRNASecData before this.
for (geo_id in names(platforms_list)){
  platform = platforms_list[[geo_id]]
  
  if(is.na(platform)){
    destination = paste0(getwd(), "/Data/NormalizedData/", geo_id, ".tsv.gz")
    if(!file.exists(destination)){
      
      print(geo_id)
      
      #download SRA toolkit if it cannot be found?
      check_sra()
      
      #prefetch the raw data
      download_raw(geo_id)

    }
  }
}
