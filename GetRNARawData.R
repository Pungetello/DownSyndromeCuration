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

library(GenomicRanges)
library(Biostrings)
library(GenomeInfoDb)
library(Rsamtools)
library(rtracklayer)

# library(jsonlite)

# library(ArrayExpress)
library(curl)


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
  
  srr = regmatches(summary_metadata$runs, gregexpr("SRR\\d+", summary_metadata$runs))
  
  return(unlist(srr))
}



#creates the file mapping all RNA GSE's to their GSM's, SRX's and SRR's. Assumes each GSM has 1 SRX, but SRX can have many SRR.
create_GSE_to_SRR = function(datasets_table){
  GSE_to_SRR = tibble(GSE = character(), GSM = character(), SRX = character(), SRR = character())
  
  for (geo_id in pull(datasets_table, Name)){
    print(geo_id)
    if(Datasets$Type[Datasets$Name == geo_id] == "RNA"){
      gse = getGEO(geo_id, GSEMatrix = FALSE)
      
      #get SRA line for each GSM for the geo_id
      for(gsm in GSMList(gse)) {
        relations = Meta(gsm)$relation
        sra_line = grep("SRA", relations, value = TRUE)
        
        #extract SRX from line, convert to SRR
        print(sra_line)#debug
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
download_raw_geo = function(geo_id){
  
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



#get the raw data for E-MTAB id types. Works, but takes a long time.
#TODO: (if needed) make faster and also compressed
download_raw_emtab = function(id){
  if (!dir.exists("fastq")){
    dir.create("fastq", recursive = TRUE)
  }
  
  #get list of ERR's
  ERRs = get_ERRs(id)
  #ERRs = c("ERR3305905")#debug
  print(ERRs)
  
  for(ERR in ERRs){
    file_report = curl_download(paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", ERR, "&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,bam_ftp&format=tsv&download=true&limit=0"), "temp")
    print(file_report)
    
    file_report = read_tsv(file_report)
    file.remove("temp")
    print(file_report)
    
    fastq_urls = strsplit(file_report$fastq_ftp[1], ";")
    print(fastq_urls)

    for (i in 1:length(fastq_urls[[1]])) {
      url = fastq_urls[[1]][i]
      print(basename(url))
      curl_download(url, destfile = paste0(getwd(), "/fastq/", basename(url)))
      h = new_handle(dirlistonly=TRUE)
      con = curl(url, "r", h)
    }
  }
}



#get a list of all ERRs from a given E-MTAB project to query download of raw data
get_ERRs = function(emtab){

  link = paste0("https://www.ebi.ac.uk/biostudies/files/", emtab, "/", emtab, ".sdrf.txt")
  
  #print("downloading SDRF")
  sdrf = read.delim(link, check.names = FALSE)
  
  errs = unique(sdrf$`Comment[ENA_RUN]`)
  
  #print("RETURNING ERRS")
  return(errs)
}


# library(jsonlite)
# library(readr)
# library(stringr)
# library(dplyr)

# get_ERRs = function(emtab){
#   
#   ## Get BioStudies metadata
#   meta <- fromJSON(
#     sprintf(
#       "https://www.ebi.ac.uk/biostudies/api/v1/studies/%s",
#       emtab
#     ),
#     simplifyVector = FALSE
#   )
#   
#   ## Find the SDRF file
#   sdrf <- Filter(
#     \(x) grepl("\\.sdrf\\.txt$", x$path),
#     meta$files
#   )[[1]]$path
#   
#   ## Download it
#   sdrf_url <- paste0(
#     "https://www.ebi.ac.uk/biostudies/files/",
#     emtab,
#     "/",
#     sdrf
#   )
#   
#   dat <- read_tsv(sdrf_url, show_col_types = FALSE)
#   
#   ## Find the column containing ERRs
#   err_col <- names(dat)[
#     sapply(dat, function(x)
#       any(str_detect(x, "^ERR[0-9]+$"), na.rm = TRUE))
#   ][1]
#   
#   unique(dat[[err_col]])
# }



safe_download = function(link, destination){
  #download file from link only if nothing at destination yet
  if(!file.exists(destination)){
    options(timeout = Inf)
    download.file(link, destination)
  }
}



#download reference genomes for GRCm39 and human
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
  link = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"
  destination = paste0(getwd(), "/RefGenomes/GRCh38_ref.fna.gz")
  safe_download(link, destination)
  
  link = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
  destination = paste0(getwd(), "/RefGenomes/49_ann.gtf.gz")
  safe_download(link, destination)
}



#append necessary sections of human chr21 to the mouse refrence genome
create_mac_reference = function(){
  #get coords of needed sequences within human genome
  gene_mapping = read_tsv(paste0(getwd(), "/MouseModel_GeneMapping_v0.6_ZS_TcMAC21_genes.txt"))
  region_coords = filter(gene_mapping, region_type == "MAC HSA21q")%>%
    select(start, end)
  
  deletion_coords = filter(gene_mapping, region_type == "deletion")%>%
    select(start, end)
  
  #create ranges
  mac_full = GRanges(seqnames = "chr21", ranges = IRanges(start = region_coords$start[1], end = region_coords$end[1]))
  deletions = GRanges(seqnames = "chr21", ranges = IRanges(start = pull(deletion_coords, start), end = pull(deletion_coords, end)))
  mac_fragments = setdiff(mac_full, deletions)
  
  #make sure it doesn't already exist
  if(!file.exists(paste0(getwd(), "/RefGenomes/mouse_plus_mac_sans16.fa"))){
    
    #read in human reference genome
    genome = readDNAStringSet(paste0(getwd(), "/RefGenomes/GRCh38_ref.fna.gz"))
    names(genome) <- sub(" .*", "", names(genome))
    
    #extract sequences
    seqs <- DNAStringSet(lapply(seq_along(mac_fragments), function(i){
      chr <- as.character(seqnames(mac_fragments)[i])
      start <- start(mac_fragments)[i]
      end <- end(mac_fragments)[i]
      
      subseq(genome[[chr]], start = start, end = end)
    }))
    
    names(seqs) = paste0("MAC_", c("1","2","3","4","5"))
    
    #append to copy of mouse reference genome
    mac_file = paste0(getwd(), "/RefGenomes/mac_sequences.fa")
    mouse_file = paste0(getwd(), "/RefGenomes/GRCm39_ref.fna.gz")
    combined_file = paste0(getwd(), "/RefGenomes/mouse_plus_mac.fa")
    
    writeXStringSet(seqs, mac_file)
    
    mouse <- readDNAStringSet(mouse_file)
    mac <- readDNAStringSet(mac_file)
    
    # #Detective work: filter out chr16
    # print(names(mouse))
    # mouse_filtered <- mouse[!names(mouse) %in% "chr16 16"]
    # print(names(mouse_filtered))
    
    combined <- c(mac, mouse)
    
    writeXStringSet(combined, combined_file)
  }
  
  return(mac_fragments)
}



#append relevant human genes in mac to the mouse annotation file
create_mac_annotation = function(mac_fragments){
  print(mac_fragments)
  
  #read in mouse and human annotation tables
  mouse_gtf <- import(paste0(getwd(), "/RefGenomes/M38_ann.gtf.gz"))
  human_gtf <- import(paste0(getwd(), "/RefGenomes/49_ann.gtf.gz"))
  
  #extract MAC sections of human annotation table
  mac1 <- subsetByOverlaps(human_gtf, mac_fragments[1])
  seqlevels(mac1) <- c(chr21 = "MAC_1")
  mac2 <- subsetByOverlaps(human_gtf, mac_fragments[2])
  seqlevels(mac2) <- c(chr21 = "MAC_2")
  mac3 <- subsetByOverlaps(human_gtf, mac_fragments[3])
  seqlevels(mac3) <- c(chr21 = "MAC_3")
  mac4 <- subsetByOverlaps(human_gtf, mac_fragments[4])
  seqlevels(mac4) <- c(chr21 = "MAC_4")
  mac5 <- subsetByOverlaps(human_gtf, mac_fragments[5]) 
  seqlevels(mac5) <- c(chr21 = "MAC_5")
  
  human_subset = c(mac1, mac2, mac3, mac4, mac5)
  
  print(seqlevels(human_subset))
  print(human_subset)
  
  #combine
  combined_gtf <- c(human_subset, mouse_gtf)
  
  #Detective work: remove all mouse chr16 genes
  # combined_gtf_sans16 = combined_gtf[seqnames(combined_gtf) != "chr16"]
  
  export(combined_gtf, paste0(getwd(), "/RefGenomes/mouse_plus_mac.gtf.gz"))
}



#--------------Download_RNA_data-------------

#create a file mapping all GSE's in platforms_list to their respective GSM's, SRX's and SRR's.
# create_GSE_to_SRR(Datasets[1:67, ])
# print("DONE!")
# stop()

#filter to geo_ids for RNAsec that do not have NormalizedData downloaded. Make sure to run GetRNASecData before this.
for (geo_id in pull(Datasets, Name)){

  # #skip if in a dataset being processed currently
  # if(sum(geo_id == c("GSE154418",
  #                "GSE160637",
  #                "GSE160690",
  #                "GSE166849",
  #                "GSE151282")) > 0){
  #   next()
  # }
  #skip if in a dataset being processed currently
  if(sum(geo_id == c("GSE109293","GSE109294","GSE202938","GSE210117")) == 0){
    print("NOT SELECTED GSE, SKIPPING")
    next()
  }

  if(Datasets$Type[Datasets$Name == geo_id] == "RNA"){ #&& Datasets$Organism[Datasets$Name == geo_id] == "mouse"){
    print(geo_id)
    human_destination = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_gene_counts.tsv")
    mouse_destination = paste0(getwd(), "/Data/NormalizedData/", geo_id, "_MAC_fixed_gene_counts.tsv")
    # if(!file.exists(human_destination) && !file.exists(mouse_destination)){ #skip those already processed

      print("DOWNLOADING RAW DATA")
      #prefetch the raw data
      if(startsWith(geo_id, "GSE")){
        #make sure SRA toolkit is downloaded
        check_sra()
        download_raw_geo(geo_id)

      }else if(startsWith(geo_id, "E-MTAB")){
        download_raw_emtab(geo_id)
      }
    # }
  }
}

# #download reference genomes needed
# download_reference()
# 
# #create MAC combined reference genome
# if(!file.exists(paste0(getwd(), "RefGenomes/mouse_plus_mac.fa"))){
#   mac_fragments = create_mac_reference()
#   create_mac_annotation(mac_fragments)
# }


# library(curl)
# url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR330/005/ERR3305905/ERR3305905_1.fastq.gz"
# curl_download(url, destfile = "testfile.gz")
# h = new_handle(dirlistonly=TRUE)
# con = curl(url, "r", h)
