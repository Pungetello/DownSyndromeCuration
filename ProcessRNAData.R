#-----------loading_libraries-----------
source("PlatformsList.R")
#source ~/.bashrc


#----------functions-------------

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



#get the true raw data using the SRA toolkit
download_raw = function(geo_id){
  
  gse = getGEO(geo_id, GSEMatrix = FALSE)
  srrs = Meta(gse)$relation
  print(ssr)
  
  system(paste(
    "fasterq-dump",
    srr,
    "--split-files",
    "--outdir fastq"
  ))
  
}



#--------------process_RNA_data-------------

#filter to geo_ids for RNAsec that do not have NormalizedData downloaded. Make sure to run GetRNASecData before this.
for (geo_id in names(platforms_list)){
  platform = platforms_list[[geo_id]]
  
  if(is.na(platform)){
    destination = paste0(getwd(), "/Data/NormalizedData/", geo_id, ".tsv.gz")
    if(!file.exists(destination)){
      
      #make sure SRA toolkit has been downloaded
      check_sra()
      
      #download the raw data
      download_raw(geo_id)
      
      #process the data
      print(geo_id)
      
      
      
    }
  }
}