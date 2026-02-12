#change library path for supercomputer
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "" || file.access(user_lib, 2) != 0) {
  user_lib <- "~/R/library"
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
}


#-----------loading_libraries-----------


#----------functions-------------



#command locations

#Optimus Prime locations
# fasterq = normalizePath(paste0(getwd(),"/sratoolkit.current-win64/sratoolkit.3.3.0-win64/bin/fasterq-dump.exe"))
# prefetch = normalizePath(paste0(getwd(),"/sratoolkit.current-win64/sratoolkit.3.3.0-win64/bin/prefetch.exe"))

#Supercomputer locations
fasterq = "fasterq-dump"
prefetch = "prefetch"


#install the raw data using the SRA toolkit
install_raw = function(srr){

  system2(
    fasterq,
    args = c(srr, "--split-files", "--outdir fastq"))
    #"-p"
    #wait = TRUE, stdout = TRUE)
  
}



#--------------process_RNA_data-------------

#filter to geo_ids for RNAsec that do not have NormalizedData downloaded. Make sure to run GetRNASecData before this.
# for (geo_id in names(platforms_list)){
#   platform = platforms_list[[geo_id]]
#   
#   if(is.na(platform)){
#     destination = paste0(getwd(), "/Data/NormalizedData/", geo_id, ".tsv.gz")
#     if(!file.exists(destination)){
#       
#       #process the data
#       print(geo_id)
#       install_raw(geo_id)
#       
#       #put mouse strain list in platformslist
#       #download the tables for strains needed
#       
#       
#       
#     }
#   }
# }

srrs = list.files("Data/RawRNA")
print(srrs)

for (srr in srrs){
  print(srr)
  install_raw(srr)
  
}