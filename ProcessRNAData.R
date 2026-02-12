#change library path for supercomputer
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "" || file.access(user_lib, 2) != 0) {
  user_lib <- "~/R/library"
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
}


#-----------loading_libraries-----------
library("Rsubread")

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
  in_path = paste0(getwd(), "/Data/RawRNA/", srr, "/", srr, ".sra")

  system2(
    fasterq,
    args = c(srr, in_path, "--split-files", "--outdir fastq"))
    #"-p"
    #wait = TRUE, stdout = TRUE)
  
}



process_data = function(srr){
  file_location = paste0(getwd(), "/fastq/", srr)#just guessing for now, fix later
  
  ref <- system.file("extdata","reference.fa",package="Rsubread")#default reference. Is probably human, so need different one.
  buildindex(basename="reference_index",reference=ref)#puts index file in current directory
  
  reads <- system.file("extdata",file_location,package="Rsubread")#maps read dataset to reference
  align.stat <- align(index="reference_index",readfile1=reads,output_file="alignResults.BAM",phredOffset=64)
  
}



#--------------process_RNA_data-------------

#get list of all srr files prefetched by previous script
srrs = list.files("Data/RawRNA")
print(srrs)

#finish installation by converting to fastq format
for (srr in srrs){
  print(srr)
  install_raw(srr)
  
  process_data(srr)
  
}





