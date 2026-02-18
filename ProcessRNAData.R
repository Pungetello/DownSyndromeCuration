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
  out_path = paste0(getwd(), "/fastq/", srr, "_1.fastq")
  
  #skip if file has already been downloaded
  if (!file.exists(out_path)){
    system2(
      fasterq,
      args = c(in_path, "--split-files", "--outdir fastq"))
  }
}



process_data = function(srr){
  input_file = paste0(getwd(), "/fastq/", srr, "_1.fastq")
  ref = paste0("strain_name.fa")
  annotation = "strain_name.gtf"
  
  buildindex(basename="strain_name_index",reference=ref)#puts index file in current directory
  featureCounts(files="*.bam", annot.ext=annotation)
  
  align.stat <- align(index="reference_index",readfile1=file_location,output_file="alignResults.BAM",phredOffset=64)
  
}



#--------------process_RNA_data-------------

#get list of all srr files prefetched by previous script
srrs = list.files("Data/RawRNA")
print(srrs)

#finish installation by converting to fastq format
for (srr in srrs){
  print(srr)
  install_raw(srr)
  
  #process_data(srr)
  
}





