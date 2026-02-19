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
  #build the index from annotation and reference genome files
  ref = paste0(getwd(), "/RefGenomes/GRCm39_ref.fna.gz")
  annotation = paste0(getwd(), "/RefGenomes/GRCm39_ann.gtf.gz")
  
  buildindex(basename="GRCm39_index",reference=ref)#puts index file in current directory, maybe move
  
  #designate the files
  input_file = paste0(getwd(), "/fastq/", srr, "_1.fastq")
  input_file_2 = paste0(getwd(), "/fastq/", srr, "_2.fastq")
  
  output_file = paste0(getwd(), "/", srr, "AlignResults.BAM") #test, move to /Data/NormalizedData eventually
  
  #map to reference genome
  if(file.exists(input_file_2)){
    #paired end
    align.stat2 = align(index="GRCm39_index",readfile1=input_file,readfile2=input_file_2,output_file=output_file,phredOffset=64)
    
  }else{
    #not paired end
    align.stat = align(index="GRCm39_index",readfile1=input_file,output_file=output_file,phredOffset=64)
  }
  
  #save feature counts
  feature_counts = featureCounts(files="*.BAM", annot.ext=annotation)
  write.csv(feature_counts$counts, file=paste0(srr, "gene_counts.csv")) #again, move eventually
  
}



#--------------process_RNA_data-------------

#get list of all srr files prefetched by previous script
srrs = list.files("Data/RawRNA")
print(srrs)

#finish installation by converting to fastq format
#for (srr in srrs){
  print(srr[1])
  install_raw(srrs[1])
  
  process_data(srrs[1])
  
#}





