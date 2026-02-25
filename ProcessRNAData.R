#change library path for supercomputer
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "" || file.access(user_lib, 2) != 0) {
  user_lib <- "~/R/library"
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
}


#-----------loading_libraries-----------
library("Rsubread")
library("tidyverse")
library("dplyr")

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



#checks if index has already been built for given reference genome, builds it if not
build_index = function(index_file, ref){
  if (!file.exists(paste0(getwd(),"/",index_file,".00.b.array"))){ #add more to check for all of them?
    print(paste0(getwd(),"/", index_file,".00.b.array"))
    print(paste0("INDEX FILES NOT FOUND FOR ", index_file, ", CREATING"))
    #build the index from reference genome file
    buildindex(basename=index_file,reference=ref)#puts index file in current directory, maybe move?
  }
}



#runs alignment, gets feature counts, and makes tpm file
process_data = function(srr, index, annotation){
  
  output_file = paste0(getwd(), "/", srr, "_AlignResults.BAM") #test, move to /Data/NormalizedData eventually
  
  #check if alignment has already been done
  if (file.exists(paste0(output_file, ".summary"))){
    print("OUTPUT ALREADY EXISTS, SKIPPING")
    return()
  }
  
  #designate the files
  input_file = paste0(getwd(), "/fastq/", srr, "_1.fastq")
  input_file_2 = paste0(getwd(), "/fastq/", srr, "_2.fastq")
  
  
  #map to reference genome and get feature counts
  if(file.exists(input_file_2)){
    #paired end
    align(index=index,readfile1=input_file,readfile2=input_file_2,output_file=output_file,phredOffset=33)
    feature_counts = featureCounts(files=output_file, annot.ext=annotation, isGTFAnnotationFile = TRUE, isPairedEnd = TRUE)
    
  }else{
    #not paired end
    align(index=index,readfile1=input_file,output_file=output_file,phredOffset=33)
    feature_counts = featureCounts(files=output_file, annot.ext=annotation, isGTFAnnotationFile = TRUE)
  }
  #delete BAM files
  file.remove(output_file)
  file.remove(paste0(output_file, ".indel.vcf"))
  
  #save feature counts
  counts_df = as.data.frame(feature_counts$counts)
  counts_df$gene_id = rownames(counts_df)
  counts_df = rename(counts_df, srr = paste0(srr,"_AlignResults.BAM")) %>%
    select("gene_id", srr)
  print(counts_df)
  
  write_tsv(counts_df, file=paste0(getwd(), "/Data/NormalizedData/", srr, "_gene_counts.csv"))
  
  #calculate tpm file
  counts <- feature_counts$counts
  gene_length <- feature_counts$annotation$Length
  
  length_kb <- gene_length / 1000
  rpk <- counts / length_kb
  tpm <- t( t(rpk) / colSums(rpk) ) * 1e6
  
  tpm_df = as.data.frame(tpm)
  tpm_df$gene_id = rownames(tpm_df)
  tpm_df = rename(tpm_df, srr = paste0(srr,"_AlignResults.BAM")) %>%
    select("gene_id", srr)
  print(tpm_df)
  
  write_tsv(tpm_df, file=paste0(getwd(), "/Data/NormalizedData/", srr, "_TPM.txt"))
  
}



#--------------process_RNA_data-------------

#get list of all srr files prefetched by previous script
srrs = list.files("Data/RawRNA")
print(srrs)

#finish installation by converting to fastq format
#for (srr in srrs){
  print(srrs[1])
  
  ref = paste0(getwd(), "/RefGenomes/GRCm39_ref.fna.gz")
  annotation = paste0(getwd(), "/RefGenomes/M38_ann.gtf.gz")
  index = "GRCm39_index"
  
  install_raw(srrs[1])
  
  build_index(index, ref)
  
  process_data(srrs[1], index, annotation)
  
#}





