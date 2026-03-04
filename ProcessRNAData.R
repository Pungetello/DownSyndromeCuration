#change library path for supercomputer
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "" || file.access(user_lib, 2) != 0) {
  user_lib <- "~/R/library"
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
}


#-----------loading_libraries-----------
library(Rsubread)
library(tidyverse)
library(dplyr)

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
  
  output_file = paste0(srr, "_AlignResults.BAM") #test, move to /Data/NormalizedData eventually
  
  #check if alignment has already been done
  if (file.exists(paste0(getwd(), "/Data/NormalizedData/", srr, "_TPM.txt"))){
    print("OUTPUT ALREADY EXISTS, SKIPPING")
    return()
  }
  
  #designate the files
  input_file = paste0(getwd(), "/fastq/", srr, "_1.fastq")
  input_file_2 = paste0(getwd(), "/fastq/", srr, "_2.fastq")
  
  #make sure data is downloaded
  if (!file.exists(input_file)){
    print(paste0("INPUT FILE NOT FOUND FOR ", srr))
    return()
  }
  
  
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
  counts_df = rename(counts_df, !!srr := paste0(srr,"_AlignResults.BAM")) %>%
    select("gene_id", srr)
  
  write_tsv(counts_df, file=paste0(getwd(), "/Data/NormalizedData/", srr, "_gene_counts.csv"))
  
  #calculate tpm file
  counts <- feature_counts$counts
  gene_length <- feature_counts$annotation$Length
  
  length_kb <- gene_length / 1000
  rpk <- counts / length_kb
  tpm <- t( t(rpk) / colSums(rpk) ) * 1e6
  
  tpm_df = as.data.frame(tpm)
  tpm_df$gene_id = rownames(tpm_df)
  tpm_df = rename(tpm_df, !!srr := paste0(srr,"_AlignResults.BAM")) %>%
    select("gene_id", srr)
  
  write_tsv(tpm_df, file=paste0(getwd(), "/Data/NormalizedData/", srr, "_TPM.txt"))
  
}



combine_results_per_GSE = function(){
  
  GSE_to_SRR = read_tsv(paste0(getwd(), "/Data/RNA_GSE_to_SRR.tsv"))
  
  gses = pull(GSE_to_SRR, GSE)%>%
    unique()
  
  print(gses)#debug
  
  for(gse in gses){
    combined_gene_counts = tibble(gene_id = character(), count = numeric())
    srrs = filter(GSE_to_SRR, GSE == gse)%>%
      pull(SRR)
    print(srrs)#debug
    
    gene_count_files = paste0(getwd(), "Data/NormalizedData/", srrs, "_gene_counts.csv")
    TPM_files = TPM = paste0(getwd(), "Data/NormalizedData/", srrs, "_TMP.txt")
    
    gene_counts_filename = paste0(getwd(), "Data/NormalizedData/", gse, "_gene_counts.tsv")
    combine_files(gene_count_files, gene_counts_filename)
    
    TPM_filename = paste0(getwd(), "Data/NormalizedData/", gse, "_TMP.tsv")
    combine_files(gene_count_files, gene_counts_filename)
  }
  
}


combine_files = function(infiles, outfile){
  combined_tibble = tibble()
  for (file in infiles){
    if(file.exists(file)){
      file_tibble = read_tsv(file)
      combined_tibble <<- full_join(combined_tibble, file_tibble)
    }
  }
  print(combined_tibble)#debug
  
  write_tsv(combined_tibble, outfile)
}



#--------------process_RNA_data-------------

#get list of all srr files prefetched by previous script
srrs = list.files("Data/RawRNA")
#print(srrs)

#finish installation by converting to fastq format
# for (srr in srrs){
#   print(srr)
#   
#   ref = paste0(getwd(), "/RefGenomes/GRCm39_ref.fna.gz")
#   annotation = paste0(getwd(), "/RefGenomes/M38_ann.gtf.gz")
#   index = "GRCm39_index"
#   
#   install_raw(srr)
#   
#   build_index(index, ref)
#   
#   process_data(srr, index, annotation)
#   
# }

combine_results_per_GSE()





