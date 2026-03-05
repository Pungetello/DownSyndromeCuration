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
library(edgeR)

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
  #skip if already done
  if (!file.exists(paste0(getwd(),"/",index_file,".00.b.array"))){
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
  
  calculate_data_files(feature_counts)
}



#Use edgeR to calculate cpm and rpkm, calculate tpm manually
calculate_data_files = function(feature_counts){
  counts = feature_counts$counts
  gene_length = feature_counts$annotation$Length
  
  dge = DGEList(counts = counts, genes = data.frame(Length = gene_length))
  
  #calculate cpm
  cpm = cpm(dge)
  cpm_df = rownames_to_column(as.data.frame(cpm), "gene_id")
  colnames(cpm_df)[2] <- srr
  
  write_tsv(cpm_df, file=paste0(getwd(), "/Data/NormalizedData/", srr, "_CPM.tsv"))
  
  #calculate rpkm
  rpkm = rpkm(dge, gene.length = dge$genes$Length)
  
  rpkm_df = rownames_to_column(as.data.frame(rpkm), "gene_id")
  colnames(rpkm_df)[2] <- srr
  
  write_tsv(rpkm_df, file=paste0(getwd(), "/Data/NormalizedData/", srr, "_RPKM.tsv"))
  
  #calculate tpm
  length_kb <- gene_length / 1000
  rpk <- counts / length_kb
  tpm <- t( t(rpk) / colSums(rpk) ) * 1e6
  
  tpm_df = rownames_to_column(as.data.frame(tpm), "gene_id")
  colnames(tpm_df)[2] <- srr
  
  write_tsv(tpm_df, file=paste0(getwd(), "/Data/NormalizedData/", srr, "_TPM.tsv"))
}



#For each GSE in the GSE_to_SRR file, combine all the result data from its SRR's into one file
combine_results_per_GSE = function(){
  
  GSE_to_SRR = read_tsv(paste0(getwd(), "/Data/RNA_GSE_to_SRR.tsv"))
  
  gses = pull(GSE_to_SRR, GSE)%>%
    unique()
  
  #for(gse in gses){
    #print(gse)#debug
  gse = "GSE184771"
    combined_gene_counts = tibble(gene_id = character(), count = numeric())
    srrs = filter(GSE_to_SRR, GSE == gse)%>%
      pull(SRR)
    
    
    #combine_files(gse, srrs, "_gene_counts.csv")
    combine_files(gse, srrs, "_TPM.tsv")
    combine_files(gse, srrs, "_CPM.tsv")
    combine_files(gse, srrs, "_RPKM.tsv")
    
  #}
}



#takes in a list of file paths, reads them in, full joins them, and writes them to the out file path.
combine_files = function(gse, srrs, suffix){
  infiles = paste0(getwd(), "/Data/NormalizedData/", srrs, suffix)
  if(file.exists(infiles[1])){
    outfile = paste0(getwd(), "/Data/NormalizedData/", gse, suffix)
    combined_tibble = read_tsv(infiles[1])
    
    for (file in infiles[-1]){
      if(file.exists(file)){
        file_tibble = read_tsv(file)
        print(file_tibble)
        combined_tibble = full_join(combined_tibble, file_tibble, by = "gene_id")
      } else{
        print("FILE MISSING!")
        print(file)
      }
    }
    #remove .num at the end of gene_id's
    combined_tibble = mutate(combined_tibble, gene_id = str_remove(gene_id, "\\..*"))
    print(combined_tibble)#debug
    write_tsv(combined_tibble, outfile)
  }
}



#--------------process_RNA_data-------------

#get list of all srr files prefetched by previous script
srrs = list.files("Data/RawRNA")
print(srrs)

#for (srr in srrs){

srr = srrs[1]#debug
  print(srr)

  ref = paste0(getwd(), "/RefGenomes/GRCm39_ref.fna.gz")
  annotation = paste0(getwd(), "/RefGenomes/M38_ann.gtf.gz")
  index = "GRCm39_index"

  #finish installation by converting to fastq format
  install_raw(srr)

  build_index(index, ref)

  feature_counts = process_data(srr, index, annotation)

#}

combine_results_per_GSE()

