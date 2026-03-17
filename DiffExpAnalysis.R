#change library path for supercomputer
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "" || file.access(user_lib, 2) != 0) {
  user_lib <- "~/R/library"
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
}


#-----------loading_libraries-----------

library(DESeq2)
library(tidyverse)


#----------functions-------------


#creates a metadata file for deseq2 for the GSE given
create_metadata = function(gse){
  GSE_to_SRR = read_tsv("Data/RNA_GSE_to_SRR.tsv")
  sample_metadata = read_tsv("Data/Metadata/SampleMetadata.tsv")
  
  metadata = filter(sample_metadata, Dataset_ID == gse)%>%
    select(ID, Value)%>%
    rename(GSM = ID)%>%
    full_join(filter(GSE_to_SRR, GSE==gse), by = "GSM")%>%
    select(SRR, Value)%>%
    as.data.frame()
  
  rownames(metadata) = metadata$SRR
  metadata$SRR <- NULL
  
  metadata$Value <- factor(metadata$Value)
  
  return(metadata)
}


#----------Differential Expression Analysis-------------

files = list.files(path = "/Data/NormalizedData/", pattern = "GSE[0-9]+_gene_counts\\.csv")
print("FILES:")
print(files)#debug

file = "Data/NormalizedData/GSE202938_gene_counts.csv" #debug
#for (file in files){
  #get gene_counts for the GRE
  counts = read_tsv(file)
  counts = as.data.frame(counts)
  rownames(counts) = counts$gene_id
  counts$gene_id = NULL
  counts = as.matrix(counts)
  print(head(counts, 10)) #debug
  
  #create tibble mapping each sample to 'control_group' or 'affected_group'
  gse = strsplit(basename(file), "_")[[1]][1]
  metadata = create_metadata(gse)
  print(metadata)#debug
  
  
  #set up input
  dds = DESeqDataSetFromMatrix( # MOST RECENT ERROR: Error in relevel.default(dds$condition, ref = "control_group") :
    #'relevel' only for (unordered) factors
    #Calls: relevel -> relevel.default
    countData = counts,
    colData = metadata,
    design = ~ Value
  )
  dds$Value = as.factor(metadata$Value)
  dds$Value = relevel(dds$Value, ref="control_group")
  
  #run the analysis
  dds = DESeq(dds)
  
  results = results(dds)
  
  #filter results to adjusted p-value < 0.05
  sig_genes = subset(results, padj < 0.05)%>%
    arrange(padj)
  
  #write to file
  write_tsv(sig_genes, file=paste0(getwd(), "/Data/NormalizedData/", gse, "_DE.tsv"))
  
#}







