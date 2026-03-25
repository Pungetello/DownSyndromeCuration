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
library(ggrepel)

#----------functions-------------


#creates a metadata file for deseq2 for the GSE given
create_metadata = function(gse, column_type){
  sample_metadata = read_tsv("Data/Metadata/SampleMetadata.tsv")
  
  metadata = filter(sample_metadata, Dataset_ID == gse)%>%
    print()%>%#debug
    select(ID, Value)%>%
    rename(GSM = ID)
  
  if(column_type=="srr"){
    GSE_to_SRR = read_tsv("Data/RNA_GSE_to_SRR.tsv")
    metadata = full_join(metadata, filter(GSE_to_SRR, GSE==gse), by = "GSM")%>%
      select(SRR, Value)%>%
      as.data.frame()
    
    #format correctly for deseq2
    rownames(metadata) = metadata$SRR
    metadata$SRR = NULL
    metadata$Value = factor(metadata$Value)
  }else{
    metadata = select(metadata, GSM, Value)%>%
      as.data.frame()%>%
      print()#debug
    
    #format correctly for deseq2
    rownames(metadata) = metadata$GSM
    metadata$GSM = NULL
    metadata$Value = factor(metadata$Value)
  }
  
  return(metadata)
}


#create a volcano plot of the data
volcano_plot = function(data, output_prefix){
  
  top = data[order(data$padj), ][1:10, ]
  
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj)))+ #, color=significant)) +
    geom_point(alpha = 0.5) +
    #scale_color_manual(values = c("grey", "blue")) +
    geom_text_repel(data = top, aes(label = gene)) +
    theme_minimal()
  
  ggsave(filename = paste0(getwd(), "/Data/Plots/", output_prefix, "_Volcano.png"), width = 5, height = 5, units = "in")
  
  
}


#----------Differential Expression Analysis-------------

files = list.files(path = "Data/NormalizedData", pattern = "GSE[0-9]+_gene_counts\\.csv")
#files = c("GSE101942.tsv.gz", "GSE190053.tsv.gz")

file_location = "Data/Plots/"
if (!dir.exists(file_location)){dir.create(file_location, recursive = TRUE)}

#file = "Data/NormalizedData/GSE109294_gene_counts.csv" #debug
#for (file in files){
file = files[1]
  print(file) #debug
  #get gene_counts for the GSE
  counts = read_tsv(paste0("Data/NormalizedData/",file))
  counts = as.data.frame(counts)
  rownames(counts) = counts$gene_id
  counts$gene_id = NULL
  counts = as.matrix(counts)
  print(head(counts))#debug
  
  #create tibble mapping each sample to 'control_group' or 'affected_group'
  gse = strsplit(basename(file), "_")[[1]][1]
  metadata = create_metadata(gse, "srr")
  
  print(metadata)#debug
  if(length(unique(metadata$Value)) < 2){
    print("Only one variable, skipping dataset")
    next()
  }
  
  #set up input
  dds = DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ Value
  )
  dds$Value = as.factor(metadata$Value)
  dds$Value = relevel(dds$Value, ref="control_group")
  
  #run the analysis
  dds = DESeq(dds)
  
  results = results(dds)
  
  #filter results to adjusted p-value < 0.05 and sort by padj
  #sig_genes = subset(results, padj < 0.05)%>%
  sig_genes = as_tibble(results, rownames='gene')%>%
    drop_na()%>%
    arrange(padj)
  
  #write to file
  write_tsv(sig_genes, file=paste0(getwd(), "/Data/NormalizedData/", gse, "_DE.tsv"))
  
  #create and save volcano plot
  volcano_plot(sig_genes, gse)
  
#}
  
  