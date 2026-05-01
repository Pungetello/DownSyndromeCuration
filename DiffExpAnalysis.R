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
      as.data.frame()
    
    #format correctly for deseq2
    rownames(metadata) = metadata$GSM
    metadata$GSM = NULL
    metadata$Value = factor(metadata$Value)
  }
  
  return(metadata)
}

#TODO: highlight genes on anomalous chormosome
#create a volcano plot of the data
volcano_plot = function(data, output_prefix){
  
  top = data[order(data$padj), ][1:10, ]
  
  genes = read_tsv(paste0(getwd(), "/Data/Metadata/GeneMetadata/", output_prefix, ".tsv.gz"))
  data = inner_join(data, select(rename(genes, gene = ensembl_gene_id), c("gene", "chromosome_name")), by = "gene")
  
  print(head(data))
  
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = factor(chromosome_name == "21",
                                                                        labels = c("Other", "Chr21")))) +
    theme(plot.background = element_rect(fill = "white"))+
    geom_point(alpha = 0.5) +
    #scale_color_manual(values = c("grey", "blue")) +
    geom_text_repel(data = top, aes(label = gene)) +
    theme_minimal()
  
  ggsave(filename = paste0(getwd(), "/Data/Plots/", output_prefix, "_Volcano.png"), width = 5, height = 5, units = "in")
  
}


#----------Differential Expression Analysis-------------

#TODO: maybe reimplement this to itterate over Datasets list instead of the files in the folder, so other traits can be found
#files = list.files(path = "Data/NormalizedData", pattern = "GSE[0-9]+\\w*_gene_counts\\.tsv")
files = c("GSE190053_gene_counts.csv")
#files = c("GSE101942_old.tsv.gz", "GSE190053_old.tsv.gz")
print(files)

file_location = "Data/Plots/"
if (!dir.exists(file_location)){dir.create(file_location, recursive = TRUE)}

for (file in files){
  #get gene_counts for the GSE
  counts = read_tsv(paste0("Data/NormalizedData/",file))
  counts = as.data.frame(counts)
  rownames(counts) = counts$gene_id
  counts$gene_id = NULL
  counts = as.matrix(counts)
  
  #create tibble mapping each sample to 'control_group' or 'affected_group'
  gse = strsplit(basename(file), "_")[[1]][1]#"\\." for human
  metadata = create_metadata(gse, "srr")#"gsm" for human
  
  if(length(unique(metadata$Value)) < 2){
    print("Only one variable, skipping dataset")
    next()
  }
  #print(metadata)
  #print(head(counts))
  
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
  
}
  
  