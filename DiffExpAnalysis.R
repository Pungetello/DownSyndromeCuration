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


#create a volcano plot of the data
volcano_plot = function(graph_data, output_prefix, file){
  
  #read in gene metadata to determine symbol and chromosome for each
  genes = read_tsv(paste0(getwd(), "/Data/Metadata/GeneMetadata/", output_prefix, ".tsv.gz"))
  # print(head(genes, n=100))
  graph_data = full_join(graph_data, select(rename(genes, gene = ensembl_gene_id), c("gene", "chromosome_name", "external_gene_name")), by = "gene")
  
  #print(head(graph_data))
  #print(tail(graph_data))
  # print(output_prefix)
  # print("HUMAN GENES IN GRAPH DATA:")
  # print(head(filter(graph_data, startsWith(gene, "ENSG"))))
  # print(sort(unique(pull(graph_data, "chromosome_name"))))
  
  
  if(grepl("MAC", file)){
    graph_data$chr21_flag = ifelse(graph_data$chromosome_name == "21", "Chr-21", 
                                   ifelse(graph_data$chromosome_name == "16", "Chr-16", 
                                          "Other"))
    # print("GENES ON CHR16:")
    # print(head(filter(graph_data, chromosome_name == "16")))
  }else{
    graph_data$chr21_flag = ifelse(graph_data$chromosome_name == "21",
                                   "Chr-21", "Other")
  }
  
  graph_data = graph_data[
    order(
      ifelse(graph_data$chromosome_name == "21", 3,
             ifelse(graph_data$chromosome_name == "16", 2, 1))
    ),
  ]
  
  #label top 5 of each category
  top = rbind(
    filter(graph_data, chr21_flag == "Chr-21")[order(graph_data$padj), ][1:5, ], 
    filter(graph_data, chr21_flag == "Chr-16")[order(graph_data$padj), ][1:5, ],
    filter(graph_data, chr21_flag == "Other")[order(graph_data$padj), ][1:5, ])
  #top = na.omit(top)
  print(top, width=Inf)
  
  ggplot(graph_data, aes(x = log2FoldChange, y = -log10(padj), color = chr21_flag)) +
    labs(color = "Chromosome") +
    geom_point(alpha = 0.2) +
    scale_color_manual(values = c("blue","black", "green")) +
    geom_text_repel(data = top, aes(label = external_gene_name), show.legend = FALSE) +
    theme_bw()
  
  ggsave(filename = paste0(getwd(), "/Data/Plots/", output_prefix, "_Volcano.png"), width = 10, height = 5, units = "in")
  
}


#----------Differential Expression Analysis-------------

#TODO: maybe reimplement this to itterate over Datasets list instead of the files in the folder, so other traits can be found
files = list.files(path = "Data/NormalizedData", pattern = "GSE[0-9]+\\w*_gene_counts\\.tsv")
#files = c("GSE190053_gene_counts.tsv")
#files = c("GSE101942_old.tsv.gz", "GSE190053_old.tsv.gz")
print(files)

file_location = "Data/Plots/"
if (!dir.exists(file_location)){dir.create(file_location, recursive = TRUE)}

for (file in files){
  
  #get gene_counts for the GSE
  counts = read_tsv(paste0("Data/NormalizedData/",file))
  #
  print("HUMAN GENES IN COUNTS:")
  counts%>%
    filter(startsWith(gene_id, "ENSG"))%>%
    arrange(across(2), decreasing = TRUE)%>%
    head()%>%
    print()


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
  # print(metadata)
  # print(head(counts))

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
  volcano_plot(sig_genes, gse, file)
  
}
  
  