#change library path for supercomputer
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "" || file.access(user_lib, 2) != 0) {
  user_lib <- "~/R/library"
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
}


#-----------loading_libraries-----------

library(tidyverse)
library(sva)
library(janitor)
library(ggplot2)

#----------functions-------------

#use surogate variable analysis to find surrogate variables
run_sva = function(path, gse){
  gene_counts = read_tsv(path)
  
  #create full model tibble
  GSM_to_Value = read_tsv(paste0(getwd(), "/Data/Metadata/SampleMetadata.tsv"))%>%
    filter(Dataset_ID == gse)%>%
    dplyr::rename(GSM = ID)%>%
    select(GSM, Value)
  var_of_interest = read_tsv(paste0(getwd(), "/Data/RNA_GSE_to_SRR.tsv"))%>%
    filter(GSE == gse)%>%
    inner_join(GSM_to_Value, by = "GSM")%>%
    dplyr::rename(status = Value)%>% #unsure what to name this. Karyotype? Status? Control_v_Affected? Group?
    select(SRR, status)
  
  #make sure SRRs are in the same order
  var_of_interest = arrange(var_of_interest, SRR)
  gene_counts = select(gene_counts, c(gene_id, pull(var_of_interest, SRR)))
  
  #convert full model to matrix
  mod = model.matrix(~status, data = var_of_interest)
  
  #null model matrix is just indeces, since we are not including any other variables
  mod0 = model.matrix(~1, data = var_of_interest)
  
  #convert gene_counts to matrix, filter out rows with no variance
  gene_counts_matrix = column_to_rownames(gene_counts, var = "gene_id") %>% 
    as.matrix()
  
  variance = apply(gene_counts_matrix, 1, var)
  
  gene_counts_matrix = gene_counts_matrix[variance > 0, ]
  
  #determine how many surrogate variables to look for
  n.sv = 2
  
  #run sva
  svobj = sva(gene_counts_matrix, mod, mod0, n.sv = n.sv)
  
  sv = as.data.frame(svobj$sv)
  sv = cbind(SRR = colnames(gene_counts_matrix), sv)
  write_tsv(sv, paste0(getwd(), "/Data/SVAResults/", gse, "_sva.tsv"))
  
  return(sv)
}



# retrieves metadata using GEOquery function
get_metadata = function(geo_ID) {
  metadata = getGEO(geo_ID)[[1]] 
  
  metadata = as_tibble(pData(metadata))
  metadata = clean_names(metadata)
  #metadata = fix_bespoke_issues(geo_ID, metadata)
  return (metadata)
}



#make a graph comparing the sva values to existing variables
make_graph = function(gse, sv){
  if(gse == "GSE109293"||gse == "GSE109293"){
    #No other variables besides the test variable, unsure what to do
    return()
  }else if(gse == "GSE202938"){
    metadata = get_metadata(gse)%>%
      dplyr::rename("GSM" = "geo_accession")
    GSE_to_SRR = read_tsv(paste0(getwd(), "/Data/RNA_GSE_to_SRR.tsv"))
    metadata = inner_join(metadata, GSE_to_SRR, by = "GSM")%>%
      select(c("SRR","GSM","genotype_ch1"))
    
    sv = inner_join(sv, metadata, by = "SRR")%>%
      select(c(SRR, V1, V2, genotype_ch1))
    
    plot_data = pivot_longer(sv, cols = starts_with("V"),
                             names_to = "SV",
                             values_to = "SV_value")%>%
      pivot_longer(cols = c("genotype_ch1"),
                   names_to = "variable",
                   values_to = "group")
    
    #graph it!
    print("GRAPH!")
      ggplot(plot_data, aes(x = factor(group), y = SV_value)) +
        geom_boxplot() +
        facet_grid(SV ~ variable, scales = "free_x") +
        theme_bw()
    
      ggsave(filename = paste0(getwd(), "/Data/SVAResults/", gse, "_plots.png"), width = 10, height = 5, units = "in")
    
  }else if(gse == "GSE210117"){
    metadata = get_metadata(gse)%>%
      dplyr::rename("GSM" = "geo_accession")%>%
      print()
    GSE_to_SRR = read_tsv(paste0(getwd(), "/Data/RNA_GSE_to_SRR.tsv"))
    metadata = inner_join(metadata, GSE_to_SRR, by = "GSM")%>%
      select(c("SRR","GSM","genotype_ch1", "age_ch1","sex_ch1","tissue_ch1"))
    sv = inner_join(sv, metadata, by = "SRR")%>%
      select(c(SRR, V1, V2, genotype_ch1, age_ch1,sex_ch1,tissue_ch1))
    
    plot_data = pivot_longer(sv, cols = starts_with("V"),
                             names_to = "SV",
                             values_to = "SV_value")%>%
      pivot_longer(cols = c("genotype_ch1","age_ch1","sex_ch1","tissue_ch1"),
                   names_to = "variable",
                   values_to = "group")
    
    #graph it!
    print("GRAPH!")
      ggplot(plot_data, aes(x = factor(group), y = SV_value)) +
        geom_boxplot() +
        facet_grid(SV ~ variable, scales = "free_x") +
        theme_bw()
      
      ggsave(filename = paste0(getwd(), "/Data/SVAResults/", gse, "_plots.png"), width = 10, height = 5, units = "in")
  }
  
  #print(sv)
  
}



#---------Find Surrogate Variables----------

file_location = "Data/SVAResults/"
if (!dir.exists(file_location)){dir.create(file_location, recursive = TRUE)}

files = list.files(path = "Data/NormalizedData", pattern = "GSE[0-9]+\\w*_gene_counts\\.tsv")
print(files)
for(file in files){
  
  path = paste0(getwd(), "/Data/NormalizedData/", file)
  gse = strsplit(basename(file), "_")[[1]][1]
  
  sv = run_sva(path, gse)
  
  make_graph(gse, sv)
  
}
