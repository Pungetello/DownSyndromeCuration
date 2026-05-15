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

#----------functions-------------

#use surogate variable analysis to find surrogate variables
run_sva = function(path, gse){
  gene_counts = read_tsv(path)
  
  #create full model tibble
  GSM_to_Value = read_tsv(paste0(getwd(), "/Data/Metadata/SampleMetadata.tsv"))%>%
    filter(Dataset_ID == gse)%>%
    rename(GSM = ID)%>%
    select(GSM, Value)
  var_of_interest = read_tsv(paste0(getwd(), "/Data/RNA_GSE_to_SRR.tsv"))%>%
    filter(GSE == gse)%>%
    inner_join(GSM_to_Value, by = "GSM")%>%
    rename(status = Value)%>% #unsure what to name this. Karyotype? Status? Control_v_Affected? Group?
    #mutate(status)
    select(SRR, status)
  
  #make sure SRRs are in the same order
  var_of_interest = arrange(var_of_interest, SRR)
  gene_counts = select(gene_counts, c(gene_id, pull(var_of_interest, SRR)))
  
  #convert full model to matrix
  mod = model.matrix(~status, data = var_of_interest)
  
  #null model matrix is just indeces, since we are not including any other variables
  mod0 = model.matrix(~1, data = var_of_interest)
  
  #determine how many surrogate variables to look for
  n.sv = 2
  
  gene_counts_matrix = column_to_rownames(gene_counts, var = "gene_id") %>% 
    as.matrix()
  
  print(mod)
  print(mod0)
  print(head(gene_counts_matrix))
  
  print(class(gene_counts_matrix))
  print(typeof(gene_counts_matrix))
  
  print(class(mod))
  print(typeof(mod))
  
  vars <- apply(gene_counts_matrix, 1, var)
  print(sum(vars == 0))
  
  gene_counts_matrix = gene_counts_matrix[vars > 0, ]
  print(head(gene_counts_matrix))
  
  #run sva
  svobj = sva(gene_counts_matrix, mod, mod0, n.sv = n.sv)%>%
    print()
  

  stop()
  
  
  #do stuff
  #learn the package, identify latent covariates, save to file, compare against metadata variables.
  
}



#---------Find Surrogate Variables----------

files = list.files(path = "Data/NormalizedData", pattern = "GSE[0-9]+\\w*_gene_counts\\.tsv")
print(files)
for(file in files){
  
  path = paste0(getwd(), "/Data/NormalizedData/", file)
  gse = strsplit(basename(file), "_")[[1]][1]
  
  run_sva(path, gse)
  
}
