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
}



#---------Find Surrogate Variables----------

file_location = "Data/SVAResults/"
if (!dir.exists(file_location)){dir.create(file_location, recursive = TRUE)}

files = list.files(path = "Data/NormalizedData", pattern = "GSE[0-9]+\\w*_gene_counts\\.tsv")
print(files)
for(file in files){
  
  path = paste0(getwd(), "/Data/NormalizedData/", file)
  gse = strsplit(basename(file), "_")[[1]][1]
  
  run_sva(path, gse)
  
}
