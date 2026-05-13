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
  
  #do stuff
  #learn the package, identify latent covariates, save to file, compare against metadata variables.
  
}



#---------Find Surrogate Variables----------

files = list.files(path = "Data/NormalizedData", pattern = "GSE[0-9]+\\w*_gene_counts\\.tsv")

for(file in files){
  
  path = paste0(getwd(), "/Data/NormalizedData/", file)
  gse = strsplit(basename(file), "_")[[1]][1]
  
  run_sva(path, gse)
  
}
