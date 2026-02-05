#replace with Conda?
#check if can do it one by one

#-----------installing_libraries-----------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")
BiocManager::install(c("tidyverse",
                       "SCAN.UPC",
                       "GEOquery",
                       "arrayQualityMetrics",
                       "janitor",
                       "biomaRt",
                       "Rcurl",
                       "rentrez",
                       "Rsubread"
                       )
)

#-----------Previously_included_libraries-----------
#in case we need to go back to a previous implementation

#"pdInfoBuilder", 
#"affxparser", 
#"Biostrings",
#"doParallel",
#"osfr",
#"affy",
#"BiocManager"

#TODO: go through this periodically to see if there are any we are no longer using
