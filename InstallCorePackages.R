#-----------installing_libraries-----------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21")
BiocManager::install(c("tidyverse",
                       "BiocManager",#
                       "SCAN.UPC",
                       "GEOquery",
                       "arrayQualityMetrics",
                       "affy",#
                       "janitor",
                       "biomaRt",
                       "Rcurl"
                       )
)

#-----------Previously_included_libraries-----------

#"pdInfoBuilder", 
#"affxparser", 
#"Biostrings",
#"doParallel",
#"osfr",

#TODO: go through this periodically to see if there are any we are no longer using
