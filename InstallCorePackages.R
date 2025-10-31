#-----------installing_libraries-----------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21")
BiocManager::install(c("tidyverse",
                       "arrayQualityMetrics",
                       "BiocManager",
                       "SCAN.UPC",
                       "GEOquery",
                       "doParallel",
                       "arrayQualityMetrics",
                       "affy",
                       "pdInfoBuilder", 
                       "affxparser", 
                       "Biostrings",
                       "janitor"
                       )
)

#-----------libraries_likely_needed_for_future_step-----------
#"clariomshumanhsentrezgprobe",
# "pd.clariom.s.human",
# "clariomshumancdf"
