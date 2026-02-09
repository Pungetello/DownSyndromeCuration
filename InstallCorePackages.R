#replace with Conda?
#check if can do it one by one

#-----------installing_libraries-----------

#check to see that the library path is writable; if not, change it.
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "" || file.access(user_lib, 2) != 0) {
  user_lib <- "~/R/library"
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
}

print(user_lib)



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22", lib = user_lib)
BiocManager::install(c("tidyverse",
                       "SCAN.UPC",
                       "GEOquery",
                       "arrayQualityMetrics",
                       "janitor",
                       "biomaRt",
                       "Rcurl",
                       "rentrez",
                       "Rsubread"
                       ), lib = user_lib
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
