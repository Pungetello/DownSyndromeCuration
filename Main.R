#function sources
source("PlatformsList.R")
source("DownloadFiles.R")

#installs packages that are necessary for later steps. These only need to be executed once
# source("InstallCorePackages.R")
# source("InstallArrayPackages.R")

#execute every time:
options(timeout = 10000) # allows very big tar files to still download
DownloadData(platforms_list, "Data/Files")

#creates Sample- and DatasetMetadata
source("StandardizeMetadata.R")
 
#creates NormalizedData for Affymetrix
source("QualityControlTests.R")

#reads NormalizedData to create GeneMetadata
source("GeneMetadata.R")

#downloads Data and GeneData for RNA
source("GetRNASecData.R")