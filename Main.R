#function sources
source("PlatformsList.R")
source("DownloadFiles.R")

#installs packages that are necessary for later steps. These only need to be executed once
# source("InstallCorePackages.R")
# source("InstallArrayPackages.R")

#execute every time:
options(timeout = 10000) # allows very big tar files to still download
DownloadData(platforms_list, "Data/Files")
source("StandardizeMetadata.R")
source("QualityControlTests.R")
source("GeneMetadata")
