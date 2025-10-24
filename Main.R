#installs packages that are necessary for later steps. These only need to be executed once
source("PlatformsList.R")
# source("InstallCorePackages.R")

#Download Brainarray Examples
source("DownloadTarFiles.R")
DownloadData(create_unique_vector(platforms_list), "BrainArrayExamples")

# source("InstallArrayPackages.R")  #TODO: finish this once website is back up

source("PlatformsList.R")
#execute every time:
DownloadData(platforms_list, "Data/Files")
source("StandardizeMetadata.R")
# source("RunStatisticalTests.R")
