#function sources
source("PlatformsList.R")
source("DownloadFiles.R")

#installs packages that are necessary for later steps. These only need to be executed once
# source("InstallCorePackages.R")
##Download Brainarray Examples
DownloadData(create_unique_vector(platforms_list), "BrainArrayExamples")
# source("InstallArrayPackages.R")  #TODO: finish this once website is back up           #TODO: ask if we should just download all the data at the beginning


#execute every time:
DownloadData(names(platforms_list), "Data/Files")
source("StandardizeMetadata.R")
# source("RunStatisticalTests.R")

#TODO: add col headers or something.
