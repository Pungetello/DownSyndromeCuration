#installs packages that are necessary for later steps. These only need to be executed once
source("InstallCorePackages.R")
source("DownloadBrainarrayExamples.R")
source("InstallArrayPackages.R")  #TODO: finish this once website is back up

#execute every time:
source("StandardizeMetadata.R")
source("RunStatisticalTests.R")
