#packages
source(PlatformsList)

#installs packages that are necessary for later steps. These two scripts only need to be executed once
source("InstallCorePackages.R")
source("InstallArrayPackages.R")

#execute every time:
#source(LibraryInstaller.R)
source("Downloader.R")
source("ColumnCleaner.R")
#source(Statistical_tester)





