#packages

#installs packages that are necessary for later steps. These two scripts only need to be executed once
#source("install_core_packages.R")
#source("install_array_packages.R")

#execute every time:
#source(LibraryInstaller.R)
source("Downloader.R")
#source(Metadata_fixer)
#source(Statistical_tester)
#source(platforms_list)

platforms_list = list(
  "GSE110064" = "GPL570",
  "GSE48611"= "GPL570",
  "GSE1789" = "GPL96"
)

#for (i in length(names(platforms_list))){
  download_data((names(platforms_list)), "C:/Users/punge/DownSyndromeCuration/Platforms")
#}
