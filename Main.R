#function sources
source("PlatformsList.R")
source("DownloadFiles.R")

#installs packages that are necessary for later steps. These only need to be executed once
# source("InstallCorePackages.R")
# source("InstallArrayPackages.R")

#execute every time:
options(timeout = 10000) # allows very big tar files to still download
DownloadData(platforms_list, "Data/Files") #TODO: maybe change this to source, if we never use this function again?

#creates Sample- and DatasetMetadata
source("StandardizeMetadata.R")
 
#creates NormalizedData for Affymetrix
source("QualityControlTests.R")

#downloads gene_data for RNA
source("GetRNASecData.R")

#processes RNA data that is not downloadable
source("ProcessRNAData.R")

#creates GeneMetadata by reading NormalizedData for affymetrix, gene_data for RNA.
source("GeneMetadata.R")