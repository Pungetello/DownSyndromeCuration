#-----------loading_libraries-----------
source("PlatformsList.R")


#-----------functions------------------
quality_control_removal(
  #do the three statistical tests -- remove array data that fails all three
)

#line 286 uses SCAN with the name of the Brainarray package to normalize the data

#-----------run statistical tests-----------


filtered_file_list = quality_control_removal(file_list, platform, geo_id)