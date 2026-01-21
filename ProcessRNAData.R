#-----------loading_libraries-----------
source("PlatformsList.R")



#----------functions-------------



#--------------process_RNA_data-------------

for (geo_id in names(platforms_list)){
  platform = platforms_list[[geo_id]]
  
  if(is.na(platform)){
    destination = paste0(getwd(), "/Data/NormalizedData/", geo_id, ".tsv.gz")
    if(!file.exists(destination)){
      
      
      
      #process the data
      print(geo_id)
      
      
      
    }
  }
}