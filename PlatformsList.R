platforms_list = list(
  "GSE110064" = "GPL570",
  "GSE48611"= "GPL570",
  "GSE1789" = "GPL96"
)

#Populated by InstallArrayPackages
platforms_to_package_list = list()    #TODO: decide if we should save this somewhere somehow do it doesn't get lost once we finish running it

#----------functions-------------

create_unique_vector = function(platforms_list){
  # returns a vector of one key for every unique value in the list
  
  platforms_vector = unlist(platforms_list)
  value_to_keys = split(names(platforms_vector), platforms_vector)
  keys_vector = unname(unlist(lapply(value_to_keys, `[`, 1)))
  return(keys_vector)
}