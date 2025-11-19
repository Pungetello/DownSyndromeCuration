platforms_list_new = list(
  "GSE110064" = "hgu133plus2hsensgprobe",
  "GSE48611"= "hgu133plus2hsensgprobe",
  "GSE1789" = "hgu133ahsensgprobe",
  "GSE11877" = "hgu133plus2hsensgprobe",
  "GSE16176" = "hgu133plus2hsensgprobe",
  "GSE16677" = "hgu133plus2hsensgprobe",
  "GSE17459" = "hgu133plus2hsensgprobe", #c("GPL570", "GPL96"),
  "GSE19680" = "hgu133plus2hsensgprobe",
  "GSE19681" = "hgu133plus2hsensgprobe",
  "GSE20910" = "hgu133plus2hsensgprobe",
  "GSE30517"= "hgu133ahsensgprobe"
)

#Populated by InstallArrayPackages
platforms_to_package_list = list(
  "GPL570" = "hgu133plus2hsensgprobe",
  "GPL96" = "hgu133ahsensgprobe"
) 

#----------functions-------------

create_unique_vector = function(platforms_list){
  # returns a vector of one key for every unique value in the list
  
  platforms_vector = unlist(platforms_list)
  value_to_keys = split(names(platforms_vector), platforms_vector)
  keys_vector = unname(unlist(lapply(value_to_keys, `[`, 1)))
  return(keys_vector)
}