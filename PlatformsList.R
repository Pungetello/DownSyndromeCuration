platforms_list = list(
  
  #RNA datasets
  "GSE101942" = NA, #downloadable
  "GSE109293" = NA,
  "GSE109294" = NA,
  "GSE184771" = NA, 
  "GSE190053" = NA, #downloadable
  "GSE202938" = NA, 
  "GSE210117" = NA, #big, finish prefetch for later
  
  #Affymetrix datasets
  "GSE110064" = "hgu133plus2hsentrezgprobe", #change names
  "GSE48611"= "hgu133plus2hsentrezgprobe",
  "GSE1789" = "hgu133ahsentrezgprobe",
  # #"GSE11877" = "hgu133plus2hsensgprobe", #the big one
  "GSE16176" = "hgu133plus2hsentrezgprobe"
  # "GSE16677" = "hgu133plus2hsensgprobe",
  # #"GSE17459" = "hgu133plus2hsensgprobe", #c("GPL570", "GPL96"), #does indeed need both, if we use this we'll need to redo stuff
  # "GSE19680" = "hgu133plus2hsensgprobe",
  # "GSE19681" = "hgu133plus2hsensgprobe",
  # "GSE20910" = "hgu133plus2hsensgprobe",
  # "GSE30517"= "hgu133ahsensgprobe",
  # "GSE138861" = "clariomshumanhsensgprobe", #"[Clariom_S_Human] Affymetrix Clariom S Assay, Human (Includes Pico Assay) (GPL23159)",
  # "GSE143885" = "clariomshumanhsensgprobe", #"[Clariom_S_Human] Affymetrix Clariom S Assay, Human (Includes Pico Assay) (GPL23159)",
  # #"GSE149464" = "mogene10stmmensgprobe",
  # #"GSE149465" = "mogene10stmmensgprobe"
  # "GSE16676" = "mouse4302mmensgprobe"
)

#Populated by InstallArrayPackages
platforms_to_package_list = list(
  "GPL570" = "hgu133plus2hsensgprobe",
  "GPL96" = "hgu133ahsensgprobe"
) 

platforms_to_refrence_genomes = list(
  "GSE109293" = "C57BL/6J.129P2",
  "GSE109294" = "C57BL/6J.129P2",
  "GSE184771" = "C57BL/6.129/SvJ", 
  "GSE202938" = "Ts65Dn, TcMAC21"
  #"GSE210117" = "DBA/2J x B6EiC3Sn/J" #(Ts65Dn),
  
)

#----------functions-------------

create_unique_vector = function(platforms_list){
  # returns a vector of one key for every unique value in the list
  
  platforms_vector = unlist(platforms_list)
  value_to_keys = split(names(platforms_vector), platforms_vector)
  keys_vector = unname(unlist(lapply(value_to_keys, `[`, 1)))
  return(keys_vector)
}