#### downloading occurrences data ####

# This scripts pulls the 10,000 latests occurrences from GBIF starting from 2025-02-01 
# so that later dowload pulls the same occurrences.
# The selected are every species indentified at least once to the species level in the full
# International Tundra EXperiment (ITEX), and store them in RData format in the raw_occ folder

#### packages ####

library(spocc)
library(data.table)
library(occTest)
library(stringr)
library(foreach)
library(doFuture)

#### occurrence download ####

# table_species_name <- fread(file.path("species_list","Taxonomy_ITEX_gbif_ClimNicheHub_2025.csv"))
table_species_name <- fread(file.path("species_list","GBIF_ITEX_query_list_05_2025.csv"))
table_species_name <- table_species_name[order(accepted_name),]

## function to get raw_occ file path
get_raw_file_path <- function(sp_name){
  sp_name <- str_remove_all(sp_name,"[.]")
name_file <- paste0("raw_occ_",str_replace_all(sp_name," ","_"),".RData")
name_file <- file.path("raw_occ",name_file)

return(name_file)
}


table_species_name[,raw_occ_path := get_raw_file_path(accepted_name)]
table_species_name[!duplicated(accepted_name),]

t1 <- Sys.time()

foreach(name_to_dl = table_species_name[,unique(accepted_name)] ) %do%{ # raw_occ_path %in% toget
  
          names_to_query <- table_species_name[accepted_name == name_to_dl,unique(species)]
          print(names_to_query)
          occ_dl <- occ(names_to_query,limit=10000,has_coords = T,date = c('1900-01-01', '2025-02-01'))
          occ_dl <- data.table(occ2df(occ_dl))
          occ_dl[,name_gbif := name_to_dl]
          occ_dl[,name_itex := name_to_dl]
          
          saveRDS(occ_dl,get_raw_file_path(name_to_dl))
          
        }


t2 <- Sys.time()
t2-t1

## done ! ##

to_hist <- table_species_name$n_occ
to_hist <- ifelse(to_hist>50000,50000,to_hist)

hist(to_hist,nc=100)
abline(v=c(1000,10000),col="red")

# this taxonomy matching is deprecated as ITEX and GBIF backbone are merged
# 
# # table_species_name <- readRDS(file.path("species_list","ITEX_species_names.RData"))
# old_table_species_name <- fread(file.path("species_list","Taxonomy_ITEX_gbif_ClimNicheHub_2025.csv"))
# old_table_species_name <- old_table_species_name[!duplicated(GBIF_name),]
# # #### taxonmy matching between ITEX and GBIF names, no need to run again ####
# 
# old_table_species_name <- fread(file.path("species_list","Taxonomy_ITEX_ClimNicheHub_2025.csv"))
# old_table_species_name
# table_species_name <- fread(file.path("species_list","ITEX_GBIF_backbone_taxonomy.csv"))
# table_species_name <- fread(file.path("species_list","GBIF_ITEX_query_list_05_2025.csv"))
# 
# table_species_name[,GBIF_name :=  species ]
# table_species_name <- table_species_name[full_species_name == T,]
# 
# table_species_name[SPECIES_NAME == "Ranunculus acris var. monticola",SPECIES_NAME:= "Ranunculus paishanensis"]
# table_species_name[SPECIES_NAME == "Rumex aquaticus subsp. arcticus",SPECIES_NAME:= "Rumex arcticus"]
# 
# table_species_name[,ITEX_name := SPECIES_NAME]
# table_species_name[,GBIF_name := SPECIES_NAME]
# 
# table_species_name[ITEX_name == "Ledum palustre","GBIF_name"] <- "Ledum palustre L."
# table_species_name[ITEX_name == "Polygonum viviparum","GBIF_name"] <- "Polygonum viviparum L."
# table_species_name[ITEX_name == "Salix fuscescens","GBIF_name"] <- "Salix fuscescens Andersson"
# table_species_name[ITEX_name == "Trientalis europaea","GBIF_name"] <- "Trientalis europaea L."
# table_species_name[ITEX_name == "Gnaphalium supinum","GBIF_name"] <- "Gnaphalium supinum L."
# table_species_name[ITEX_name == "Saxifraga nivalis","GBIF_name"] <- "Saxifraga nivalis L."
# table_species_name[ITEX_name == "Trientalis europaea","GBIF_name"] <- "Trientalis europaea L."
# table_species_name[ITEX_name == "Saxifraga hieraciifolia","GBIF_name"] <- "Saxifraga hieraciifolia Waldst. & Kit."
# table_species_name[ITEX_name == "Alopecurus alpinus","GBIF_name"] <- "Alopecurus alpinus Vill."
# table_species_name[ITEX_name == "Epilobium angustifolium","GBIF_name"] <- "Epilobium angustifolium L."
# 
# table_species_name <- table_species_name[order(ITEX_name),]
# 
# 
# how_many <- foreach(i = unique(table_species_name$GBIF_name ),.combine = rbind)%do%{
#   test_one_sp <- occ(i,limit = 1,has_coords = T)
#   res <- data.table(SPECIES_NAME = i,
#                     gbif_name = test_one_sp$gbif$meta$opts$scientificName  ,
#                     how_many = test_one_sp$gbif$meta$found)
#   cat(paste0(i, " // "))
#   res
# 
# 
# }
# 
# how_many <- foreach(i = unique(table_species_name$accepted_name),.combine = rbind)%do%{
#   
#   sp_to_get <- table_species_name[accepted_name == i, unique(species)]
#   
#   test_one_sp <- occ(sp_to_get,limit = 1,has_coords = T)
#   res <- data.table(SPECIES_NAME = i,
#                     #gbif_name = test_one_sp$gbif$meta$opts$scientificName  ,
#                     how_many = test_one_sp$gbif$meta$found)
#   cat(paste0(i, " // "))
#   res
# 
# 
# }
# 
# 
# 
# how_many
# 
# dt_species[SPECIES_NAME %in% how_many[how_many == 0,SPECIES_NAME]]
# 
# occ("Alopecurus alpinus Sm.",limit = 1)
# 
# how_many[order(-how_many),]
# 
# how_many[,n_occ:=how_many]
# 
# how_many[how_many == 0 ,]
# 
# table_species_name <- merge(table_species_name,how_many[,c("gbif_name","n_occ")],by.x ="GBIF_name", by.y ="gbif_name" )
# 
# ### some species will not be downloaded with just the binomial name, need the authors also 
# # library(taxize)
# not_found_itex_name <- how_many[how_many == 0  ,SPECIES_NAME ]
# # not_found <- get_gbifid(not_found_itex_name)
# gbif_name <- c("Anemone narcissiflora L.","Anemone obtusiloba D.Don","Arctophila fulva (Trin.) Andersson","Bryum teres Warnst.","Sarmentypnum tundrae", "Craspedia aurantia var. jamesii",
#   "Deyeuxia angustifolia Vickery","Epilobium latifolium L.","Pappochroma nitidum","Euphrasia picta","Hypnum hamulosum Schimp.","Hypnum revolutum (Mitt.) Lindb.",
#   "Lophozia serpens","Scleranthus biflorus", "Parasenecio auriculata", "Saxifraga nelsoniana D.Don",
#   "Saxifraga reflexa Hook.","Scleranthus singuliflorus")
# 
# foreach(itex_name = not_found_itex_name, gbif_name = gbif_name) %do%{
#   table_species_name[ITEX_name == itex_name,"GBIF_name"] <- gbif_name
#   
# }
# 
# 
# 
# write.table(table_species_name,file.path("species_list","Taxonomy_ITEX_gbif_ClimNicheHub_2025.csv"),row.names = F)
# # 
# # itex.all <- data.table(itex.all)
# # itex.all.2 <- itex.all[!grepl("XXX",Name),]
# # nrow(itex.all.2)
# # nrow(itex.all.2[Name %in% table_species_name[n_occ >= 500,ITEX_name],])/nrow(itex.all.2)
# # 
# # saveRDS(table_species_name,"ITEX_species_names.RData")
# # table_species_name <- readRDS("ITEX_species_names.RData")
# 
# raw_sp <- list.files("raw_occ",full.names = T)
# raw_sp_theory <- unique(table_species_name$raw_occ_path)
# 
# 
# file.copy(sp_to_remove_file_from,sp_to_remove_file_to)
# file.remove(sp_to_remove_file_from)
# 
# 
# raw_sp[!raw_sp%in%raw_sp_theory]
# raw_sp_theory[!raw_sp_theory %in%raw_sp]
# 
# toget <- raw_sp_theory[!raw_sp_theory %in%raw_sp]
# 
# ## dl
# 
# sp_to_remove <- final_count[! name %in% table_species_name$accepted_name,]
# sp_to_remove <- old_table_species_name[! GBIF_name %in% table_species_name$accepted_name,]
# 
# sp_to_remove_file <- sp_to_remove
# 
# sp_to_remove_file <- str_remove_all(sp_to_remove$GBIF_name,"[.]")
# sp_to_remove_file <- paste0("raw_occ_",str_replace_all(sp_to_remove_file," ","_"),".RData")
# list.files("raw_occ",full.names = T)
# 
# sp_to_remove_file_from <- file.path("raw_occ",sp_to_remove_file)
# sp_to_remove_file_to<- file.path("raw_occ_2",sp_to_remove_file)
# 
# file.copy(sp_to_remove_file_from,sp_to_remove_file_to)
# file.remove(sp_to_remove_file_from)
# 
# 
# sp_already_computed <- table_species_name[!duplicated(accepted_name),][accepted_name %in% old_table_species_name$GBIF_name,]
# 
# sp_already_computed_old <- table_species_name[!duplicated(accepted_name),][accepted_name %in% old_table_species_name$GBIF_name,]
# 
# sp_to_compute <- table_species_name[!duplicated(accepted_name),][! accepted_name %in% old_table_species_name$GBIF_name,]
# 
