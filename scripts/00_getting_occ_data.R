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
table_species_name <- readRDS(file.path("species_list","ITEX_species_names.RData"))
table_species_name <- fread(file.path("species_list","Taxonomy_ITEX_gbif_ClimNicheHub_2025.csv"))

table_species_name_no_dup <- table_species_name[!duplicated(GBIF_name),]

row_to_dl <- 1:993

t1 <- Sys.time()
foreach(gbif_name = table_species_name_no_dup[row_to_dl,GBIF_name],
        itex_name = table_species_name_no_dup[row_to_dl,ITEX_name])%do%{
  occ_dl <- occ(gbif_name,limit=10000,has_coords = T,date = c('1900-01-01', '2025-02-01'))
  occ_dl <- data.table(occ2df(occ_dl))
  occ_dl[,name_gbif := name]
  occ_dl[,name := itex_name]
  tmp_name <- str_remove_all(itex_name,"[.]")
  name_file <- paste0("raw_occ_",str_replace_all(tmp_name," ","_"),".RData")
  
  saveRDS(occ_dl,file.path("raw_occ",name_file))
  
}

t2 <- Sys.time()
t2-t1

to_hist <- table_species_name$n_occ
to_hist <- ifelse(to_hist>50000,50000,to_hist)

hist(to_hist,nc=100)
abline(v=c(1000,10000),col="red")


#### taxonmy matching between ITEX and GBIF names, no need to run again ####

table_species_name <- fread(file.path("species_list","Taxonomy_ITEX_ClimNicheHub_2025.csv"))

table_species_name[SPECIES_NAME == "Ranunculus acris var. monticola",SPECIES_NAME:= "Ranunculus paishanensis"]
table_species_name[SPECIES_NAME == "Rumex aquaticus subsp. arcticus",SPECIES_NAME:= "Rumex arcticus"]

table_species_name[,ITEX_name := SPECIES_NAME]
table_species_name[,GBIF_name := SPECIES_NAME]

table_species_name[ITEX_name == "Ledum palustre","GBIF_name"] <- "Ledum palustre L."
table_species_name[ITEX_name == "Polygonum viviparum","GBIF_name"] <- "Polygonum viviparum L."
table_species_name[ITEX_name == "Salix fuscescens","GBIF_name"] <- "Salix fuscescens Andersson"
table_species_name[ITEX_name == "Trientalis europaea","GBIF_name"] <- "Trientalis europaea L."
table_species_name[ITEX_name == "Gnaphalium supinum","GBIF_name"] <- "Gnaphalium supinum L."
table_species_name[ITEX_name == "Saxifraga nivalis","GBIF_name"] <- "Saxifraga nivalis L."
table_species_name[ITEX_name == "Trientalis europaea","GBIF_name"] <- "Trientalis europaea L."
table_species_name[ITEX_name == "Saxifraga hieraciifolia","GBIF_name"] <- "Saxifraga hieraciifolia Waldst. & Kit."
table_species_name[ITEX_name == "Alopecurus alpinus","GBIF_name"] <- "Alopecurus alpinus Vill."
table_species_name[ITEX_name == "Epilobium angustifolium","GBIF_name"] <- "Epilobium angustifolium L."

table_species_name <- table_species_name[order(ITEX_name),]


how_many <- foreach(i = unique(table_species_name$GBIF_name ),.combine = rbind)%do%{
  test_one_sp <- occ(i,limit = 1,has_coords = T)
  res <- data.table(SPECIES_NAME = i,
                    gbif_name = test_one_sp$gbif$meta$opts$scientificName  ,
                    how_many = test_one_sp$gbif$meta$found)
  cat(paste0(i, " // "))
  res


}

how_many

dt_species[SPECIES_NAME %in% how_many[how_many == 0,SPECIES_NAME]]

how_many[order(-how_many),]

how_many[,n_occ:=how_many]

how_many[how_many == 0 ,]

table_species_name <- merge(table_species_name,how_many[,c("gbif_name","n_occ")],by.x ="GBIF_name", by.y ="gbif_name" )

### some species will not be downloaded with just the binomial name, need the authors also 
# library(taxize)
not_found_itex_name <- how_many[how_many == 0  ,SPECIES_NAME ]
# not_found <- get_gbifid(not_found_itex_name)
gbif_name <- c("Anemone narcissiflora L.","Anemone obtusiloba D.Don","Arctophila fulva (Trin.) Andersson","Bryum teres Warnst.","Sarmentypnum tundrae", "Craspedia aurantia var. jamesii",
  "Deyeuxia angustifolia Vickery","Epilobium latifolium L.","Pappochroma nitidum","Euphrasia picta","Hypnum hamulosum Schimp.","Hypnum revolutum (Mitt.) Lindb.",
  "Lophozia serpens","Scleranthus biflorus", "Parasenecio auriculata", "Saxifraga nelsoniana D.Don",
  "Saxifraga reflexa Hook.","Scleranthus singuliflorus")

foreach(itex_name = not_found_itex_name, gbif_name = gbif_name) %do%{
  table_species_name[ITEX_name == itex_name,"GBIF_name"] <- gbif_name
  
}



write.table(table_species_name,file.path("species_list","Taxonomy_ITEX_gbif_ClimNicheHub_2025.csv"),row.names = F)
# 
# itex.all <- data.table(itex.all)
# itex.all.2 <- itex.all[!grepl("XXX",Name),]
# nrow(itex.all.2)
# nrow(itex.all.2[Name %in% table_species_name[n_occ >= 500,ITEX_name],])/nrow(itex.all.2)
# 
# saveRDS(table_species_name,"ITEX_species_names.RData")
# table_species_name <- readRDS("ITEX_species_names.RData")
