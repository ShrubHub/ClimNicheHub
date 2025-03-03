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
table_species_name[1:10,]
table_species_name[n_occ< 100,]

table_species_name_no_dup <- table_species_name[!duplicated(GBIF_name),]

row_to_dl <- 701:993
row_to_dl[44]
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
## Lophozia serpens n_occ = 0  Persicaria laxmannii n_occ = 1 
table_species_name_no_dup[644,]

to_hist <- table_species_name$n_occ
to_hist <- ifelse(to_hist>50000,50000,to_hist)

hist(to_hist,nc=100)
abline(v=c(1000,10000),col="red")

plan(sequential)
occ_to_test2 <- occ( "Salix arctica",limit=1000,date = c('1900-01-01', '2025-02-01'))

readRDS(file.path("raw_occ","raw_occ_Agrostis_vinealis.RData"))
readRDS(file.path("raw_occ","raw_occ_Antennaria_alpina.RData"))
readRDS(file.path("raw_occ","raw_occ_Achillea_atrata.RData"))


occ("Anemone parviflora",limit = 1,has_coords = T,date = c('1900-01-01', '2025-02-01'))
occ("Achillea atrata",limit = 1,has_coords = T,date = c('1900-01-01', '2025-02-01'))



t1 <- Sys.time()
test_occ <- occ("Salix arctica",limit = 30000)
test_occ <- occ("Eriophorum vaginatum",limit = 30000)
t2 <- Sys.time()
t2-t1

test_occ_df <- occ2df(test_occ)

test_occ_df <- data.table(test_occ_df)

test_occ_sf <- st_as_sf(test_occ_df[!is.na(longitude),],coords=c("longitude","latitude"),crs=st_crs(4326))
mapview(test_occ_sf)
library(rgbif)
test_occ_2 <- occ_search(scientificName = "Salix arctica",limit =30000 )
test_occ_2 <- occ_search(scientificName = "Eriophorum vaginatum",limit =30000 )
test_occ_2$data

library(occTest)

occ_to_test$gbif$meta

occ(ids = "4900876566")
occ_to_test2 <- get_gbif(sp_name = "Salix arctica",occ_samp=100000,grain = 10000)
occ_to_test <- get_gbif(sp_name = "Chamaenerion angustifolium",occ_samp=1000)
occ_to_test <- get_gbif(sp_name = "Epilobium angustifolium L.",search = F ,occ_samp=1000)
occ_to_test3 <- get_gbif(sp_name = "Equisetum hyemale",search = F ,occ_samp=100000)

table_species_name

occ("Delphinus delphis", limit = 6000)

occ_to_test2 <- occ( "Salix arctica",limit=2000)

eco.terra <- read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)

range_test <- get_range(occ_to_test2,eco.terra,"WWF_MHTNAM",res = 10,buffer_width_point = 1,buffer_width_polygon = 2)
mapview(range_test$rangeOutput) + mapview(test_occ_sf_3_sf)

range_nex_ext <- project(range_test$rangeOutput,test)

test_occ_sf_3_sf <- st_as_sf(occ_to_test2,coords = c("decimalLongitude","decimalLatitude"),crs = st_crs (4326))

occ_to_test <- occ("Viola altaica", start = 100,limit = 776,has_coords = T)  # date = c('2010-01-01', '2015-12-31')
occ_to_test  <- occ2df(occ_to_test)
settings <- readRDS (system.file('ext/exSettings.rds',package = 'occTest'))

settings$tableSettings$x.field <- "longitude"
settings$tableSettings$y.field <- "latitude"

settings$tableSettings$x.field <- "decimalLongitude"
settings$tableSettings$y.field <- "decimalLatitude"

settings$analysisSettings$geoOutliers$doGeoOutliers <- F
settings$analysisSettings$envOutliers$doEnvOutliers <- F
settings$analysisSettings$countryStatusRange$doRangeAnalysis <- F
settings$analysisSettings$filterAtlas <- F



test_salix <- occTest(sp.name ="Salix arctica Pall." ,sp.table =  test_occ_df,r.env = test,tableSettings = settings$tableSettings)
test_salix <- occTest(sp.name ="Salix arctica Pall." ,
                      sp.table =  occ_to_test,
                      r.env = test,
                      tableSettings = settings$tableSettings,
                      analysisSettings = settings$analysisSettings,
                      doParallel = T,
                      mc.cores = 10)
test_salix$centroidDetection_CoordinateCleaner_test
test_salix_filter <- occFilter(test_salix,errorAcceptance = "strict")
plot(test_salix,test_salix_filter)[[1]]
plot(test_salix,test_salix_filter)[[2]]
plot(test_salix,test_salix_filter)[[3]]

#### taxonmy matching between ITEX and GBIF names ####

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
