library(spocc)
library(data.table)
library(sf)
library(mapview)

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

readRDS("raw_occ/raw_occ_Aconitum_tschangbaischanense.RData")

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

library(climenv)
library(stringr)
chelsa("chelsa_data",var = "all" )

wget_dist <- fread("envidatS3paths.txt",header = F)
path_url <- wget_dist$V1
path_local <- str_extract(path_url,"CHELSA.*")
path_local <- file.path("chelsa_data",path_local)
path_local <- str_replace(path_local,"-","_")
download.file(path_url[c(1,2)],path_local[c(1,2)],mode="wb")
download.file(path_url[c(3:32)],path_local[c(3:32)],mode="wb")
command <- paste0("download.file(path_url[",3:32,"],path_local[",3:32,"],mode='wb')")
parse(text  = command[1])

for (i in 1:length(command))eval(parse(text = command[i]))

options(timeout=1600)
library(foreach)
foreach(url = path_url, write = path_local,.errorhandling = "remove")%do%{
  
  download.file(url,path_local,mode="wb")
  
}

library(terra)
test <- rast(path_local[2])
test <- rast(file.path("chelsa_data/CHELSA_bio10_1981_2010_V.2.1.tif"))

all_raster <- list.files("chelsa_data",full.names = T)
all_ras <- rast(all_raster)

cell_id <- all_ras$CHELSA_bio1_1981_2010_V.2.1
names(cell_id) <- "cell_id"
values(cell_id) <- 1:ncell(cell_id)

all_ras <- c(all_ras,cell_id)

test_occ_sf_2 <- data.table(test_salix_filter$filteredDataset)
test_occ_sf_2_sf <- st_as_sf(test_occ_sf_2,coords=c("longitude","latitude"),crs=st_crs(4326))


range_test <- mask(test,range_nex_ext)
global(range_test,mean,na.rm=T)
global(range_test,quantile,na.rm=T,probs= c(0.05,0.95))

extract_clim <- extract(all_ras,vect(test_occ_sf_2_sf))
extract_clim <- data.table(extract_clim)
extract_clim <- extract_clim[!duplicated(cell_id),]


hist(extract_clim$CHELSA_bio10_1981_2010_V.2.1,breaks = 20)
median(extract_clim$CHELSA_bio10_1981_2010_V.2.1)
mean(extract_clim$CHELSA_bio10_1981_2010_V.2.1)
quantile(extract_clim$CHELSA_bio10_1981_2010_V.2.1,probs= c(0.05,0.95))

set.seed(0)
id_sampled <- sample(1:nrow(extract_clim),2000,replace = F)
extract_clim_sample <- extract_clim[id_sampled,]
plot(test)
library(lubridate)

test_occ_sf_2_sf$date <- ymd(test_occ_sf_2_sf$date)
test_occ_sf_2_sf <- test_occ_sf_2_sf[!is.na(test_occ_sf_2_sf$date),]
test_occ_sf_2_sf$date_num <- as.numeric(test_occ_sf_2_sf$date)

mapview(test_occ_sf_2_sf,zcol="date_num")


#### taxonmy matching between ITEX and GBIF names ####

load("data/ITEX_ALL.RData")
species_all <- unique(itex.all$Name)
species_all <- species_all[!grepl("XXX",species_all)]

table_species_name <- data.table(ITEX_name =species_all  , GBIF_name = species_all)

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

saveRDS(table_species_name,"ITEX_species_names.RData")
table_species_name <- readRDS("ITEX_species_names.RData")

how_many <- foreach(i = table_species_name$GBIF_name ,.combine = rbind)%do%{
  test_one_sp <- occ(i,limit = 1,has_coords = T)
  res <- data.table(SPECIES_NAME = i, 
                    gbif_name = test_one_sp$gbif$meta$opts$scientificName  ,
                    how_many = test_one_sp$gbif$meta$found)
  cat(paste0(i, " // "))
  res
  
  
}

hist(how_many$how_many,nc=40)

summary(how_many)
how_many[how_many == 0,SPECIES_NAME ]

dt_species[SPECIES_NAME %in% how_many[how_many == 0,SPECIES_NAME]]

how_many[order(-how_many),]

library(taxize)
not_found_itex_name <- how_many[how_many == 0,SPECIES_NAME ]
not_found <- get_gbifid(not_found_itex_name)
gbif_name <- c("Anemone narcissiflora L.","Anemone obtusiloba D.Don","Arabis petraea (L.) Lam","Arctophila fulva (Trin.) Andersson","Arenaria lychnidea M.Bieb.",
  "Arnica alpina (L.)","Arnica stricta A.Nelson","Chamaeorchis alpinus (L.) Rich.","Cladina mitis (Sandst.) Hustich","Comastoma ferratum",
  "Epilobium latifolium L.","Gagea fistulosa (Ramond ex DC.) Ker Gawl.","Gentiana algida Pall.","Leucorchis albida (L.) E.Mey.","Orthothecium chryseon (Schwägr.) Schimp.",
  "Potentilla gelida C.A.Mey.","Riccardia latifrons (Lindb.) Lindb.","Saxifraga caespitosa L.","Saxifraga nelsoniana D.Don",
  "Saxifraga reflexa Hook.","Sieversia pentapetala (L.) Greene")

foreach(itex_name = not_found_itex_name, gbif_name = gbif_name) %do%{
  table_species_name[ITEX_name == itex_name,"GBIF_name"] <- gbif_name
  
}

how_many_full <- foreach(i = table_species_name$GBIF_name ,.combine = rbind)%do%{
  test_one_sp <- occ(i,limit = 1,has_coords = T)
  res <- data.table(SPECIES_NAME = i, 
                    gbif_name = test_one_sp$gbif$meta$opts$scientificName  ,
                    how_many = test_one_sp$gbif$meta$found)
  cat(paste0(i, " // "))
  res
  
  
}

how_many_full[how_many == 0,]

table_species_name[,n_occ := how_many_full$how_many]

table_species_name[n_occ >10000,]
table_species_name[n_occ >5000,]

itex.all <- data.table(itex.all)
itex.all.2 <- itex.all[!grepl("XXX",Name),]
nrow(itex.all.2)
nrow(itex.all.2[Name %in% table_species_name[n_occ >= 500,ITEX_name],])/nrow(itex.all.2)

saveRDS(table_species_name,"ITEX_species_names.RData")
table_species_name <- readRDS("ITEX_species_names.RData")
99000

table_species_name[n_occ>99000,]
