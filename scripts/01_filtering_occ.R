#### Filtering occurrences ####

# This script uses the occTest package to perform simple filtering operations.
# The occurrences are thinned so that there is only one occurrence per CHELSA grid,
# duplicated occurrences are also removed, as well as occurrences corresponding to country
# centroids.

#### packages ####

library(data.table)
library(occTest)
library(stringr)
library(foreach)
library(terra)
library(sf)

#### meta data ####

table_species_name <- readRDS(file.path("species_list","ITEX_species_names.RData"))
table_species_name <- fread(file.path("species_list","Taxonomy_ITEX_gbif_ClimNicheHub_2025.csv"))

chelsa_resolution <- rast(file.path("chelsa_data_V2.1","CHELSA_bio10_1981_2010_V.2.1.tif"))

#### occTest settings ####

settings <- readRDS (system.file('ext/exSettings.rds',package = 'occTest'))

settings$tableSettings$x.field <- "longitude"
settings$tableSettings$y.field <- "latitude"

settings$analysisSettings$geoOutliers$doGeoOutliers <- F
settings$analysisSettings$envOutliers$doEnvOutliers <- F
settings$analysisSettings$countryStatusRange$doRangeAnalysis <- F
settings$analysisSettings$filterAtlas <- F

library(rnaturalearth)
shp_coastline <- ne_download(10)
sf_use_s2(FALSE)
shp_coastline <- st_buffer(shp_coastline,dist=0.01)

settings$analysisSettings$landSurfacePol<- shp_coastline


mapview::mapview(shp_coastline)
mapview::mapview(shp_coastline[shp_coastline$ADMIN=="Canada",])

#### loading the occurrences ####

path_occ <- list.files("raw_occ",full.names = T,pattern = ".RData")

path_occ <- path_occ[c(177,665,647,648,649:670)]


#### occurrence testing and filtering ####
readRDS(path_occ[177])
readRDS(path_occ[665])
readRDS("raw_occ/raw_occ_Siversia_pentapetala.RData")
path_occ[17]

248

foreach(path = path_occ[600:991])%do%{
  occ_tmp <- readRDS(path)
  
  sp_name <- occ_tmp$name[1]
  
  occTest_tmp <- occTest(sp.name = sp_name,
                         sp.table =  occ_tmp,
                         r.env = chelsa_resolution,
                         tableSettings = settings$tableSettings,
                         analysisSettings = settings$analysisSettings)
  
  occFilter_tmp <- occFilter(occTest_tmp,errorAcceptance = "strict")
  
  final_occ <- data.table(occFilter_tmp$filteredDataset)
  
  name_file <- paste0("filtered_occ_",str_replace_all(sp_name," ","_"),".RData")
  
  print(paste0( " -|-|- ",sp_name," -|-|- ","done ✔"))
  
  saveRDS(final_occ,file.path("filtered_occ",name_file))
  
  
}

result_occ_path <-  list.files("filtered_occ",full.names = T,pattern = ".RData")

readRDS(path_occ[663])

occTest_tmp <- occTest(sp.name = occ_tmp$name[1],
                      sp.table =  occ_tmp,
                      r.env = chelsa_resolution,
                      tableSettings = settings$tableSettings,
                      analysisSettings = settings$analysisSettings)

plot(occTest_tmp)

occTest_tmp_sf <- st_as_sf(occTest_tmp,coords = c("longitude","latitude"),crs=st_crs(4326))

mapview::mapview(occTest_tmp_sf,zcol = "Reason")+mapview::mapview(shp_coastline)

occFilter_tmp <- occFilter(occTest_tmp,errorAcceptance = "strict")

final_occ <- occFilter_tmp$filteredDataset


plot(occTest_tmp,occFilter_tmp)[[2]]
occTest_tmp$Reason
occTest_tmp$

#### occurrence filtering ####


