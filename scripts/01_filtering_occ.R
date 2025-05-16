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

table_species_name <- fread(file.path("species_list","GBIF_ITEX_query_list_05_2025.csv"))
table_species_name <- table_species_name[order(accepted_name),]

## function to get filtered occ file path
get_raw_file_path <- function(sp_name){
  sp_name <- str_remove_all(sp_name,"[.]")
  name_file <- paste0("filtered_occ_",str_replace_all(sp_name," ","_"),".RData")
  name_file <- file.path("filtered_occ",name_file)
  
  return(name_file)
}

table_species_name[,filtered_occ_path := get_raw_file_path(accepted_name)]

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
shp_coastline <- st_buffer(shp_coastline,dist=0.01) ## buffer a more detailed coastline shapefile, otherwise some 
## coastal occurrences are filtered but should be kept

settings$analysisSettings$landSurfacePol<- shp_coastline

#### loading the occurrences ####
library(lubridate)
date <- ymd_hms( file.info(path_occ)$ctime)
date< ymd_hms( "2025-05-10 10:37:40 UTC" )
path_occ[date> ymd_hms( "2025-05-10 10:37:40 UTC" )]
path_occ <- list.files("raw_occ",full.names = T,pattern = ".RData")
path_occ <- path_occ[date> ymd_hms( "2025-05-10 10:37:40 UTC" )]

#### occurrence testing and filtering ####

foreach(path = path_occ)%do%{
  occ_tmp <- readRDS(path)
  
  sp_name <- occ_tmp$name_itex[1]
  
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

## done ! ##

