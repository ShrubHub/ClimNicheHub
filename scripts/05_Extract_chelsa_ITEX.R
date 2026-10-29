library(sf)
library(terra)
library(data.table)
library(stringr)

ITEX_place <- fread("ITEX_METADATA_2026(post_review)_12March2026.csv")
ITEX_place[SITE == "LATNJA" & SUBSITE == "WET_MEADOW", LAT :=  LAT/10e6]
ITEX_place[SITE == "LATNJA" & SUBSITE == "WET_MEADOW", LONG :=   LONG /10e6]

ITEX_place_sf <- st_as_sf(ITEX_place[complete.cases(ITEX_place[,c("LONG","LAT")]),],coords = c("LONG","LAT"),crs = st_crs(4326))


all_chelsa_path <- list.files("chelsa_data_V2.1",full.names = T,pattern = ".tif")
all_chelsa_path <- all_chelsa_path[c(1:19,21:33,20)] ## put cell_id raster at the end
chelsa_full <- rast(all_chelsa_path)

names(chelsa_full) <- str_remove(names(chelsa_full),"_1981_2010_V.2.1")
#names(chelsa_full) <- str_remove(names(chelsa_full),"CHELSA_")


clim_itex <- data.table(extract(chelsa_full,vect(ITEX_place_sf)))
clim_itex <- cbind(st_drop_geometry(ITEX_place_sf)[,c("SITE","SUBSITE")],clim_itex[,-c("ID","cell_id")])

write.table(clim_itex,"ITEX_climate.csv",row.names = F,sep = ",")

cavm_legend <- fread("Raster CAVM legend.csv")

CAVM <- rast("raster_cavm_v1.tif")
CAVM <- classify(CAVM,matrix(c(91,92,93,99,NA,NA,NA,NA),ncol = 2))
CAVM <- project(CAVM,chelsa_full)
plot(CAVM)

chelsa_full_masked <- mask(chelsa_full,CAVM)

chelsa_full_masked <- trim(chelsa_full_masked)

chelsa_full_masked[[c("CHELSA_bio10","CHELSA_bio18")]]

plot(chelsa_full_masked$CHELSA_bio10)

sampled_raster_every <- spatSample(chelsa_full_masked,2000,na.rm = T)

write.table(sampled_raster_every,"CHELSA_sample_CAVM.csv",sep = ",",row.names = F)
