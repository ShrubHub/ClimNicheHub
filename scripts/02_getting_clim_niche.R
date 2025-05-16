#### computing climatic niche ####

# This script uses the filtered occurrence and overlay them on the rasters of climatologies
# obtained in CHELSA V2.1. 

#### packages ####

library(data.table)
library(terra)
library(sf)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(ggpubr)

#### meta data and climate data ####

table_species_name <- fread(file.path("species_list","GBIF_ITEX_query_list_05_2025.csv"))

table_species_name[,raw_occ_path := get_raw_file_path(accepted_name)]
table_species_name[,filtered_occ_path := get_raw_file_path(accepted_name)]

## clim data
all_chelsa_path <- list.files("chelsa_data_V2.1",full.names = T,pattern = ".tif")
all_chelsa_path <- all_chelsa_path[c(1:19,21:33,20)] ## put cell_id raster at the end
chelsa_full <- rast(all_chelsa_path)

names(chelsa_full) <- str_remove(names(chelsa_full),"_1981_2010_V.2.1")
names(chelsa_full) <- str_remove(names(chelsa_full),"CHELSA_")

# run once after downloading Chelsa data 
# cell_id <- chelsa_full$bio1
# names(cell_id) <- "cell_id"
# values(cell_id) <- 1:ncell(cell_id)
# writeRaster(cell_id,file.path("chelsa_data_V2.1","CHELSA_cell_id_1981_2010_V.2.1.tif"))
#### loading occurrences ####

result_occ_path <-  list.files("filtered_occ",full.names = T,pattern = ".RData")

list_of_occ <- lapply(result_occ_path,readRDS)

occ_dt <- rbindlist(list_of_occ,fill = T);rm(list_of_occ);gc()

occ_dt <- occ_dt[,c(1:6,147)]

occ_dt[,accepted_name := ifelse(is.na(name_itex),name,name_itex)]

## some names aren't stored in the same column because of some early version of the dowload script
occ_dt[,name := accepted_name ]
occ_dt[,name_itex := NULL ]
occ_dt[,accepted_name := NULL ]

final_count <- occ_dt[,.N,by= name]

enough_occ <- final_count[N>= 100,]
## selecting species with at least 100 occ
occ_dt <- occ_dt[name%in%enough_occ$name,]
## spatial object
occ_dt_sf <- st_as_sf(occ_dt,coords=c("longitude","latitude"),crs=st_crs(4326))

bugged_sp <- table_species_name_count[is.na(N) ,accepted_name]
bugged_sp
##no occurrences found for "Lophozia serpens" 

#### extraction of climate data ####

all_clim <- extract(chelsa_full,vect(occ_dt_sf))
all_clim <- data.table(all_clim)
occ_dt_clim <- cbind(occ_dt,all_clim[,-"ID"])

rm(occ_dt);rm(occ_dt_sf);rm(all_clim);gc()

#### summary statistics computation ####

final_count_all <- occ_dt_clim[,.N,by= name]
final_count_boreal_tundra <- occ_dt_clim[kg5<=7,.N,by= name]
final_count_tundra <- occ_dt_clim[kg5<=4,.(`N_tund` =.N),by= name]

final_count_all <- merge(final_count_all,final_count_boreal_tundra,all.x=T,by="name",suffixes = c("_all","_bor_tund"))
final_count_all <- merge(final_count_all,final_count_tundra,all.x=T,by="name")
# final_count_all[,N_bor_tund := ifelse(is.na(N_bor_tund),0,N_bor_tund)]
# final_count_all[,N_tund := ifelse(is.na(N_tund),0,N_tund)]

sampling_summary <- merge(final_count,final_count_all[,-"N_all"],all.x =T)
colnames(sampling_summary)[2] <- "N_all"
sampling_summary[,included_in_all := as.numeric(N_all>=100)]
sampling_summary[,included_in_bor_tund := as.numeric(N_bor_tund>=100)]
sampling_summary[,included_in_bor_tund := ifelse(is.na(included_in_bor_tund),0,included_in_bor_tund)]

write.table(sampling_summary,file.path("ClimNiche_database","sampling_summary.csv"),row.names = F,sep = ",")

#### ClimNiche computation####

col_to_compute <- c(names(chelsa_full)[-33],"longitude","latitude")
col_to_compute <- col_to_compute[c(1,12:19,2:11,20:30,32:34)]  ### ordering, not interested in kg5

## run this for the complete computation
occ_dt_clim_full <- occ_dt_clim
occ_dt_clim <- occ_dt_clim_full

## run this for the bor_tund computation
occ_dt_clim <- occ_dt_clim[kg5<=7,]
occ_dt_clim <- occ_dt_clim[name %in% final_count_all[N_bor_tund>=100,name],]
## precision of CHELSA
round_digit <- c(2,1,3,1,2,2,1,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,2,0,3,3)

## mean computation and rounding to the nearest CHELSA decimal pre
occ_dt_clim[ , metric := "mean" ]
mean_clim_dt <- occ_dt_clim[ , lapply(.SD,mean,na.rm=T), by = .(name,metric),.SDcols = col_to_compute]
mean_clim_dt <- mean_clim_dt[ , mapply(round,.SD,round_digit,SIMPLIFY = F), by = .(name,metric),.SDcols = col_to_compute] ##rounding to the decimal precision of CHELSA

## median computation
occ_dt_clim[ , metric := "median" ]
med_clim_dt <- occ_dt_clim[ , lapply(.SD,median,na.rm=T), by = .(name,metric),.SDcols = col_to_compute]

## Quantile 5 and 95 computation
occ_dt_clim[ , metric := "Q" ]
Q_clim_dt <- occ_dt_clim[ , lapply(.SD,quantile,probs=c(0.05,0.95),na.rm=T), by = .(name,metric),.SDcols = col_to_compute]
Q_clim_dt[,metric := rep(c("Q05","Q95"),times= nrow(Q_clim_dt)/2)]

## interquantile range computation
Q_clim_dt_tmp <- Q_clim_dt
Q_clim_dt_tmp$metric <- "range"
range_clim_dt <- Q_clim_dt_tmp[,lapply(.SD,diff),by = .(name,metric),.SDcols = col_to_compute];rm(Q_clim_dt_tmp)

## optimum of the density computation 
get_density_max <- function(x,na.rm=T){
  tmp <- try(density(x,na.rm=na.rm,n = 512*2),silent = T)
  
  if(class(tmp)=="density") return(tmp$x[which.max(tmp$y)]) else return(as.numeric(NA))
  
}
occ_dt_clim[ , metric := "optimum" ]
opt_clim_dt <- occ_dt_clim[ , lapply(.SD,get_density_max,na.rm=T), by = .(name,metric),.SDcols = col_to_compute]

## rounding
opt_clim_dt <- opt_clim_dt[ , mapply(round,.SD,round_digit,SIMPLIFY = F), by = .(name,metric),.SDcols = col_to_compute] ##rounding to the decimal precision of CHELSA

## getting the database together
ClimNicheHub <- rbind(mean_clim_dt,med_clim_dt,Q_clim_dt,range_clim_dt,opt_clim_dt)
ClimNicheHub <- ClimNicheHub[order(name),]

ClimNicheHub[,print(summary(.SD)),by = metric]

#### Export ####

write.table(ClimNicheHub,file.path("ClimNiche_database","ClimNiche_climate_summary.csv"),row.names = F,sep = ",")

write.table(ClimNicheHub,file.path("ClimNiche_database","ClimNiche_boreal_tundra_climate_summary.csv"),row.names = F,sep = ",")


#### Export of the whole chelsa sampling 
grid_clim <- occ_dt_clim[!duplicated(cell_id),-c("name","prov","date","key","metric")]


saveRDS(grid_clim,file.path("Complete_sampling","all_CHELSA_cells.RData"))
saveRDS(occ_dt_clim[,c("name","cell_id")],file.path("Complete_sampling","species_cells.RData"))

saveRDS(grid_clim,file.path("Complete_sampling","boreal_tundra_CHELSA_cells.RData"))
saveRDS(occ_dt_clim[,c("name","cell_id")],file.path("Complete_sampling","boreal_tundra_species_cells.RData"))

#### PCA-based approach ####

grid_clim <- occ_dt_clim[!duplicated(cell_id),-c("name","prov","date","key","metric")]
set.seed(0)
sample_of_grid <- sample(grid_clim$cell_id,20000)

plot(grid_clim[cell_id%in% sample_of_grid,.(bio1,bio4)])

library(ade4)

apply(grid_clim,2,function(x)sum(is.na(x)))

to_pca <- grid_clim[cell_id%in%sample_of_grid,-c("longitude" , "latitude","kg5")]
rownames(to_pca) <- to_pca$cell_id
to_pca[,cell_id:=NULL]

to_pca <- to_pca[,-c("gddlgd0","gdgfgd0","fgd","gsp","gsl")]
to_pca <- to_pca[complete.cases(to_pca),]
pca_clim_grid <- dudi.pca(to_pca,nf = 4, scannf = F)
inertia.dudi(pca_clim_grid)
s.arrow(pca_clim_grid$co)

col_to_get <- colnames(to_pca)

test_salix <- suprow(pca_clim_grid,occ_dt_clim[name == "Betula pubescens",..col_to_get])

projected_species <- suprow(pca_clim_grid,occ_dt_clim[,..col_to_get])
projected_species <- data.table(projected_species$lisup)
projected_species[,name:=occ_dt_clim$name]

pca_species <- projected_species[,.(Axis1 = mean(Axis1,na.rm = T),Axis2 = mean(Axis2,na.rm = T),
                                    Axis3 = mean(Axis3,na.rm = T),Axis4 = mean(Axis4,na.rm = T)),by=name]

pca_ClimNicheHub <- pca_species

pca_species <- merge(pca_species,ClimNicheHub[metric == "median",])

ggplot(pca_species,aes(x = Axis3, y = Axis4, color = bio10,size = bio4))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_c()

plot(test_salix$lisup)
write.table(pca_ClimNicheHub,file.path("ClimNiche_database","ClimNiche_all_pca_summary.csv"),row.names = F,sep = ",")

write.table(pca_ClimNicheHub,file.path("ClimNiche_database","ClimNiche_boreal_tundra_pca_summary.csv"),row.names = F,sep = ",")

#### Plotting ####
country_shape <- ne_download(50)

occ_dt_clim <- occ_dt_clim[kg5<=4,]
occ_dt_clim <- occ_dt_clim_full

plot_one_sp <- function(name_sp, what = "bio10",labs = NULL,occ = occ_dt_clim){
  set.seed(0)
  if(is.null(labs))labs <- what
  
  occ_dt_clim_one_sp <- occ[name ==  name_sp,]
  occ_dt_clim_one_sp_subset <- occ_dt_clim_one_sp
  if(nrow(occ_dt_clim_one_sp)>1000) occ_dt_clim_one_sp_subset <- occ_dt_clim_one_sp[sample(1:nrow(occ_dt_clim_one_sp),1000,F),]
  
  sf_occ <- st_as_sf(occ_dt_clim_one_sp_subset,coords=c("longitude" , "latitude"),crs=st_crs(4326))
  
  bbox <- st_bbox(sf_occ)
  
  print(bbox)
  plot_1 <- ggplot()+
    geom_sf(data = country_shape)+
    geom_sf(data = sf_occ,aes(fill= get(what)),shape=21,color="grey10",alpha = 0.95)+
    scale_fill_gradient(low = "white", high = "orange")+
    coord_sf(xlim = bbox[c(1,3)],ylim = bbox[c(2,4)])+
    theme_bw()+
    labs(fill = labs,title = name_sp)+
    theme(legend.position = "bottom",legend.key.width = unit( 20,"mm"))
 
  plot_2 <- ggplot(occ_dt_clim_one_sp,aes(x = get(what)))+
    theme_bw()+
    labs(x = labs,y="")+
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+
    geom_density(fill="orange",alpha = 0.5)
  
  final_plot <- ggarrange(plot_1,plot_2,
                          ncol=1,align = "v",heights = c(1.3,1),
                          labels = c(" "," "))
  
  return(final_plot)
     
}

save_a_plot <- function(name_sp, what = "bio10",labs = NULL){
  
  tmp <- plot_one_sp(name_sp,what,labs)
  
  name_sp_2 <- str_replace_all(name_sp," ","_")
  
  file_name <- paste0(name_sp_2,"_",what,".jpg")
    
  ggsave(file.path("figures","individual_species_map",file_name),
         plot = tmp,
         width = 180,
         height = 140,
         unit = "mm",
         dpi = 150)
  
}

save_a_plot("Eriophorum vaginatum","bio10","TWQ (°C)")
save_a_plot("Salix chamissonis","bio10","TWQ (°C)")
save_a_plot("Salix arctica","bio10","TWQ (°C)")
save_a_plot("Anemone narcissiflora","bio10","TWQ (°C)")
save_a_plot("Epilobium palustre","bio10","TWQ (°C)")
save_a_plot("Achillea atrata","bio10","TWQ (°C)")
save_a_plot("Betula pubescens","bio10","TWQ (°C)")
save_a_plot("Betula nana","bio10","TWQ (°C)")
save_a_plot("Androsace obtusifolia","bio10","TWQ (°C)")
save_a_plot("Androsace chamaejasme","bio10","TWQ (°C)")

plot_one_sp("Achillea millefolium","bio10","TWQ (°C)",occ_dt_clim_full)

plot_one_sp("Salix arctica","bio10","TWQ (°C)")
plot_one_sp("Betula pubescens","bio10","TWQ (°C)",occ_dt_clim)
plot_one_sp("Betula pubescens","bio10","TWQ (°C)",occ_dt_clim[kg5<=7])
plot_one_sp("Betula pubescens","bio10","TWQ (°C)",occ_dt_clim[kg5<=4])

plot_one_sp("Andromeda polifolia","bio10","TWQ (°C)",occ_dt_clim)
plot_one_sp("Andromeda polifolia","bio10","TWQ (°C)",occ_dt_clim[kg5<=7])
plot_one_sp("Andromeda polifolia","bio10","TWQ (°C)",occ_dt_clim[kg5<=4])
Salix petrophila


plot_one_sp("Silene acaulis","bio18","TWQ (°C)",occ_dt_clim)

